"""assemble stage: fetch the selected sequencing runs, assemble them, and write
the genome contract (``genomes/``, ``selection.tsv``, the manifest).

Each run in ``reads.tsv`` is fetched from the FASTQ locations ENA reported
(a local path is copied instead), verified against its checksum, assembled
with the adapter that accepts its platform and layout, filtered to contigs of
a minimum length, and labelled with the family, genus and species the reads
stage resolved. Runs that cannot be fetched or assembled are written to
``excused_runs.tsv`` so the completeness guard of later stages does not count
them as missing. Finished assemblies carry a marker under ``assemblies/<run>/``
and are skipped on a re-run.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..assemblers.base import AssembleParams as AdapterParams
from ..assemblers.base import ReadSet, accepting_assemblers, registry, select_assembler
from ..assemblers.contigs import ContigStats, filter_contigs
from ..classifiers.base import Classification, ClassifyParams
from ..classifiers.base import registry as classifier_registry
from ..core import http
from ..core.context import WorkdirContext
from ..core.contracts import (
    ASSEMBLY_STATS_TSV,
    EXCUSED_RUNS_TSV,
    READS_TSV,
    SELECTION_TSV,
    AssemblyStatsRow,
    ExcusedRun,
    ReadRow,
    SelectionRow,
    accession_from_filename,
    atomic_replace,
    genome_filename,
    read_reads,
    sanitise_taxon_tokens,
    write_assembly_stats,
    write_excused_runs,
    write_selection,
)
from ..core.errors import (
    MissingBinaryError,
    RepGenRError,
    ToolExecutionError,
    UserInputError,
    WorkdirError,
)
from ..core.executors import parallel_map
from ..core.manifest import record_from_selection
from ..core.process import check_free_disk, link_or_copy, remove_tree, staged_dir
from ..polishers.base import PolishParams, accepting_polishers, select_polisher
from ..polishers.base import registry as polisher_registry
from .assemble_qc import checkm2_db_from_env, preflight_checkm2, run_checkm2
from .ingest import OUTGROUP_ACCESSION_TXT

GTDB_SKETCH_ENV = "REPGENR_GTDB_SKETCH"
GTDB_LINEAGES_ENV = "REPGENR_GTDB_LINEAGES"
_DONE_MARKER = "assembly.ok"
_CONTIGS_NAME = "contigs.fasta"
# Rough peak disk per run: the FASTQ files, their uncompressed form and the
# assembler's scratch.
_DISK_FACTOR = 4


@dataclass
class AssembleParams:
    assembler: str = "auto"
    threads: int = 16
    # Concurrent runs; memory, not CPU, bounds this, so it stays small. None
    # resolves to 2, or to 1 when a long-read run is pending.
    jobs: int | None = None
    memory_gb: int = 16
    min_contig_length: int = 500
    # Correct long-read assemblies with the run's reads: auto picks by platform
    # (medaka for ONT, racon for PacBio CLR, none for HiFi and Illumina).
    polisher: str = "auto"
    polish_rounds: int = 1
    # A FASTA file to set aside as the outgroup (ingest semantics).
    outgroup: str | None = None
    # Add the assemblies to a working directory that already holds a selection
    # (metadata + genome, or ingest) instead of replacing it: existing rows,
    # files and outgroup stay, rows of the same accession are replaced.
    append: bool = False
    keep_reads: bool = False
    keep_files: bool = False
    # Quality: CheckM2 runs when a database is configured (flag or CHECKM2DB);
    # assemblies outside the gate are excused rather than selected.
    checkm2_db: str | None = None
    min_completeness: float = 50.0
    max_contamination: float = 10.0
    # Classification: verifies the submitted organism against a GTDB sketch
    # (flag or REPGENR_GTDB_SKETCH / REPGENR_GTDB_LINEAGES). "auto" runs the
    # sourmash classifier when a sketch is configured, "none" never.
    classifier: str = "auto"
    gtdb_sketch: str | None = None
    gtdb_lineages: str | None = None
    # Tool tuning from ``--tool-arg``; each adapter declares the keys it reads.
    extra: dict = field(default_factory=dict)


@dataclass
class _Outcome:
    row: ReadRow
    assembler: str | None = None
    polisher: str | None = None
    polish_rounds: int = 0
    stats: ContigStats | None = None
    tool_stats: dict = field(default_factory=dict)
    excused: ExcusedRun | None = None
    quality: tuple[float, float] | None = None
    gtdb: Classification | None = None
    # Resolved filename tokens and where they came from.
    label: tuple[str, str, str] | None = None
    label_source: str = "metadata"
    taxonomy_flag: str = ""
    # Tool versions a reused marker recorded (its run needs no preflight).
    versions: dict = field(default_factory=dict)


def run(ctx: WorkdirContext, params: AssembleParams) -> int:
    logger = ctx.logger
    reads_path = ctx.workdir / READS_TSV
    if not reads_path.exists():
        raise WorkdirError(f"{READS_TSV} not found in {ctx.workdir}. Run the reads stage first.")
    rows = read_reads(reads_path)
    if not rows:
        raise WorkdirError(f"{READS_TSV} lists no runs.")

    assemblies = ctx.workdir / "assemblies"
    scratch = ctx.scratch_dir / "assemble"
    assemblies.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)

    plan = _plan(rows, params, assemblies, logger=logger)
    if params.assembler == "auto":
        _excuse_missing_assemblers(plan, logger)
        _require_something_to_assemble(plan)
    if params.polisher == "auto":
        _warn_missing_polishers(plan, logger)
    versions = {k: v for o in plan for k, v in o.versions.items()}
    versions.update(_preflight(plan, logger))
    pending = [o for o in plan if o.excused is None and o.stats is None]
    check_free_disk(
        ctx.workdir,
        sum(sum(o.row.fastq_bytes) for o in pending) * _DISK_FACTOR,
        logger,
        what=f"assemble {len(pending)} sequencing runs",
    )
    requested = params.jobs if params.jobs is not None else _default_jobs(pending)
    jobs = max(1, min(requested, len(pending) or 1))
    threads_each = max(1, params.threads // jobs)
    logger.info(
        "Assembling %d runs (%d already done, %d excused) with %d concurrent jobs, %d threads each",
        len(pending),
        sum(1 for o in plan if o.stats is not None),
        sum(1 for o in plan if o.excused is not None),
        jobs,
        threads_each,
    )

    def work(outcome: _Outcome) -> _Outcome:
        return _fetch_and_assemble(
            outcome,
            params,
            threads_each,
            assemblies / outcome.row.run_accession,
            scratch / outcome.row.run_accession,
            versions,
            logger,
        )

    done = parallel_map(work, pending, jobs, logger=logger)
    by_run = {o.row.run_accession: o for o in [*plan, *done]}
    outcomes = [by_run[r.run_accession] for r in rows]

    assembled = [o for o in outcomes if o.stats is not None and o.stats.n_contigs > 0]
    for o in assembled:
        o.label = _taxa(o.row)

    # Every run's contigs file is named alike; QC and classification see them
    # through links named by run accession, which are unique.
    named = named_links(
        {o.row.run_accession: assemblies / o.row.run_accession / _CONTIGS_NAME for o in assembled},
        scratch / "named",
    )
    checkm2_db = params.checkm2_db or checkm2_db_from_env()
    gtdb_sketch = params.gtdb_sketch or os.environ.get(GTDB_SKETCH_ENV)
    gtdb_lineages = params.gtdb_lineages or os.environ.get(GTDB_LINEAGES_ENV)
    classifier_name = classifier_for(params.classifier, gtdb_sketch)
    quality, classified = assess(
        named,
        checkm2_db=checkm2_db,
        classifier=classifier_name,
        gtdb_sketch=gtdb_sketch,
        gtdb_lineages=gtdb_lineages,
        threads=params.threads,
        extra=params.extra,
        scratch=scratch,
        versions=versions,
        logger=logger,
    )
    if quality is not None:
        apply_quality(assembled, quality, params.min_completeness, params.max_contamination, logger)
        assembled = [o for o in assembled if o.excused is None]
    n_disagree = 0
    if classified is not None:
        n_disagree = apply_classification(assembled, classified, versions, logger)

    excused = [o.excused for o in outcomes if o.excused is not None]
    excused_path = ctx.workdir / EXCUSED_RUNS_TSV
    if excused:
        write_excused_runs(excused_path, excused)
    else:
        excused_path.unlink(missing_ok=True)
    if not assembled:
        raise WorkdirError(
            f"None of the {len(rows)} runs produced an accepted assembly; see {EXCUSED_RUNS_TSV}."
        )

    new_rows = [_selection_row(o) for o in assembled]
    if params.append:
        selection_rows = _append_rows(ctx, new_rows, params, logger)
        ctx.genomes_dir.mkdir(parents=True, exist_ok=True)
        for o in assembled:
            _unlink_previous(ctx.genomes_dir, o.row.run_accession)
            link_or_copy(
                assemblies / o.row.run_accession / _CONTIGS_NAME, ctx.genomes_dir / _name(o)
            )
        ctx.manifest.upsert_many([record_from_selection(r, "sra") for r in new_rows])
    else:
        selection_rows = list(new_rows)
        outgroup_row = _stage_outgroup(ctx, params.outgroup, logger)
        if outgroup_row is not None:
            selection_rows.append(outgroup_row)
        with staged_dir(ctx.genomes_dir) as genomes_dir:
            for o in assembled:
                link_or_copy(
                    assemblies / o.row.run_accession / _CONTIGS_NAME, genomes_dir / _name(o)
                )
        ctx.manifest.replace_genomes([record_from_selection(r, "sra") for r in selection_rows])
    write_selection(ctx.workdir / SELECTION_TSV, selection_rows)
    write_assembly_stats(ctx.workdir / ASSEMBLY_STATS_TSV, [_stats_row(o) for o in assembled])
    outgroup_row = next((r for r in selection_rows if r.is_outgroup), None)

    assemblers_used = sorted({o.assembler for o in assembled if o.assembler})
    ctx.config.record_stage(
        "assemble",
        # The assemblers that built the accepted genomes; params keep the request.
        tool=",".join(assemblers_used) or params.assembler,
        params={
            **asdict(params),
            # The databases as resolved (flag or environment variable).
            "checkm2_db": checkm2_db,
            "gtdb_sketch": gtdb_sketch,
            "gtdb_lineages": gtdb_lineages,
            "classifier_effective": classifier_name,
            "assemblers_used": assemblers_used,
            "polishers_used": sorted({o.polisher for o in assembled if o.polisher}),
            "n_assembled": len(assembled),
            "n_excused": len(excused),
            "n_disagree": n_disagree,
            "outgroup_accession": outgroup_row.accession if outgroup_row else None,
        },
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info(
        "Assembled %d of %d runs (%d excused, %d classifier disagreements); wrote %s",
        len(assembled),
        len(rows),
        len(excused),
        n_disagree,
        SELECTION_TSV,
    )
    return len(assembled)


_LONG_READ_PLATFORMS = frozenset({"OXFORD_NANOPORE", "PACBIO_SMRT"})


def _why(exc: RepGenRError) -> str:
    """The failure for excused_runs.tsv: a tool failure keeps its output tail."""
    if isinstance(exc, ToolExecutionError) and exc.output:
        return f"{exc}: {exc.output}"
    return str(exc)


def _default_jobs(pending: list[_Outcome]) -> int:
    """Two short-read assemblies fit side by side; a long-read one wants the machine."""
    return 1 if any(o.row.platform in _LONG_READ_PLATFORMS for o in pending) else 2


# --- planning -------------------------------------------------------------------


def _plan(
    rows: list[ReadRow],
    params: AssembleParams,
    assemblies: Path,
    *,
    check_settings: bool = True,
    logger: logging.Logger | None = None,
) -> list[_Outcome]:
    """Decide, per run, whether it is done, excused up front, or to be assembled.

    A finished run is reused when its marker's settings agree with ``params``
    (``check_settings``; ``reads-gather`` only collects results and passes
    False). A higher contig floor is applied to the finished contigs; any
    other difference assembles the run again.
    """
    if params.assembler != "auto" and params.assembler not in registry.names():
        raise UserInputError(
            f"Unknown assembler {params.assembler!r}; available: {', '.join(registry.names())}."
        )
    if params.polisher not in ("auto", "none") and params.polisher not in polisher_registry.names():
        raise UserInputError(
            f"Unknown polisher {params.polisher!r}; available: "
            f"{', '.join(polisher_registry.names())}, none."
        )
    plan = []
    for row in rows:
        outcome = _Outcome(row=row)
        if _reuse_finished(outcome, params, assemblies / row.run_accession, check_settings, logger):
            pass
        elif not row.fastq_urls:
            outcome.excused = ExcusedRun(row.run_accession, "fetch", "no_fastq_mirror")
        else:
            reads = ReadSet(
                row.run_accession, row.platform, row.instrument_model, row.layout, (), row.bases
            )
            if params.assembler == "auto":
                outcome.assembler = select_assembler(registry, reads)
            elif registry.create(params.assembler).accepts(reads):
                outcome.assembler = params.assembler
            if outcome.assembler is None:
                outcome.excused = ExcusedRun(row.run_accession, "assemble", "unsupported_platform")
            elif params.polisher == "auto":
                outcome.polisher = select_polisher(polisher_registry, reads)
            elif params.polisher != "none":
                probe_cls = polisher_registry.get(params.polisher)
                if probe_cls.__new__(probe_cls).accepts(reads):
                    outcome.polisher = params.polisher
        plan.append(outcome)
    return plan


# --- reuse of finished runs ---------------------------------------------------------


def _read_marker(path: Path) -> dict | None:
    """A finished run's marker, or None when absent or unreadable (cut short by a kill)."""
    try:
        done = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(done, dict) or "assembler" not in done or "stats" not in done:
        return None
    return done


def _settings(params: AssembleParams) -> dict:
    """The result-affecting settings a marker records, as requested."""
    return {
        "assembler": params.assembler,
        "polisher": params.polisher,
        "polish_rounds": params.polish_rounds,
        "min_contig_length": params.min_contig_length,
        "extra": dict(sorted(params.extra.items())),
    }


def _accepted_extras(name: str | None, reg) -> frozenset[str]:  # noqa: ANN001
    if not name or name not in reg.names():
        return frozenset()
    return reg.get(name).capabilities.accepted_extras


def _polisher_request(name: str, row: ReadRow) -> str:
    """A polisher request as it applies to one run: 'none' when nothing would polish it."""
    reads = ReadSet(row.run_accession, row.platform, row.instrument_model, row.layout, (), 0)
    if name == "none" or not accepting_polishers(polisher_registry, reads):
        return "none"
    if name != "auto" and name in polisher_registry.names():
        cls = polisher_registry.get(name)
        if not cls.__new__(cls).accepts(reads):
            return "none"
    return name


def _setting_changes(done: dict, params: AssembleParams, row: ReadRow) -> list[str]:
    """How ``params`` differ from the settings a finished run was built with.

    A request that names the tool the marker records agrees with it, so
    ``--assembler skesa`` reuses a run that ``auto`` assembled with skesa.
    """
    old = done["settings"]
    changes = []
    if old["assembler"] != params.assembler and params.assembler != done["assembler"]:
        changes.append(f"assembler {old['assembler']} -> {params.assembler}")
    old_pol = _polisher_request(old["polisher"], row)
    new_pol = _polisher_request(params.polisher, row)
    if old_pol != new_pol and new_pol != (done.get("polisher") or "none"):
        changes.append(f"polisher {old['polisher']} -> {params.polisher}")
    elif done.get("polisher") and old["polish_rounds"] != params.polish_rounds:
        changes.append(f"polish_rounds {old['polish_rounds']} -> {params.polish_rounds}")
    keys = _accepted_extras(done["assembler"], registry) | _accepted_extras(
        done.get("polisher"), polisher_registry
    )
    old_extra = {k: str(v) for k, v in old["extra"].items() if k in keys}
    new_extra = {k: str(v) for k, v in params.extra.items() if k in keys}
    if old_extra != new_extra:
        changes.append(f"--tool-arg {old_extra} -> {new_extra}")
    return changes


def _reuse_finished(
    outcome: _Outcome,
    params: AssembleParams,
    run_dir: Path,
    check_settings: bool,
    logger: logging.Logger | None,
) -> bool:
    """Fill ``outcome`` from a finished run in ``run_dir``; False when it must be assembled."""
    contigs = run_dir / _CONTIGS_NAME
    done = _read_marker(run_dir / _DONE_MARKER)
    if done is None or not contigs.exists():
        return False
    row = outcome.row
    old_floor = None
    # Markers written before the settings were recorded are reused as they are.
    if check_settings and "settings" in done:
        changes = _setting_changes(done, params, row)
        old_floor = done["settings"]["min_contig_length"]
        if old_floor > params.min_contig_length:
            changes.append(f"min_contig_length {old_floor} -> {params.min_contig_length}")
        if changes:
            if logger is not None:
                logger.info(
                    "%s: finished with other settings (%s); assembling again",
                    row.run_accession,
                    "; ".join(changes),
                )
            return False
    outcome.assembler = done["assembler"]
    outcome.polisher = done.get("polisher") or None
    outcome.polish_rounds = done.get("polish_rounds", 0)
    outcome.stats = ContigStats(**done["stats"])
    outcome.tool_stats = done.get("tool_stats", {})
    outcome.versions = dict(done.get("tool_versions") or {})
    if not outcome.versions and done.get("version"):
        outcome.versions = {done["assembler"]: done["version"]}
    if old_floor is not None and old_floor < params.min_contig_length:
        _refilter(outcome, done, run_dir, params.min_contig_length, logger)
    return True


def _refilter(
    outcome: _Outcome,
    done: dict,
    run_dir: Path,
    min_length: int,
    logger: logging.Logger | None,
) -> None:
    """Apply a higher contig floor to finished contigs, which hold every contig
    above the old floor in assembly order, so the result equals filtering the
    raw assembly. A floor that keeps nothing excuses the run and leaves its
    files as they are, for a later, lower floor."""
    run = outcome.row.run_accession
    contigs = run_dir / _CONTIGS_NAME
    trial = run_dir / f"{_CONTIGS_NAME}.refilter"
    stats = filter_contigs(contigs, trial, min_length=min_length, prefix=run)
    if stats.n_contigs == 0:
        trial.unlink(missing_ok=True)
        outcome.stats = None
        outcome.excused = ExcusedRun(
            run, "assemble", f"assembly_failed: no contig of {min_length} bp or more"
        )
        return
    trial.replace(contigs)
    done["stats"] = asdict(stats)
    done["settings"]["min_contig_length"] = min_length
    _write_marker(run_dir / _DONE_MARKER, done)
    outcome.stats = stats
    if logger is not None:
        logger.info("%s: contigs filtered again at %d bp", run, min_length)


def _write_marker(path: Path, marker: dict) -> None:
    with atomic_replace(path) as fo:
        fo.write(json.dumps(marker, indent=1))


def _marker_versions(versions: dict[str, str], names: list[tuple[str, object]]) -> dict[str, str]:
    """The entries of ``versions`` that belong to the named adapters (by tool and binary)."""
    keys: set[str] = set()
    for name, reg in names:
        caps = reg.get(name).capabilities  # type: ignore[attr-defined]
        keys.add(caps.name)
        keys.update(b.name for b in caps.required_binaries)
    return {k: v for k, v in versions.items() if k in keys}


ASSEMBLER_NOT_INSTALLED = "assembler_not_installed"


def _excuse_missing_assemblers(plan: list[_Outcome], logger: logging.Logger) -> None:
    """Re-excuse runs an adapter accepts but whose tool is not installed.

    Under ``auto``, ``_plan`` excuses every run without an available adapter
    as ``unsupported_platform``. That is right when no adapter takes the
    platform; when one would, the tool is missing, so the run is excused as
    ``assembler_not_installed`` and one warning per platform names the
    adapters. The planner itself is unchanged because ``reads-gather`` calls
    it on hosts without assemblers.
    """
    needed: dict[str, list[str]] = {}
    for o in plan:
        if o.excused is None or o.excused.reason != "unsupported_platform":
            continue
        row = o.row
        reads = ReadSet(
            row.run_accession, row.platform, row.instrument_model, row.layout, (), row.bases
        )
        names = accepting_assemblers(registry, reads)
        if names:
            o.excused = ExcusedRun(row.run_accession, "assemble", ASSEMBLER_NOT_INSTALLED)
            needed.setdefault(row.platform, names)
    for platform, names in sorted(needed.items()):
        n = sum(
            1
            for o in plan
            if o.excused is not None
            and o.excused.reason == ASSEMBLER_NOT_INSTALLED
            and o.row.platform == platform
        )
        logger.warning(
            "%d %s run(s) excused as %s: none of %s is installed. Put one on PATH, run "
            "with --container, or narrow the runs with reads --platform, then rerun "
            "assemble with --force.",
            n,
            platform,
            ASSEMBLER_NOT_INSTALLED,
            ", ".join(names),
        )


def _warn_missing_polishers(plan: list[_Outcome], logger: logging.Logger) -> None:
    """Warn when ``--polisher auto`` leaves runs unpolished for want of a tool.

    A run to be assembled with no polisher chosen is either one no adapter
    takes (Illumina, HiFi) or one an adapter would take whose tool is not
    installed. The latter is still assembled, unpolished, and one warning per
    platform names the adapters, as for a missing assembler.
    """
    needed: dict[str, list[str]] = {}
    counts: dict[str, int] = {}
    for o in plan:
        if o.excused is not None or o.stats is not None or o.polisher is not None:
            continue
        row = o.row
        reads = ReadSet(
            row.run_accession, row.platform, row.instrument_model, row.layout, (), row.bases
        )
        names = accepting_polishers(polisher_registry, reads)
        if names:
            needed.setdefault(row.platform, names)
            counts[row.platform] = counts.get(row.platform, 0) + 1
    for platform, names in sorted(needed.items()):
        logger.warning(
            "%d %s run(s) will be assembled without polishing: none of %s is installed. "
            "Put one on PATH or run with --container to polish them, or pass "
            "--polisher none to assemble unpolished without this warning.",
            counts[platform],
            platform,
            ", ".join(names),
        )


def _require_something_to_assemble(plan: list[_Outcome]) -> None:
    """Exit 4 when nothing can be assembled and a missing assembler is the cause."""
    missing = [
        o for o in plan if o.excused is not None and o.excused.reason == ASSEMBLER_NOT_INSTALLED
    ]
    if not missing or any(o.excused is None for o in plan):
        return
    by_platform: dict[str, list[str]] = {}
    for o in missing:
        reads = ReadSet(
            o.row.run_accession, o.row.platform, o.row.instrument_model, o.row.layout, (), 0
        )
        by_platform.setdefault(o.row.platform, accepting_assemblers(registry, reads))
    detail = "; ".join(
        f"{platform} runs need one of {', '.join(names)}"
        for platform, names in sorted(by_platform.items())
    )
    raise MissingBinaryError(
        f"No run can be assembled: no assembler is installed for --assembler auto ({detail}). "
        "Put one on PATH or run with --container."
    )


def _preflight(plan: list[_Outcome], logger: logging.Logger) -> dict[str, str]:
    """Check the assemblers the pending runs need; finished runs need none."""
    versions: dict[str, str] = {}
    pending = [o for o in plan if o.assembler and o.excused is None and o.stats is None]
    for name in sorted({o.assembler for o in pending if o.assembler}):
        versions.update(registry.create(name).preflight())
    for name in sorted({o.polisher for o in pending if o.polisher}):
        versions.update(polisher_registry.create(name).preflight())
    return versions


# --- one run --------------------------------------------------------------------


def _fetch_and_assemble(
    outcome: _Outcome,
    params: AssembleParams,
    threads: int,
    out_dir: Path,
    run_scratch: Path,
    versions: dict[str, str],
    logger: logging.Logger,
) -> _Outcome:
    """Fetch and assemble one run into ``out_dir`` (contigs and the done marker)."""
    row = outcome.row
    # The marker names finished contigs; a run assembled again has none until it ends.
    (out_dir / _DONE_MARKER).unlink(missing_ok=True)
    if run_scratch.exists():
        remove_tree(run_scratch)
    run_scratch.mkdir(parents=True)
    try:
        files = _fetch(row, run_scratch, logger)
    except RepGenRError as exc:
        logger.warning("%s: fetch failed (%s)", row.run_accession, exc)
        outcome.excused = ExcusedRun(row.run_accession, "fetch", f"download_failed: {exc}")
        remove_tree(run_scratch)
        return outcome

    assert outcome.assembler is not None
    adapter = registry.create(outcome.assembler)
    reads = ReadSet(
        row.run_accession, row.platform, row.instrument_model, row.layout, files, row.bases
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        result = adapter.assemble(
            reads,
            run_scratch / "asm",
            AdapterParams(threads=threads, memory_gb=params.memory_gb, extra=dict(params.extra)),
            logger,
        )
    except RepGenRError as exc:
        logger.warning("%s: assembly failed (%s)", row.run_accession, exc)
        outcome.excused = ExcusedRun(row.run_accession, "assemble", f"assembly_failed: {_why(exc)}")
        if not params.keep_files:
            remove_tree(run_scratch)
        return outcome
    contigs_in = result.contigs
    polish_stats: dict = {}
    if outcome.polisher is not None:
        try:
            polished = polisher_registry.create(outcome.polisher).polish(
                reads,
                result.contigs,
                run_scratch / "polish",
                PolishParams(
                    threads=threads, rounds=params.polish_rounds, extra=dict(params.extra)
                ),
                logger,
            )
        except RepGenRError as exc:
            logger.warning("%s: polishing failed (%s)", row.run_accession, exc)
            outcome.excused = ExcusedRun(
                row.run_accession, "assemble", f"polish_failed: {_why(exc)}"
            )
            if not params.keep_files:
                remove_tree(run_scratch)
            return outcome
        contigs_in = polished.contigs
        outcome.polish_rounds = polished.rounds
        polish_stats = dict(polished.tool_stats)
    try:
        stats = filter_contigs(
            contigs_in,
            out_dir / _CONTIGS_NAME,
            min_length=params.min_contig_length,
            prefix=row.run_accession,
        )
    except RepGenRError as exc:
        logger.warning("%s: contig filtering failed (%s)", row.run_accession, exc)
        outcome.excused = ExcusedRun(row.run_accession, "assemble", f"assembly_failed: {_why(exc)}")
        if not params.keep_files:
            remove_tree(run_scratch)
        return outcome
    if stats.n_contigs == 0:
        outcome.excused = ExcusedRun(
            row.run_accession,
            "assemble",
            f"assembly_failed: no contig of {params.min_contig_length} bp or more",
        )
        (out_dir / _CONTIGS_NAME).unlink(missing_ok=True)
        if not params.keep_files:
            remove_tree(run_scratch)
        return outcome

    outcome.stats = stats
    outcome.tool_stats = dict(result.tool_stats)
    tools: list[tuple[str, object]] = [(outcome.assembler, registry)]
    if outcome.polisher is not None:
        tools.append((outcome.polisher, polisher_registry))
    marker = {
        "assembler": outcome.assembler,
        "version": versions.get(adapter.capabilities.name, ""),
        "polisher": outcome.polisher or "",
        "polish_rounds": outcome.polish_rounds,
        "polish_stats": polish_stats,
        "stats": asdict(stats),
        "tool_stats": outcome.tool_stats,
        "tool_versions": _marker_versions(versions, tools),
        "settings": _settings(params),
    }
    _write_marker(out_dir / _DONE_MARKER, marker)
    if params.keep_files:
        pass
    elif params.keep_reads:
        remove_tree(run_scratch / "asm")
    else:
        remove_tree(run_scratch)
    logger.info(
        "%s: %d contigs, %d bp, N50 %d (%s)",
        row.run_accession,
        stats.n_contigs,
        stats.total_length,
        stats.n50,
        outcome.assembler,
    )
    return outcome


def _fetch(row: ReadRow, run_scratch: Path, logger: logging.Logger) -> tuple[Path, ...]:
    """Bring the run's FASTQ files into scratch, verified when a checksum is known."""
    files = []
    md5s = list(row.fastq_md5) + [""] * (len(row.fastq_urls) - len(row.fastq_md5))
    for url, md5 in zip(row.fastq_urls, md5s, strict=True):
        dest = run_scratch / Path(url).name
        source = Path(url)
        if source.exists():
            shutil.copy2(source, dest)
        else:
            http.download(url, dest, logger=logger)
        if md5:
            http.verify_md5(dest, md5)
        files.append(dest)
    return tuple(files)


# --- contract rows ----------------------------------------------------------------


def _taxa(row: ReadRow) -> tuple[str, str, str]:
    return row.family or "unknown", row.genus or "unknown", row.species or "unknown"


def named_links(contigs: dict[str, Path], link_dir: Path) -> dict[str, Path]:
    """Run accession -> a link ``<run>.fasta`` to its contigs (unique names for
    tools that key their reports by file name)."""
    if link_dir.exists():
        remove_tree(link_dir)
    link_dir.mkdir(parents=True)
    links = {}
    for run, path in contigs.items():
        link = link_dir / f"{run}.fasta"
        os.symlink(path.resolve(), link)
        links[run] = link
    return links


def classifier_for(classifier: str, gtdb_sketch: str | None) -> str | None:
    """The classifier to run: explicit, or sourmash when a sketch is configured."""
    if classifier == "none":
        return None
    sketch = gtdb_sketch or os.environ.get(GTDB_SKETCH_ENV)
    if classifier == "auto":
        return "sourmash" if sketch else None
    if not sketch:
        raise UserInputError(
            f"--classifier {classifier} needs a reference sketch (--gtdb-sketch or "
            f"{GTDB_SKETCH_ENV})."
        )
    return classifier


def assess(
    named: dict[str, Path],
    *,
    checkm2_db: str | None,
    classifier: str | None,
    gtdb_sketch: str | None,
    gtdb_lineages: str | None,
    threads: int,
    extra: dict[str, str],
    scratch: Path,
    versions: dict[str, str],
    logger: logging.Logger,
) -> tuple[dict[str, tuple[float, float]] | None, dict[str, Classification] | None]:
    """Score (CheckM2) and classify the assemblies in ``named`` (run -> FASTA).

    Either result is None when the corresponding database is not configured.
    Both are keyed by run accession.
    """
    quality: dict[str, tuple[float, float]] | None = None
    classified: dict[str, Classification] | None = None
    if checkm2_db and named:
        versions.update(preflight_checkm2())
        by_name = run_checkm2(
            list(named.values()),
            scratch / "checkm2",
            db=Path(checkm2_db),
            threads=threads,
            logger=logger,
        )
        quality = {run: by_name[link.name] for run, link in named.items() if link.name in by_name}
        for run, link in named.items():
            if link.name not in by_name:
                logger.warning(
                    "%s: CheckM2 reported no quality; the assembly is kept without quality values",
                    run,
                )
    elif not checkm2_db:
        logger.info(
            "No CheckM2 database configured (--checkm2-db or %s); assemblies are not "
            "quality-scored and --keeper quality will fall back to the tool's pick.",
            "CHECKM2DB",
        )
    if classifier and named:
        adapter = classifier_registry.create(classifier)
        versions.update(adapter.preflight())
        sketch = Path(gtdb_sketch or os.environ[GTDB_SKETCH_ENV])
        lineages = gtdb_lineages or os.environ.get(GTDB_LINEAGES_ENV)
        by_name_cls = adapter.classify(
            list(named.values()),
            scratch / "classify",
            ClassifyParams(
                db=sketch,
                lineages=None if lineages is None else Path(lineages),
                threads=threads,
                extra=dict(extra),
            ),
            logger,
        )
        classified = {
            run: by_name_cls[link.name] for run, link in named.items() if link.name in by_name_cls
        }
    return quality, classified


def apply_quality(
    outcomes: list[_Outcome],
    quality: dict[str, tuple[float, float]],
    min_completeness: float,
    max_contamination: float,
    logger: logging.Logger,
) -> None:
    """Attach CheckM2 values and excuse assemblies outside the gate."""
    for o in outcomes:
        o.quality = quality.get(o.row.run_accession)
        if o.quality is None:
            continue
        completeness, contamination = o.quality
        if completeness < min_completeness or contamination > max_contamination:
            o.excused = ExcusedRun(
                o.row.run_accession,
                "qc",
                f"qc_failed: completeness {completeness:.1f} "
                f"(min {min_completeness:g}), contamination {contamination:.1f} "
                f"(max {max_contamination:g})",
            )
            logger.warning("%s: excused, %s", o.row.run_accession, o.excused.reason)


def apply_classification(
    outcomes: list[_Outcome],
    classified: dict[str, Classification],
    versions: dict[str, str],
    logger: logging.Logger,
) -> int:
    """Name by GTDB tokens where the classifier agrees at genus; flag the rest.

    Returns the number of disagreements.
    """
    n_disagree = 0
    for o in outcomes:
        o.gtdb = classified.get(o.row.run_accession)
        if o.gtdb is None:
            continue
        versions["gtdb_sketch"] = o.gtdb.db_version or versions.get("gtdb_sketch", "")
        gtdb_tokens = _gtdb_tokens(o.gtdb.taxonomy)
        assert o.label is not None
        if gtdb_tokens[1] and gtdb_tokens[1] == o.label[1]:
            o.label = gtdb_tokens
            o.label_source = "classifier"
        else:
            o.taxonomy_flag = "classifier_disagrees"
            n_disagree += 1
            logger.warning(
                "%s: submitted as genus %s (%s) but classified as %s, %s; keeping the "
                "submitted name and flagging classifier_disagrees",
                o.row.run_accession,
                o.label[1],
                o.row.organism,
                next((c for c in o.gtdb.taxonomy.split(";") if c.startswith("g__")), "g__?"),
                o.gtdb.taxonomy.split(";")[-1],
            )
    return n_disagree


def _gtdb_tokens(lineage: str) -> tuple[str, str, str]:
    """Family, genus and species tokens from a GTDB lineage string."""
    ranks = {}
    for chunk in lineage.split(";"):
        chunk = chunk.strip()
        if len(chunk) > 3 and chunk[1:3] == "__":
            ranks[chunk[0]] = chunk[3:]
    return sanitise_taxon_tokens(ranks.get("f", ""), ranks.get("g", ""), ranks.get("s", ""))


def _name(o: _Outcome) -> str:
    assert o.label is not None
    return genome_filename(*o.label, o.row.run_accession)


def _selection_row(o: _Outcome) -> SelectionRow:
    assert o.label is not None
    family, genus, species = o.label
    completeness, contamination = o.quality if o.quality else (None, None)
    return SelectionRow(
        o.row.run_accession,
        family,
        genus,
        species,
        False,
        _name(o),
        completeness,
        contamination,
    )


def _stats_row(o: _Outcome) -> AssemblyStatsRow:
    assert o.stats is not None and o.label is not None
    coverage = o.row.bases / o.stats.total_length if o.stats.total_length else None
    completeness, contamination = o.quality if o.quality else (None, None)
    return AssemblyStatsRow(
        run_accession=o.row.run_accession,
        filename=_name(o),
        assembler=o.assembler or "",
        n_contigs=o.stats.n_contigs,
        total_length=o.stats.total_length,
        n50=o.stats.n50,
        largest_contig=o.stats.largest_contig,
        est_coverage=None if coverage is None else round(coverage, 2),
        completeness=completeness,
        contamination=contamination,
        ncbi_taxonomy=";".join(_taxa(o.row)),
        gtdb_taxonomy=o.gtdb.taxonomy if o.gtdb else "",
        label_source=o.label_source,
        taxonomy_flag=o.taxonomy_flag,
        polisher=o.polisher or "",
    )


def _append_rows(
    ctx: WorkdirContext, new_rows: list[SelectionRow], params: AssembleParams, logger
) -> list[SelectionRow]:
    """The existing selection with the new rows added (same accession replaced)."""
    from ..core.contracts import read_selection

    selection = ctx.workdir / SELECTION_TSV
    if not selection.exists():
        raise UserInputError(
            f"--append needs a working directory that already holds {SELECTION_TSV} "
            "(from metadata and genome, or ingest); run without --append to start one."
        )
    replaced = {r.accession for r in new_rows}
    kept = [r for r in read_selection(selection) if r.accession not in replaced]
    if params.outgroup is not None:
        kept = [r for r in kept if not r.is_outgroup]
        outgroup_row = _stage_outgroup(ctx, params.outgroup, logger)
        if outgroup_row is not None:
            kept.append(outgroup_row)
    logger.info("Appending %d assemblies to a selection of %d genomes", len(new_rows), len(kept))
    return [*kept, *new_rows]


def _unlink_previous(genomes_dir: Path, accession: str) -> None:
    """Remove an earlier genome file of the same run (its label may have changed)."""
    for f in genomes_dir.iterdir():
        if f.is_file() and accession_from_filename(f.name) == accession:
            f.unlink()


def _stage_outgroup(
    ctx: WorkdirContext, outgroup: str | None, logger: logging.Logger
) -> SelectionRow | None:
    acc_file = ctx.workdir / OUTGROUP_ACCESSION_TXT
    if outgroup is None:
        if ctx.outgroup_dir.exists():
            remove_tree(ctx.outgroup_dir)
        acc_file.unlink(missing_ok=True)
        return None
    source = Path(outgroup).expanduser()
    if not source.is_file():
        raise UserInputError(f"--outgroup {outgroup} is not a file.")
    family, genus, species, accession = _parse(source.name)
    if ctx.outgroup_dir.exists():
        remove_tree(ctx.outgroup_dir)
    ctx.outgroup_dir.mkdir(parents=True)
    shutil.copy2(source, ctx.outgroup_dir / source.name)
    acc_file.write_text(accession + "\n", encoding="utf-8")
    logger.info("Outgroup: %s (%s)", accession, source.name)
    return SelectionRow(accession, family, genus, species, True, source.name)


def _parse(name: str) -> tuple[str, str, str, str]:
    from ..core.contracts import parse_genome_filename

    family, genus, species, accession = parse_genome_filename(name)
    return family, genus, species, accession or accession_from_filename(name) or name
