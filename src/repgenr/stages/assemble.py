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
from ..assemblers.base import ReadSet, registry, select_assembler
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
    genome_filename,
    read_reads,
    sanitise_taxon_tokens,
    write_assembly_stats,
    write_excused_runs,
    write_selection,
)
from ..core.errors import RepGenRError, UserInputError, WorkdirError
from ..core.executors import parallel_map
from ..core.manifest import record_from_selection
from ..core.process import check_free_disk, link_or_copy, remove_tree, staged_dir
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
    stats: ContigStats | None = None
    tool_stats: dict = field(default_factory=dict)
    excused: ExcusedRun | None = None
    quality: tuple[float, float] | None = None
    gtdb: Classification | None = None
    # Resolved filename tokens and where they came from.
    label: tuple[str, str, str] | None = None
    label_source: str = "metadata"
    taxonomy_flag: str = ""


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

    plan = _plan(rows, params, assemblies)
    versions = _preflight(plan, logger)
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
    classifier_name = classifier_for(params.classifier, params.gtdb_sketch)
    quality, classified = assess(
        named,
        checkm2_db=checkm2_db,
        classifier=classifier_name,
        gtdb_sketch=params.gtdb_sketch,
        gtdb_lineages=params.gtdb_lineages,
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

    ctx.config.record_stage(
        "assemble",
        tool=params.assembler,
        params={
            **asdict(params),
            "checkm2_db": checkm2_db,
            "classifier_effective": classifier_name,
            "assemblers_used": sorted({o.assembler for o in assembled if o.assembler}),
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


def _default_jobs(pending: list[_Outcome]) -> int:
    """Two short-read assemblies fit side by side; a long-read one wants the machine."""
    return 1 if any(o.row.platform in _LONG_READ_PLATFORMS for o in pending) else 2


# --- planning -------------------------------------------------------------------


def _plan(rows: list[ReadRow], params: AssembleParams, assemblies: Path) -> list[_Outcome]:
    """Decide, per run, whether it is done, excused up front, or to be assembled."""
    if params.assembler != "auto" and params.assembler not in registry.names():
        raise UserInputError(
            f"Unknown assembler {params.assembler!r}; available: {', '.join(registry.names())}."
        )
    plan = []
    for row in rows:
        outcome = _Outcome(row=row)
        marker = assemblies / row.run_accession / _DONE_MARKER
        if marker.exists() and (assemblies / row.run_accession / _CONTIGS_NAME).exists():
            done = json.loads(marker.read_text(encoding="utf-8"))
            outcome.assembler = done["assembler"]
            outcome.stats = ContigStats(**done["stats"])
            outcome.tool_stats = done.get("tool_stats", {})
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
        plan.append(outcome)
    return plan


def _preflight(plan: list[_Outcome], logger: logging.Logger) -> dict[str, str]:
    versions: dict[str, str] = {}
    for name in sorted({o.assembler for o in plan if o.assembler and o.excused is None}):
        versions.update(registry.create(name).preflight())
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
        stats = filter_contigs(
            result.contigs,
            out_dir / _CONTIGS_NAME,
            min_length=params.min_contig_length,
            prefix=row.run_accession,
        )
    except RepGenRError as exc:
        logger.warning("%s: assembly failed (%s)", row.run_accession, exc)
        outcome.excused = ExcusedRun(row.run_accession, "assemble", f"assembly_failed: {exc}")
        if not params.keep_files:
            remove_tree(run_scratch)
        return outcome
    if stats.n_contigs == 0:
        outcome.excused = ExcusedRun(
            row.run_accession,
            "assemble",
            f"assembly_failed: no contig of {params.min_contig_length} bp or more",
        )
        return outcome

    outcome.stats = stats
    outcome.tool_stats = dict(result.tool_stats)
    marker = {
        "assembler": outcome.assembler,
        "version": versions.get(adapter.capabilities.name, ""),
        "stats": asdict(stats),
        "tool_stats": outcome.tool_stats,
    }
    (out_dir / _DONE_MARKER).write_text(json.dumps(marker, indent=1), encoding="utf-8")
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
                logger.warning("%s: CheckM2 reported no quality", run)
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
                "%s: submitted as %s but classified as %s; keeping the submitted name",
                o.row.run_accession,
                o.row.organism,
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
