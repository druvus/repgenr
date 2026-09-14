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
import shutil
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..assemblers.base import AssembleParams as AdapterParams
from ..assemblers.base import ReadSet, registry, select_assembler
from ..assemblers.contigs import ContigStats, filter_contigs
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
    write_assembly_stats,
    write_excused_runs,
    write_selection,
)
from ..core.errors import RepGenRError, UserInputError, WorkdirError
from ..core.executors import parallel_map
from ..core.manifest import record_from_selection
from ..core.process import check_free_disk, link_or_copy, remove_tree, staged_dir
from .ingest import OUTGROUP_ACCESSION_TXT

_DONE_MARKER = "assembly.ok"
_CONTIGS_NAME = "contigs.fasta"
# Rough peak disk per run: the FASTQ files, their uncompressed form and the
# assembler's scratch.
_DISK_FACTOR = 4


@dataclass
class AssembleParams:
    assembler: str = "auto"
    threads: int = 16
    # Concurrent runs; memory, not CPU, bounds this, so it stays small.
    jobs: int = 2
    memory_gb: int = 16
    min_contig_length: int = 500
    # A FASTA file to set aside as the outgroup (ingest semantics).
    outgroup: str | None = None
    keep_reads: bool = False
    keep_files: bool = False
    # Tool tuning from ``--tool-arg``; each adapter declares the keys it reads.
    extra: dict = field(default_factory=dict)


@dataclass
class _Outcome:
    row: ReadRow
    assembler: str | None = None
    stats: ContigStats | None = None
    tool_stats: dict = field(default_factory=dict)
    excused: ExcusedRun | None = None


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
    jobs = max(1, min(params.jobs, len(pending) or 1))
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
            outcome, params, threads_each, assemblies, scratch, versions, logger
        )

    done = parallel_map(work, pending, jobs, logger=logger)
    by_run = {o.row.run_accession: o for o in [*plan, *done]}
    outcomes = [by_run[r.run_accession] for r in rows]

    assembled = [o for o in outcomes if o.stats is not None and o.stats.n_contigs > 0]
    excused = [o.excused for o in outcomes if o.excused is not None]
    excused_path = ctx.workdir / EXCUSED_RUNS_TSV
    if excused:
        write_excused_runs(excused_path, excused)
    else:
        excused_path.unlink(missing_ok=True)
    if not assembled:
        raise WorkdirError(
            f"None of the {len(rows)} runs produced an assembly; see {EXCUSED_RUNS_TSV}."
        )

    selection_rows = [_selection_row(o.row) for o in assembled]
    outgroup_row = _stage_outgroup(ctx, params.outgroup, logger)
    if outgroup_row is not None:
        selection_rows.append(outgroup_row)

    with staged_dir(ctx.genomes_dir) as genomes_dir:
        for o in assembled:
            link_or_copy(
                assemblies / o.row.run_accession / _CONTIGS_NAME,
                genomes_dir / genome_filename(*_taxa(o.row), o.row.run_accession),
            )
    write_selection(ctx.workdir / SELECTION_TSV, selection_rows)
    write_assembly_stats(ctx.workdir / ASSEMBLY_STATS_TSV, [_stats_row(o) for o in assembled])
    ctx.manifest.replace_genomes([record_from_selection(r, "sra") for r in selection_rows])

    ctx.config.record_stage(
        "assemble",
        tool=params.assembler,
        params={
            **asdict(params),
            "assemblers_used": sorted({o.assembler for o in assembled if o.assembler}),
            "n_assembled": len(assembled),
            "n_excused": len(excused),
            "outgroup_accession": outgroup_row.accession if outgroup_row else None,
        },
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info(
        "Assembled %d of %d runs (%d excused); wrote %s",
        len(assembled),
        len(rows),
        len(excused),
        SELECTION_TSV,
    )
    return len(assembled)


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
    assemblies: Path,
    scratch: Path,
    versions: dict[str, str],
    logger: logging.Logger,
) -> _Outcome:
    row = outcome.row
    run_scratch = scratch / row.run_accession
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
    out_dir = assemblies / row.run_accession
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


def _selection_row(row: ReadRow) -> SelectionRow:
    family, genus, species = _taxa(row)
    return SelectionRow(
        row.run_accession,
        family,
        genus,
        species,
        False,
        genome_filename(family, genus, species, row.run_accession),
    )


def _stats_row(o: _Outcome) -> AssemblyStatsRow:
    assert o.stats is not None
    family, genus, species = _taxa(o.row)
    coverage = o.row.bases / o.stats.total_length if o.stats.total_length else None
    return AssemblyStatsRow(
        run_accession=o.row.run_accession,
        filename=genome_filename(family, genus, species, o.row.run_accession),
        assembler=o.assembler or "",
        n_contigs=o.stats.n_contigs,
        total_length=o.stats.total_length,
        n50=o.stats.n50,
        largest_contig=o.stats.largest_contig,
        est_coverage=None if coverage is None else round(coverage, 2),
        ncbi_taxonomy=";".join((family, genus, species)),
        label_source="metadata",
    )


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
