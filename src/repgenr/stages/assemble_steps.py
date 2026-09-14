"""Stateless reads-chain steps for data-channel orchestration.

The shared-workdir ``assemble`` stage fetches, assembles, scores, classifies
and publishes every run in one process. The Nextflow layer runs the same work
as discrete tasks, one per run for the assembly and one batch for the quality
and classification pass, so each is exposed here as a step that reads explicit
inputs and writes a result directory: no manifest, no shared workdir.

``assemble_run``   one run from ``reads.tsv`` -> ``<out>/contigs.fasta`` and
                   ``<out>/assembly.ok``, or ``<out>/excused_runs.tsv``.
``genome_qc``      a directory of such run directories -> ``quality.tsv``
                   (CheckM2) and ``classification.tsv`` (classifier).
``reads_gather``   ``reads.tsv`` + the run directories + the optional QC
                   directory -> ``genomes/``, ``selection.tsv``,
                   ``assembly_stats.tsv``, ``excused_runs.tsv`` and an empty
                   ``outgroup_accession.txt`` (the genome contract).
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path

from ..classifiers.base import Classification
from ..core.contracts import (
    ASSEMBLY_STATS_TSV,
    EXCUSED_RUNS_TSV,
    SELECTION_TSV,
    ExcusedRun,
    read_excused_runs,
    read_reads,
    write_assembly_stats,
    write_excused_runs,
    write_selection,
)
from ..core.errors import UserInputError, WorkdirError
from ..core.process import link_or_copy, remove_tree
from ..core.versions import write_versions_fragment
from . import assemble as stage
from .assemble import AssembleParams
from .assemble_qc import checkm2_db_from_env
from .ingest import OUTGROUP_ACCESSION_TXT

QUALITY_TSV = "quality.tsv"
CLASSIFICATION_TSV = "classification.tsv"
_CONTIGS = stage._CONTIGS_NAME
_MARKER = stage._DONE_MARKER


# --- assemble-run -------------------------------------------------------------------


@dataclass
class AssembleRunParams:
    reads_tsv: Path
    run: str
    out_dir: Path
    assembler: str = "auto"
    threads: int = 16
    memory_gb: int = 16
    min_contig_length: int = 500
    keep_reads: bool = False
    keep_files: bool = False
    extra: dict[str, str] = field(default_factory=dict)
    versions_out: Path | None = None


def assemble_run(params: AssembleRunParams, logger: logging.Logger) -> bool:
    """Fetch and assemble one run of ``reads.tsv`` into ``out_dir``.

    Returns True when contigs were written, False when the run was excused
    (recorded in ``<out_dir>/excused_runs.tsv``).
    """
    if not params.reads_tsv.exists():
        raise WorkdirError(f"assemble-run: reads file not found: {params.reads_tsv}")
    rows = [r for r in read_reads(params.reads_tsv) if r.run_accession == params.run]
    if not rows:
        raise UserInputError(f"assemble-run: run {params.run} is not listed in {params.reads_tsv}")
    stage_params = AssembleParams(
        assembler=params.assembler,
        threads=params.threads,
        memory_gb=params.memory_gb,
        min_contig_length=params.min_contig_length,
        keep_reads=params.keep_reads,
        keep_files=params.keep_files,
        extra=dict(params.extra),
    )
    out_dir = params.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    [outcome] = stage._plan(rows, stage_params, out_dir.parent)
    versions = stage._preflight([outcome], logger)
    if outcome.excused is None and outcome.stats is None:
        scratch = out_dir.parent / f"{out_dir.name}.scratch"
        outcome = stage._fetch_and_assemble(
            outcome, stage_params, params.threads, out_dir, scratch, versions, logger
        )
        if scratch.exists() and not (params.keep_files or params.keep_reads):
            remove_tree(scratch)
    if params.versions_out is not None:
        write_versions_fragment(params.versions_out, versions)
    if outcome.excused is not None:
        write_excused_runs(out_dir / EXCUSED_RUNS_TSV, [outcome.excused])
        logger.info("assemble-run: %s excused (%s)", params.run, outcome.excused.reason)
        return False
    return True


# --- genome-qc --------------------------------------------------------------------------


@dataclass
class GenomeQcParams:
    assemblies_dir: Path
    out_dir: Path
    threads: int = 16
    checkm2_db: str | None = None
    classifier: str = "auto"
    gtdb_sketch: str | None = None
    gtdb_lineages: str | None = None
    extra: dict[str, str] = field(default_factory=dict)
    versions_out: Path | None = None


def _assembled_contigs(assemblies_dir: Path) -> dict[str, Path]:
    """Run accession -> contigs file, for every run directory holding a done marker."""
    if not assemblies_dir.is_dir():
        raise WorkdirError(f"assemblies directory not found: {assemblies_dir}")
    found = {}
    for run_dir in sorted(p for p in assemblies_dir.iterdir() if p.is_dir()):
        if (run_dir / _MARKER).exists() and (run_dir / _CONTIGS).exists():
            found[run_dir.name] = run_dir / _CONTIGS
    return found


def genome_qc(params: GenomeQcParams, logger: logging.Logger) -> int:
    """Score and classify every finished assembly under ``assemblies_dir``.

    Writes ``quality.tsv`` when a CheckM2 database is configured and
    ``classification.tsv`` when a classifier runs. Returns the number of
    assemblies assessed.
    """
    checkm2_db = params.checkm2_db or checkm2_db_from_env()
    classifier = stage.classifier_for(params.classifier, params.gtdb_sketch)
    if not checkm2_db and not classifier:
        raise UserInputError(
            "genome-qc needs a CheckM2 database (--checkm2-db or CHECKM2DB) and/or a "
            f"reference sketch (--gtdb-sketch or {stage.GTDB_SKETCH_ENV})."
        )
    contigs = _assembled_contigs(params.assemblies_dir)
    out = params.out_dir
    out.mkdir(parents=True, exist_ok=True)
    scratch = out / "scratch"
    named = stage.named_links(contigs, scratch / "named")
    versions: dict[str, str] = {}
    quality, classified = stage.assess(
        named,
        checkm2_db=checkm2_db,
        classifier=classifier,
        gtdb_sketch=params.gtdb_sketch,
        gtdb_lineages=params.gtdb_lineages,
        threads=params.threads,
        extra=params.extra,
        scratch=scratch,
        versions=versions,
        logger=logger,
    )
    if quality is not None:
        write_quality(out / QUALITY_TSV, quality)
    if classified is not None:
        write_classification(out / CLASSIFICATION_TSV, classified)
    remove_tree(scratch)
    if params.versions_out is not None:
        write_versions_fragment(params.versions_out, versions)
    logger.info(
        "genome-qc: %d assemblies, CheckM2 %s, classifier %s",
        len(contigs),
        "on" if quality is not None else "off",
        classifier or "off",
    )
    return len(contigs)


def write_quality(path: Path, quality: dict[str, tuple[float, float]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["run_accession", "completeness", "contamination"])
        for run, (completeness, contamination) in sorted(quality.items()):
            w.writerow([run, f"{completeness:.2f}", f"{contamination:.2f}"])


def read_quality(path: Path) -> dict[str, tuple[float, float]]:
    with path.open(encoding="utf-8", newline="") as fi:
        return {
            row["run_accession"]: (float(row["completeness"]), float(row["contamination"]))
            for row in csv.DictReader(fi, delimiter="\t")
        }


def write_classification(path: Path, classified: dict[str, Classification]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["run_accession", "taxonomy", "rank", "score", "db_version"])
        for run, c in sorted(classified.items()):
            w.writerow([run, c.taxonomy, c.rank, "" if c.score is None else c.score, c.db_version])


def read_classification(path: Path) -> dict[str, Classification]:
    with path.open(encoding="utf-8", newline="") as fi:
        return {
            row["run_accession"]: Classification(
                taxonomy=row["taxonomy"],
                rank=row["rank"],
                score=float(row["score"]) if row["score"] else None,
                db_version=row["db_version"],
            )
            for row in csv.DictReader(fi, delimiter="\t")
        }


# --- reads-gather -----------------------------------------------------------------------


@dataclass
class ReadsGatherParams:
    reads_tsv: Path
    assemblies_dir: Path
    out_dir: Path
    qc_dir: Path | None = None
    min_completeness: float = 50.0
    max_contamination: float = 10.0


def reads_gather(params: ReadsGatherParams, logger: logging.Logger) -> int:
    """Assemble the genome contract from the per-run results and the QC tables.

    Returns the number of accepted genomes.
    """
    if not params.reads_tsv.exists():
        raise WorkdirError(f"reads-gather: reads file not found: {params.reads_tsv}")
    rows = read_reads(params.reads_tsv)
    if not rows:
        raise WorkdirError(f"reads-gather: {params.reads_tsv} lists no runs.")
    if not params.assemblies_dir.is_dir():
        raise WorkdirError(f"reads-gather: assemblies directory not found: {params.assemblies_dir}")

    # Rebuild the per-run outcomes from the markers and the excuse files.
    outcomes = stage._plan(rows, AssembleParams(), params.assemblies_dir)
    for o in outcomes:
        excuse_file = params.assemblies_dir / o.row.run_accession / EXCUSED_RUNS_TSV
        if o.stats is None and excuse_file.exists():
            o.excused = read_excused_runs(excuse_file)[0]
        elif o.stats is None and o.excused is None:
            o.excused = ExcusedRun(
                o.row.run_accession, "assemble", "assembly_failed: no result directory"
            )
    assembled = [o for o in outcomes if o.stats is not None and o.stats.n_contigs > 0]
    for o in assembled:
        o.label = stage._taxa(o.row)

    versions: dict[str, str] = {}
    n_disagree = 0
    if params.qc_dir is not None:
        quality_path = params.qc_dir / QUALITY_TSV
        if quality_path.exists():
            stage.apply_quality(
                assembled,
                read_quality(quality_path),
                params.min_completeness,
                params.max_contamination,
                logger,
            )
            assembled = [o for o in assembled if o.excused is None]
        classification_path = params.qc_dir / CLASSIFICATION_TSV
        if classification_path.exists():
            n_disagree = stage.apply_classification(
                assembled, read_classification(classification_path), versions, logger
            )

    out = params.out_dir
    out.mkdir(parents=True, exist_ok=True)
    excused = [o.excused for o in outcomes if o.excused is not None]
    if excused:
        write_excused_runs(out / EXCUSED_RUNS_TSV, excused)
    if not assembled:
        raise WorkdirError(
            f"None of the {len(rows)} runs produced an accepted assembly; see {EXCUSED_RUNS_TSV}."
        )
    genomes = out / "genomes"
    if genomes.exists():
        remove_tree(genomes)
    genomes.mkdir(parents=True)
    for o in assembled:
        link_or_copy(
            params.assemblies_dir / o.row.run_accession / _CONTIGS, genomes / stage._name(o)
        )
    write_selection(out / SELECTION_TSV, [stage._selection_row(o) for o in assembled])
    write_assembly_stats(out / ASSEMBLY_STATS_TSV, [stage._stats_row(o) for o in assembled])
    # The standalone reads chain has no outgroup here; downstream steps read an
    # empty accession file as "none".
    (out / OUTGROUP_ACCESSION_TXT).write_text("", encoding="utf-8")
    logger.info(
        "reads-gather: %d of %d runs accepted (%d excused, %d classifier disagreements)",
        len(assembled),
        len(rows),
        len(excused),
        n_disagree,
    )
    return len(assembled)
