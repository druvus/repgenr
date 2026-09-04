"""Discrete dereplication steps for scatter-gather orchestration.

The shared-workdir ``dereplicate`` stage runs the whole two-stage pipeline in one
process. For horizontal scaling (Nextflow scatter-gather across nodes) the same
work is exposed as two stateless, file-in / file-out steps:

* :func:`dereplicate_chunk` dereplicates one chunk of genomes and writes a chunk
  result directory in the canonical contract (``representatives/`` +
  ``clusters.tsv`` + ``genome_status.tsv``).
* :func:`dereplicate_merge` takes several chunk result directories, dereplicates
  the union of their representatives with the final thresholds, and composes the
  two-stage membership into a final contract directory.

Neither step touches the SQLite manifest or the run config: their entire state is
the directories they read and write, so Nextflow can stage them between tasks and
cache them independently. A chunk result directory is itself a valid contract, so
``dereplicate_merge`` consumes exactly what ``dereplicate_chunk`` produces. This
mirrors the in-process ``_dereplicate_chunked`` reduce-tree, split into steps the
orchestrator scatters; multi-level reduction is expressed by feeding a merge
output back through another chunk/merge round.
"""

from __future__ import annotations

import logging
import os
import shutil
from dataclasses import dataclass
from pathlib import Path

from ..core.contracts import (
    CLUSTERS_TSV,
    GENOME_STATUS_TSV,
    read_clusters,
    read_genome_status,
    read_selection,
    write_clusters,
    write_genome_status,
)
from ..core.errors import WorkdirError
from ..core.plugins import warn_unconsumed_extras
from ..core.process import link_or_copy
from ..core.versions import write_versions_fragment
from ..dereplicators.base import DerepParams, DerepResult, check_result_complete, registry
from .derep_keeper import rescore_representatives
from .dereplicate import _compose_two_stage

_REPRESENTATIVES_DIR = "representatives"


def _maybe_write_versions(path: Path | None, versions: dict[str, str]) -> None:
    if path is not None:
        write_versions_fragment(path, versions)


def _quality_from_selection(selection_tsv: Path) -> dict[str, tuple[float, float]]:
    """filename -> (completeness, contamination) for rows carrying both values."""
    return {
        r.filename: (r.completeness, r.contamination)
        for r in read_selection(selection_tsv)
        if r.completeness is not None and r.contamination is not None
    }


@dataclass
class ChunkParams:
    tool: str
    genomes: list[Path]
    out_dir: Path
    primary_ani: float = 0.90
    secondary_ani: float = 0.99
    aligned_fraction: float = 0.50
    threads: int = 16
    extra: dict | None = None
    versions_out: Path | None = None
    # A selection.tsv (with quality columns) that enables the quality-aware
    # keeper below. Unlike at the merge step, every genome in the chunk has a
    # real file here (params.genomes), so any promotion is always resolvable.
    selection_tsv: Path | None = None
    keeper: str = "quality"  # quality | tool


@dataclass
class MergeParams:
    tool: str
    chunk_dirs: list[Path]
    out_dir: Path
    primary_ani: float = 0.90
    secondary_ani: float = 0.99
    aligned_fraction: float = 0.50
    threads: int = 16
    extra: dict | None = None
    versions_out: Path | None = None
    # A selection.tsv (with quality columns) that enables the quality-aware
    # keeper below; without it (or with keeper="tool") the adapter's own
    # merge-level pick stands, as before.
    selection_tsv: Path | None = None
    keeper: str = "quality"  # quality | tool


def dereplicate_chunk(params: ChunkParams, logger: logging.Logger) -> DerepResult:
    """Dereplicate one chunk of genomes; write its contract to ``out_dir``."""
    if not params.genomes:
        raise WorkdirError("dereplicate-chunk: no genome paths provided.")
    missing = [g for g in params.genomes if not g.exists()]
    if missing:
        raise WorkdirError(
            f"dereplicate-chunk: {len(missing)} genome file(s) not found, e.g. {missing[0]}"
        )

    adapter = registry.create(params.tool)
    caps = adapter.capabilities
    warn_unconsumed_extras(caps, params.extra or {}, logger, family="Dereplicator")
    versions = adapter.preflight()
    _maybe_write_versions(params.versions_out, versions)
    derep_params = DerepParams(
        primary_ani=params.primary_ani,
        secondary_ani=params.secondary_ani,
        aligned_fraction=params.aligned_fraction,
        threads=params.threads,
        extra={**caps.default_params, **(params.extra or {})},
    )
    scratch = _fresh(params.out_dir / "scratch")
    result = adapter.dereplicate(params.genomes, scratch, derep_params, logger)

    if params.keeper == "quality" and params.selection_tsv is not None:
        quality = _quality_from_selection(params.selection_tsv)
        result, keeper_swaps = rescore_representatives(result, quality, logger)
        if keeper_swaps:
            logger.info(
                "dereplicate-chunk: quality-aware keeper changed %d representative(s)",
                keeper_swaps,
            )
    elif params.keeper == "tool" and params.selection_tsv is not None:
        logger.info(
            "dereplicate-chunk: --selection-tsv given but --keeper is 'tool'; "
            "keeping the adapter's own representatives."
        )

    check_result_complete(result, [g.name for g in params.genomes])
    fallbacks = sorted({g.parent for g in params.genomes})
    _write_step_contract(params.out_dir, result, fallbacks)
    shutil.rmtree(scratch, ignore_errors=True)  # drop tool intermediates from the output
    logger.info(
        "dereplicate-chunk: %d genomes -> %d representatives (%s)",
        len(params.genomes), len(result.representatives), params.tool,
    )
    return result


def dereplicate_merge(params: MergeParams, logger: logging.Logger) -> DerepResult:
    """Dereplicate the union of chunk representatives; compose the final contract."""
    if not params.chunk_dirs:
        raise WorkdirError("dereplicate-merge: no chunk directories provided.")
    stage1 = [_load_chunk(d) for d in params.chunk_dirs]
    union = [rep for r in stage1 for rep in r.representatives]
    if not union:
        raise WorkdirError("dereplicate-merge: the chunk directories hold no representatives.")

    adapter = registry.create(params.tool)
    caps = adapter.capabilities
    warn_unconsumed_extras(caps, params.extra or {}, logger, family="Dereplicator")
    versions = adapter.preflight()
    _maybe_write_versions(params.versions_out, versions)
    derep_params = DerepParams(
        primary_ani=params.primary_ani,
        secondary_ani=params.secondary_ani,
        aligned_fraction=params.aligned_fraction,
        threads=params.threads,
        extra={**caps.default_params, **(params.extra or {})},
    )
    scratch = _fresh(params.out_dir / "scratch")
    stage2 = adapter.dereplicate(union, scratch, derep_params, logger)
    final = _compose_two_stage(stage1, stage2)

    if params.keeper == "quality" and params.selection_tsv is not None:
        # Only stage-1 representatives are staged at the merge step (each
        # chunk's representatives/ directory); a chunk-level CONTAINED member's
        # file never reaches this step, even though _compose_two_stage's
        # expanded membership lists it in final.clusters. Restrict the quality
        # map to resolvable names so such a member can never be promoted to a
        # representative _write_step_contract cannot then find a file for.
        resolvable = {rep.name for r in stage1 for rep in r.representatives}
        quality = {
            name: qual
            for name, qual in _quality_from_selection(params.selection_tsv).items()
            if name in resolvable
        }
        final, keeper_swaps = rescore_representatives(final, quality, logger)
        if keeper_swaps:
            logger.info(
                "dereplicate-merge: quality-aware keeper changed %d representative(s)",
                keeper_swaps,
            )
    elif params.keeper == "tool" and params.selection_tsv is not None:
        logger.info(
            "dereplicate-merge: --selection-tsv given but --keeper is 'tool'; "
            "keeping the adapter's own representatives."
        )

    # Every genome the chunks saw: cluster members plus the genomes that carry a
    # status without a cluster (QC rejects), so the completeness check covers the
    # whole input set rather than only the clustered part of it.
    all_names = {
        name
        for r in stage1
        for rep, members in r.clusters.items()
        for name in (rep, *members)
    }
    all_names |= {genome for r in stage1 for genome in r.genome_status}
    check_result_complete(final, all_names)

    # The final representatives are stage-2 representative paths, which live in the
    # chunk representatives/ directories; fall back to those when resolving files.
    fallbacks = [d / _REPRESENTATIVES_DIR for d in params.chunk_dirs]
    _write_step_contract(params.out_dir, final, fallbacks)
    shutil.rmtree(scratch, ignore_errors=True)  # drop tool intermediates from the output
    logger.info(
        "dereplicate-merge: %d chunks, union of %d reps -> %d representatives (%s)",
        len(params.chunk_dirs), len(union), len(final.representatives), params.tool,
    )
    return final


def _load_chunk(chunk_dir: Path) -> DerepResult:
    """Read a chunk result directory back into a DerepResult.

    Reads the per-genome statuses as well as the clusters: genomes the chunk
    rejected on QC appear only in ``genome_status.tsv``, and dropping them here
    would lose them from the merged contract.
    """
    clusters_path = chunk_dir / CLUSTERS_TSV
    if not clusters_path.exists():
        raise WorkdirError(f"dereplicate-merge: {clusters_path} not found (not a chunk result).")
    clusters = read_clusters(clusters_path)
    status_path = chunk_dir / GENOME_STATUS_TSV
    if not status_path.exists():
        raise WorkdirError(f"dereplicate-merge: {status_path} not found (not a chunk result).")
    status = read_genome_status(status_path)
    rep_dir = chunk_dir / _REPRESENTATIVES_DIR
    reps: list[Path] = []
    for rep_name in clusters:
        rep_path = rep_dir / rep_name
        if not rep_path.exists():
            raise WorkdirError(
                f"dereplicate-merge: representative file missing in chunk: {rep_path}"
            )
        reps.append(rep_path)
    return DerepResult(representatives=reps, clusters=clusters, genome_status=status)


def _write_step_contract(
    out_dir: Path, result: DerepResult, fallback_dirs: list[Path]
) -> None:
    """Write representatives/ + clusters.tsv + genome_status.tsv under ``out_dir``."""
    rep_dir = out_dir / _REPRESENTATIVES_DIR
    if rep_dir.exists():
        shutil.rmtree(rep_dir)
    rep_dir.mkdir(parents=True, exist_ok=True)
    for rep in result.representatives:
        source = rep if rep.exists() else _find(fallback_dirs, rep.name)
        if source is None:
            raise WorkdirError(f"Representative genome file missing: {rep.name}")
        if os.path.getsize(os.path.realpath(source)) == 0:
            raise WorkdirError(
                f"Representative genome is empty: {rep.name}. The source file "
                f"({source}) has zero length -- an upstream download or staging "
                "step produced an empty genome."
            )
        link_or_copy(source, rep_dir / rep.name)
    write_clusters(out_dir / CLUSTERS_TSV, result.clusters)
    write_genome_status(out_dir / GENOME_STATUS_TSV, result.genome_status)


def _find(dirs: list[Path], name: str) -> Path | None:
    for d in dirs:
        candidate = d / name
        if candidate.exists():
            return candidate
    return None


def _fresh(path: Path) -> Path:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)
    return path
