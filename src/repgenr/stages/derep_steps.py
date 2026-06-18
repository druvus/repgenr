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
import shutil
from dataclasses import dataclass
from pathlib import Path

from ..core.contracts import (
    CLUSTERS_TSV,
    GENOME_STATUS_TSV,
    read_clusters,
    write_clusters,
    write_genome_status,
)
from ..core.errors import WorkdirError
from ..core.process import link_or_copy
from ..dereplicators.base import DerepParams, DerepResult, registry
from .dereplicate import _compose_two_stage

_REPRESENTATIVES_DIR = "representatives"


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
    adapter.preflight()
    derep_params = DerepParams(
        primary_ani=params.primary_ani,
        secondary_ani=params.secondary_ani,
        aligned_fraction=params.aligned_fraction,
        threads=params.threads,
        extra=dict(params.extra or {}),
    )
    scratch = _fresh(params.out_dir / "scratch")
    result = adapter.dereplicate(params.genomes, scratch, derep_params, logger)

    fallbacks = sorted({g.parent for g in params.genomes})
    _write_step_contract(params.out_dir, result, fallbacks)
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
    adapter.preflight()
    derep_params = DerepParams(
        primary_ani=params.primary_ani,
        secondary_ani=params.secondary_ani,
        aligned_fraction=params.aligned_fraction,
        threads=params.threads,
        extra=dict(params.extra or {}),
    )
    scratch = _fresh(params.out_dir / "scratch")
    stage2 = adapter.dereplicate(union, scratch, derep_params, logger)
    final = _compose_two_stage(stage1, stage2)

    # The final representatives are stage-2 representative paths, which live in the
    # chunk representatives/ directories; fall back to those when resolving files.
    fallbacks = [d / _REPRESENTATIVES_DIR for d in params.chunk_dirs]
    _write_step_contract(params.out_dir, final, fallbacks)
    logger.info(
        "dereplicate-merge: %d chunks, union of %d reps -> %d representatives (%s)",
        len(params.chunk_dirs), len(union), len(final.representatives), params.tool,
    )
    return final


def _load_chunk(chunk_dir: Path) -> DerepResult:
    """Read a chunk result directory back into a DerepResult (clusters + rep paths)."""
    clusters_path = chunk_dir / CLUSTERS_TSV
    if not clusters_path.exists():
        raise WorkdirError(f"dereplicate-merge: {clusters_path} not found (not a chunk result).")
    clusters = read_clusters(clusters_path)
    rep_dir = chunk_dir / _REPRESENTATIVES_DIR
    reps: list[Path] = []
    for rep_name in clusters:
        rep_path = rep_dir / rep_name
        if not rep_path.exists():
            raise WorkdirError(
                f"dereplicate-merge: representative file missing in chunk: {rep_path}"
            )
        reps.append(rep_path)
    return DerepResult(representatives=reps, clusters=clusters, genome_status={})


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
