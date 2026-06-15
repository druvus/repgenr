"""Dereplication stage.

Selects a dereplicator adapter, runs it (chunking large sets when the adapter
does not scale natively), then normalizes the result into the canonical
contract: ``derep/representatives/`` + ``clusters.tsv`` + ``genome_status.tsv``.
Provenance (tool, params, versions) is recorded in ``repgenr.yaml`` and the
manifest derep status is updated.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import (
    CLUSTERS_TSV,
    GENOME_STATUS_TSV,
    write_clusters,
    write_genome_status,
)
from ..core.errors import WorkdirError
from ..dereplicators.base import DerepParams, DerepResult, registry

_FASTA_SUFFIXES = (".fasta", ".fasta.gz", ".fa", ".fna", ".fas")


@dataclass
class DereplicateParams:
    tool: str = "skder"
    primary_ani: float = 0.90
    secondary_ani: float = 0.99
    aligned_fraction: float = 0.50
    threads: int = 16
    process_size: int | None = None  # chunk size for non-native-scaling tools
    extra: dict | None = None


def run(ctx: WorkdirContext, params: DereplicateParams) -> DerepResult:
    logger = ctx.logger
    genomes = _list_genomes(ctx.genomes_dir)
    if not genomes:
        raise WorkdirError(
            f"No genome FASTAs found under {ctx.genomes_dir}. Run the genome stage first."
        )

    adapter = registry.create(params.tool)
    versions = adapter.preflight()
    logger.info("Dereplicating %d genomes with %s", len(genomes), params.tool)

    derep_params = DerepParams(
        primary_ani=params.primary_ani,
        secondary_ani=params.secondary_ani,
        aligned_fraction=params.aligned_fraction,
        threads=params.threads,
        extra=dict(params.extra or {}),
    )

    scratch = ctx.scratch_dir / "dereplicate"
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True, exist_ok=True)

    needs_chunking = (
        params.process_size
        and not adapter.capabilities.supports_native_scaling
        and len(genomes) > params.process_size
    )
    if needs_chunking:
        assert params.process_size is not None
        logger.info(
            "Tool %s does not scale natively; chunking %d genomes at size %d",
            params.tool, len(genomes), params.process_size,
        )
        result = _dereplicate_chunked(
            adapter, genomes, scratch, derep_params, params.process_size, logger
        )
    else:
        result = adapter.dereplicate(genomes, scratch, derep_params, logger)

    _write_contract(ctx, result)
    _update_manifest(ctx, result)

    ctx.config.record_stage(
        "dereplicate",
        tool=params.tool,
        params={
            "primary_ani": params.primary_ani,
            "secondary_ani": params.secondary_ani,
            "aligned_fraction": params.aligned_fraction,
            "process_size": params.process_size,
            **(params.extra or {}),
        },
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info(
        "Dereplication complete: %d representatives of %d genomes",
        len(result.representatives), len(genomes),
    )
    return result


def _list_genomes(genomes_dir: Path) -> list[Path]:
    if not genomes_dir.exists():
        return []
    return sorted(
        p for p in genomes_dir.iterdir()
        if not p.name.startswith(".") and any(p.name.endswith(s) for s in _FASTA_SUFFIXES)
    )


def _dereplicate_chunked(
    adapter,
    genomes: list[Path],
    scratch: Path,
    params: DerepParams,
    process_size: int,
    logger,
) -> DerepResult:
    """Two-stage chunked dereplication for tools that don't scale natively.

    Stage 1: dereplicate each chunk independently. Stage 2: dereplicate the
    union of stage-1 representatives. Membership is composed so each original
    genome maps to a final representative.
    """
    chunks = [genomes[i : i + process_size] for i in range(0, len(genomes), process_size)]
    if len(chunks) > 1 and len(chunks[-1]) == 1:
        chunks[-2].extend(chunks[-1])
        chunks.pop()

    stage1_results: list[DerepResult] = []
    for idx, chunk in enumerate(chunks):
        chunk_dir = scratch / "intra_chunks" / f"chunk{idx}"
        logger.info("Stage 1 chunk %d/%d (%d genomes)", idx + 1, len(chunks), len(chunk))
        stage1_results.append(adapter.dereplicate(chunk, chunk_dir, params, logger))

    if len(stage1_results) == 1:
        return stage1_results[0]

    stage1_reps = [rep for r in stage1_results for rep in r.representatives]
    stage2 = adapter.dereplicate(stage1_reps, scratch / "inter_chunks", params, logger)

    return _compose_two_stage(stage1_results, stage2)


def _compose_two_stage(stage1: list[DerepResult], stage2: DerepResult) -> DerepResult:
    # member -> stage1 representative
    member_to_s1rep: dict[str, str] = {}
    for r in stage1:
        for rep, members in r.clusters.items():
            member_to_s1rep[rep] = rep
            for m in members:
                member_to_s1rep[m] = rep

    final_clusters: dict[str, list[str]] = {}
    status: dict[str, str] = {}
    from ..dereplicators.base import STATUS_CONTAINED, STATUS_REPRESENTATIVE

    for final_rep, s1reps_contained in stage2.clusters.items():
        final_clusters[final_rep] = []
        status[final_rep] = STATUS_REPRESENTATIVE
        # absorb everything that pointed at final_rep in stage 1
        for member, s1rep in member_to_s1rep.items():
            if s1rep == final_rep and member != final_rep:
                final_clusters[final_rep].append(member)
                status[member] = STATUS_CONTAINED
        # absorb members of the other stage-1 reps that stage 2 folded in
        for s1rep in s1reps_contained:
            for member, mapped in member_to_s1rep.items():
                if mapped == s1rep:
                    final_clusters[final_rep].append(member)
                    status[member] = STATUS_CONTAINED

    return DerepResult(
        representatives=stage2.representatives,
        clusters=final_clusters,
        genome_status=status,
    )


def _write_contract(ctx: WorkdirContext, result: DerepResult) -> None:
    rep_dir = ctx.representatives_dir
    if rep_dir.exists():
        shutil.rmtree(rep_dir)
    rep_dir.mkdir(parents=True, exist_ok=True)

    for rep in result.representatives:
        source = rep if rep.exists() else ctx.genomes_dir / rep.name
        if not source.exists():
            raise WorkdirError(f"Representative genome file missing: {rep.name}")
        shutil.copy2(source, rep_dir / rep.name)

    write_clusters(ctx.derep_dir / CLUSTERS_TSV, result.clusters)
    write_genome_status(ctx.derep_dir / GENOME_STATUS_TSV, result.genome_status)


def _update_manifest(ctx: WorkdirContext, result: DerepResult) -> None:
    manifest = ctx.manifest
    rep_by_member: dict[str, str] = {}
    for rep, members in result.clusters.items():
        for m in members:
            rep_by_member[m] = rep
    for genome, status in result.genome_status.items():
        accession = _accession_from_filename(genome)
        if accession is None:
            continue
        representative = None
        if status == "contained":
            rep_file = rep_by_member.get(genome)
            representative = _accession_from_filename(rep_file) if rep_file else None
        try:
            manifest.set_derep_status(accession, status, representative)
        except Exception:  # genome may not be in manifest (e.g. tests) -- non-fatal
            pass


def _accession_from_filename(filename: str | None) -> str | None:
    """Genome files are named ``..._<GCx>_<digits>.fasta``; recover the accession."""
    if not filename:
        return None
    stem = filename
    for suffix in _FASTA_SUFFIXES:
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    parts = stem.split("_")
    if len(parts) >= 2:
        return f"{parts[-2]}_{parts[-1]}"
    return None
