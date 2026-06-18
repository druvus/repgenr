"""Dereplication stage.

Selects a dereplicator adapter, runs it (two-stage chunking when ``--process-size``
is set and exceeded, for any tool), then normalizes the result into the canonical
contract: ``derep/representatives/`` + ``clusters.tsv`` + ``genome_status.tsv``.
Provenance (tool, params, versions) is recorded in ``repgenr.yaml`` and the
manifest derep status is updated.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass, replace
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
from ..core.executors import parallel_map
from ..core.plugins import auto_select, scale_warning
from ..core.process import link_or_copy
from ..dereplicators.base import DerepParams, DerepResult, registry

_FASTA_SUFFIXES = (".fasta", ".fasta.gz", ".fa", ".fna", ".fas")


@dataclass
class DereplicateParams:
    tool: str = "skder"
    primary_ani: float = 0.90
    secondary_ani: float = 0.99
    aligned_fraction: float = 0.50
    threads: int = 16
    # Chunk size. When set and exceeded, the two-stage chunked path runs for ANY
    # tool (native-scaling tools are single-pass by default but can be chunked
    # explicitly, e.g. to bound memory/open-file counts at 1000s-10000s genomes).
    process_size: int | None = None
    num_processes: int = 1  # parallel stage-1 chunk workers (threads split across them)
    # Stage-1 (intra-chunk) ANI thresholds. Default to the main thresholds; set
    # them looser to avoid over-collapsing within a chunk before the stage-2 pass.
    pre_primary_ani: float | None = None
    pre_secondary_ani: float | None = None
    # Taxonomy-aware reduction applied after ANI dereplication: collapse the
    # representatives to one per taxon ("species" or "genus") using the manifest
    # taxonomy. "none" keeps the ANI result as-is.
    reduce: str = "none"
    extra: dict | None = None


def run(ctx: WorkdirContext, params: DereplicateParams) -> DerepResult:
    logger = ctx.logger
    genomes = _list_genomes(ctx.genomes_dir)
    if not genomes:
        raise WorkdirError(
            f"No genome FASTAs found under {ctx.genomes_dir}. Run the genome stage first."
        )

    tool = params.tool
    if tool == "auto":
        tool = auto_select(registry, len(genomes)) or "skder"
        logger.info("Auto-selected dereplicator '%s' for %d genomes", tool, len(genomes))
    else:
        warn = scale_warning(registry, tool, len(genomes))
        if warn:
            limit, alts = warn
            logger.warning(
                "Dereplicator '%s' is tuned for <=%d genomes but you have %d; consider: %s",
                tool, limit, len(genomes), ", ".join(alts) or "none",
            )

    adapter = registry.create(tool)
    versions = adapter.preflight()
    logger.info("Dereplicating %d genomes with %s", len(genomes), tool)

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

    # Chunk whenever a process size is set and exceeded -- for ANY tool, not just
    # non-native-scaling ones. ``supports_native_scaling`` only informs the
    # default (single-pass) and the auto-select/scale warnings; ``--process-size``
    # is an explicit opt-in escape hatch for very large sets.
    needs_chunking = (
        params.process_size is not None
        and params.process_size > 0
        and len(genomes) > params.process_size
    )
    if needs_chunking:
        assert params.process_size is not None
        pre_primary = (
            params.pre_primary_ani if params.pre_primary_ani is not None else params.primary_ani
        )
        pre_secondary = (
            params.pre_secondary_ani
            if params.pre_secondary_ani is not None
            else params.secondary_ani
        )
        native = adapter.capabilities.supports_native_scaling
        logger.info(
            "Chunking %d genomes at size %d (%s)%s",
            len(genomes), params.process_size,
            "native-scaling tool, explicit chunking" if native else "non-native-scaling tool",
            ""
            if (pre_primary, pre_secondary) == (params.primary_ani, params.secondary_ani)
            else f"; stage-1 ANI pri={pre_primary} sec={pre_secondary}",
        )
        result = _dereplicate_chunked(
            adapter, genomes, scratch, derep_params, params.process_size,
            params.num_processes, pre_primary, pre_secondary, logger,
        )
    else:
        result = adapter.dereplicate(genomes, scratch, derep_params, logger)

    if params.reduce != "none":
        before = len(result.representatives)
        result = _reduce_by_taxonomy(ctx, result, params.reduce, logger)
        logger.info(
            "Taxonomy reduction (one per %s): %d -> %d representatives",
            params.reduce, before, len(result.representatives),
        )

    _write_contract(ctx, result)
    _update_manifest(ctx, result)

    ctx.config.record_stage(
        "dereplicate",
        tool=tool,
        params={
            "requested_tool": params.tool,
            "primary_ani": params.primary_ani,
            "secondary_ani": params.secondary_ani,
            "aligned_fraction": params.aligned_fraction,
            "process_size": params.process_size,
            "num_processes": params.num_processes,
            "pre_primary_ani": params.pre_primary_ani,
            "pre_secondary_ani": params.pre_secondary_ani,
            "reduce": params.reduce,
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
    num_processes: int,
    pre_primary_ani: float,
    pre_secondary_ani: float,
    logger,
) -> DerepResult:
    """Two-stage chunked dereplication, usable for any tool.

    Stage 1: dereplicate each chunk independently (up to ``num_processes`` chunks
    in parallel, the thread budget split across them), using the stage-1 ANI
    thresholds. Stage 2: dereplicate the union of stage-1 representatives with the
    final thresholds (``params``). Membership is composed so each original genome
    maps to a final representative.
    """
    chunks = [genomes[i : i + process_size] for i in range(0, len(genomes), process_size)]
    if len(chunks) > 1 and len(chunks[-1]) == 1:
        chunks[-2].extend(chunks[-1])
        chunks.pop()

    workers = max(1, min(num_processes, len(chunks)))
    threads_per_worker = max(1, params.threads // workers)
    chunk_params = replace(
        params,
        threads=threads_per_worker,
        primary_ani=pre_primary_ani,
        secondary_ani=pre_secondary_ani,
    )
    logger.info(
        "Stage 1: %d chunks, %d parallel worker(s) at %d threads each",
        len(chunks), workers, threads_per_worker,
    )

    def run_chunk(indexed: tuple[int, list[Path]]) -> DerepResult:
        idx, chunk = indexed
        chunk_dir = scratch / "intra_chunks" / f"chunk{idx}"
        logger.info("Stage 1 chunk %d/%d (%d genomes)", idx + 1, len(chunks), len(chunk))
        return adapter.dereplicate(chunk, chunk_dir, chunk_params, logger)

    stage1_results = parallel_map(run_chunk, list(enumerate(chunks)), workers, logger=logger)

    if len(stage1_results) == 1:
        return stage1_results[0]

    stage1_reps = [rep for r in stage1_results for rep in r.representatives]
    stage2 = adapter.dereplicate(stage1_reps, scratch / "inter_chunks", params, logger)

    return _compose_two_stage(stage1_results, stage2)


def _compose_two_stage(stage1: list[DerepResult], stage2: DerepResult) -> DerepResult:
    from ..dereplicators.base import STATUS_CONTAINED, STATUS_REPRESENTATIVE

    # Each stage-1 representative -> every original genome in its cluster (the rep
    # itself plus its members). Built once; lets composition run in O(N) instead
    # of scanning all members for every final rep (which was O(reps * N)).
    s1rep_to_members: dict[str, list[str]] = {}
    for r in stage1:
        for rep, members in r.clusters.items():
            s1rep_to_members[rep] = [rep, *members]

    final_clusters: dict[str, list[str]] = {}
    status: dict[str, str] = {}

    for final_rep, s1reps_contained in stage2.clusters.items():
        status[final_rep] = STATUS_REPRESENTATIVE
        contained: list[str] = []
        # final_rep's own stage-1 cluster, plus every stage-1 rep stage 2 folded
        # into it. Each original genome belongs to exactly one stage-1 cluster, so
        # these groups are disjoint -- no double counting.
        for s1rep in (final_rep, *s1reps_contained):
            for genome in s1rep_to_members.get(s1rep, [s1rep]):
                if genome == final_rep:
                    continue
                contained.append(genome)
                status[genome] = STATUS_CONTAINED
        final_clusters[final_rep] = contained

    return DerepResult(
        representatives=stage2.representatives,
        clusters=final_clusters,
        genome_status=status,
    )


def _reduce_by_taxonomy(
    ctx: WorkdirContext, result: DerepResult, level: str, logger
) -> DerepResult:
    """Collapse the ANI representatives to one per taxon (species|genus).

    Representatives sharing a manifest taxon are merged into a single keeper (the
    one with the largest existing cluster, then lexical), and the others -- plus
    their cluster members -- become contained under it. Representatives whose
    taxon is unknown/empty are kept as-is (each its own group), so reduction never
    silently drops an un-annotated genome.
    """
    from ..dereplicators.base import STATUS_CONTAINED, STATUS_REPRESENTATIVE

    taxon_of = _taxon_lookup(ctx, level)

    # group rep filename -> taxon key; unknown taxon -> unique per-rep key (kept)
    groups: dict[str, list[str]] = {}
    rep_names = [rep.name for rep in result.representatives]
    for name in rep_names:
        taxon = taxon_of.get(name, "")
        key = taxon if taxon else f"\x00{name}"  # sentinel: unknown taxa stay singletons
        groups.setdefault(key, []).append(name)

    rep_by_name = {rep.name: rep for rep in result.representatives}
    new_reps: list = []
    new_clusters: dict[str, list[str]] = {}
    status: dict[str, str] = {}

    for members in groups.values():
        # keeper: largest existing cluster, tie-break by name (deterministic)
        keeper = max(members, key=lambda n: (len(result.clusters.get(n, [])), n))
        new_reps.append(rep_by_name[keeper])
        status[keeper] = STATUS_REPRESENTATIVE
        contained: list[str] = []
        for rep_name in members:
            # the keeper's own members, plus each folded rep and its members
            own = result.clusters.get(rep_name, [])
            folded = [] if rep_name == keeper else [rep_name]
            for g in (*folded, *own):
                contained.append(g)
                status[g] = STATUS_CONTAINED
        new_clusters[keeper] = contained

    return DerepResult(
        representatives=sorted(new_reps),
        clusters=new_clusters,
        genome_status=status,
    )


def _taxon_lookup(ctx: WorkdirContext, level: str) -> dict[str, str]:
    """Map each genome filename to its taxon string at ``level`` from the manifest."""
    lookup: dict[str, str] = {}
    try:
        records = ctx.manifest.all_genomes(include_outgroup=True)
    except Exception:  # no manifest (e.g. tests) -> empty taxonomy, everything stays a singleton
        return lookup
    for rec in records:
        if not rec.filename:
            continue
        taxon = rec.species if level == "species" else rec.genus
        lookup[rec.filename] = taxon or ""
    return lookup


def _write_contract(ctx: WorkdirContext, result: DerepResult) -> None:
    rep_dir = ctx.representatives_dir
    if rep_dir.exists():
        shutil.rmtree(rep_dir)
    rep_dir.mkdir(parents=True, exist_ok=True)

    for rep in result.representatives:
        source = rep if rep.exists() else ctx.genomes_dir / rep.name
        if not source.exists():
            raise WorkdirError(f"Representative genome file missing: {rep.name}")
        link_or_copy(source, rep_dir / rep.name)

    write_clusters(ctx.derep_dir / CLUSTERS_TSV, result.clusters)
    write_genome_status(ctx.derep_dir / GENOME_STATUS_TSV, result.genome_status)


def _update_manifest(ctx: WorkdirContext, result: DerepResult) -> None:
    manifest = ctx.manifest
    rep_by_member: dict[str, str] = {}
    for rep, members in result.clusters.items():
        for m in members:
            rep_by_member[m] = rep
    updates: list[tuple[str, str, str | None]] = []
    for genome, status in result.genome_status.items():
        accession = _accession_from_filename(genome)
        if accession is None:
            continue
        representative = None
        if status == "contained":
            rep_file = rep_by_member.get(genome)
            representative = _accession_from_filename(rep_file) if rep_file else None
        updates.append((accession, status, representative))
    # One batched transaction instead of a commit per genome. Accessions absent
    # from the manifest are no-op UPDATEs (not errors), so we no longer swallow
    # exceptions -- a real database error should surface, not corrupt the manifest
    # silently.
    manifest.set_derep_status_many(updates)


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
