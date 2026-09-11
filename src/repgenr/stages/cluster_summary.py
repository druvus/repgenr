"""cluster_summary stage: one row per dereplication representative.

Summarises ``derep/clusters.tsv`` into ``derep/cluster_summary.tsv`` with the
cluster size, the species it spans (parsed from the canonical filenames) and
the manifest CheckM quality of the keeper against its members. The dereplicate
stage writes the file itself; this stage regenerates it for an existing
working directory without rerunning the dereplicator.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import (
    CLUSTER_SUMMARY_TSV,
    CLUSTERS_TSV,
    ClusterSummaryRow,
    parse_genome_filename,
    read_clusters,
    write_cluster_summary,
)
from ..core.errors import WorkdirError
from .derep_keeper import quality_score

Quality = Mapping[str, tuple[float, float]]


@dataclass
class ClusterSummaryParams:
    pass


def summarise_clusters(
    clusters: Mapping[str, list[str]], quality: Quality
) -> list[ClusterSummaryRow]:
    """Build the summary rows, largest cluster first, then by representative name."""
    rows = [_summarise(rep, members, quality) for rep, members in clusters.items()]
    rows.sort(key=lambda r: (-r.n_members, r.representative))
    return rows


def _summarise(rep: str, members: list[str], quality: Quality) -> ClusterSummaryRow:
    others = [m for m in members if m != rep]
    species: list[str] = []
    for name in (rep, *others):
        sp = parse_genome_filename(name)[2]
        if sp not in species:
            species.append(sp)

    rep_q = quality.get(rep)
    scored = [(m, quality[m]) for m in others if m in quality]
    best_member = ""
    candidates = [(rep, rep_q)] if rep_q is not None else []
    candidates.extend(scored)
    if candidates:
        # ``max`` keeps the first of equal scores, so a tie leaves the keeper.
        best_member = max(candidates, key=lambda c: quality_score(*c[1]))[0]

    return ClusterSummaryRow(
        representative=rep,
        n_members=len(others),
        n_species=len(species),
        species=",".join(species),
        rep_completeness=None if rep_q is None else rep_q[0],
        rep_contamination=None if rep_q is None else rep_q[1],
        member_max_completeness=max((q[0] for _, q in scored), default=None),
        member_min_contamination=min((q[1] for _, q in scored), default=None),
        best_member=best_member,
    )


def run(ctx: WorkdirContext, params: ClusterSummaryParams) -> Path:
    logger = ctx.logger
    clusters_file = ctx.derep_dir / CLUSTERS_TSV
    if not clusters_file.exists():
        raise WorkdirError(f"Missing {clusters_file}. Run the dereplicate stage first.")
    clusters = read_clusters(clusters_file)
    from .dereplicate import _quality_lookup

    quality = _quality_lookup(ctx)
    if not quality:
        logger.info("No assembly quality in the manifest; quality columns are left blank")
    rows = summarise_clusters(clusters, quality)
    out = ctx.derep_dir / CLUSTER_SUMMARY_TSV
    write_cluster_summary(out, rows)
    ctx.config.record_stage(
        "cluster_summary",
        params={**asdict(params), "clusters": len(rows)},
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Summarised %d clusters into %s", len(rows), out)
    return out
