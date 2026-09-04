"""Quality-aware representative selection applied after any dereplicator.

Dereplicators pick representatives by connectivity, first-seen order or tool
defaults, which favours the most-sequenced genotype in a cluster. This step
re-picks each cluster's representative by assembly quality when the manifest
carries CheckM-style completeness and contamination, and leaves the adapter's
choice in place otherwise.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path

from ..dereplicators.base import STATUS_CONTAINED, STATUS_REPRESENTATIVE, DerepResult

CONTAMINATION_WEIGHT = 5.0


def quality_score(completeness: float, contamination: float) -> float:
    """dRep-style score: completeness minus five times contamination."""
    return completeness - CONTAMINATION_WEIGHT * contamination


def rescore_representatives(
    result: DerepResult,
    quality: Mapping[str, tuple[float, float]],
    logger: logging.Logger,
) -> tuple[DerepResult, int]:
    """Return a copy of ``result`` whose representatives are the best-scoring
    cluster members, and the number of clusters whose representative changed.

    A member replaces the representative only when its score is strictly higher;
    genomes without quality never replace a scored representative, and a cluster
    with no scored genome keeps the adapter's choice.
    """
    if not quality:
        return result, 0

    rep_paths = {p.name: p for p in result.representatives}
    new_reps: list[Path] = []
    new_clusters: dict[str, list[str]] = {}
    status = dict(result.genome_status)
    swaps = 0

    def score(name: str) -> float | None:
        q = quality.get(name)
        return None if q is None else quality_score(*q)

    for rep_name, members in result.clusters.items():
        candidates = [rep_name, *members]
        current = score(rep_name)
        best_name, best = rep_name, current
        for m in members:
            s = score(m)
            if s is None:
                continue
            if best is None or s > best:
                best_name, best = m, s
        if best_name == rep_name:
            new_reps.append(rep_paths[rep_name])
            new_clusters[rep_name] = list(members)
            continue
        swaps += 1
        new_reps.append(rep_paths[rep_name].with_name(best_name))
        new_clusters[best_name] = sorted(n for n in candidates if n != best_name)
        status[best_name] = STATUS_REPRESENTATIVE
        status[rep_name] = STATUS_CONTAINED
        logger.info(
            "Keeper: %s replaces %s (score %.2f vs %s)",
            best_name,
            rep_name,
            best,
            "n/a" if current is None else f"{current:.2f}",
        )

    if swaps:
        logger.info(
            "Quality-aware keeper changed %d of %d representatives", swaps, len(new_clusters)
        )
    return DerepResult(
        representatives=sorted(new_reps),
        clusters=new_clusters,
        genome_status=status,
        genome_information=result.genome_information,
    ), swaps
