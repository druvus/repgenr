"""Result metrics for the scaling/bias benchmarks (stdlib only)."""

from __future__ import annotations

import json
from collections import Counter
from math import comb
from pathlib import Path


def load_truth(set_dir: Path) -> dict[str, str]:
    """Genome filename -> true cluster id, from the generator's truth.json."""
    truth = json.loads((Path(set_dir) / "truth.json").read_text(encoding="utf-8"))
    return dict(truth["clusters"])


def partition_from_clusters_tsv(clusters_tsv: Path) -> dict[str, str]:
    """Genome filename -> representative filename, from a derep clusters.tsv."""
    from repgenr.core.contracts import read_clusters

    partition: dict[str, str] = {}
    for rep, members in read_clusters(Path(clusters_tsv)).items():
        partition[rep] = rep
        for member in members:
            partition[member] = rep
    return partition


def adjusted_rand_index(a: dict[str, str], b: dict[str, str]) -> float:
    """ARI between two labelings over their common keys (1.0 = identical)."""
    keys = sorted(set(a) & set(b))
    if len(keys) < 2:
        return 1.0
    pair_counts: Counter[tuple[str, str]] = Counter((a[k], b[k]) for k in keys)
    a_counts: Counter[str] = Counter(a[k] for k in keys)
    b_counts: Counter[str] = Counter(b[k] for k in keys)

    index = sum(comb(c, 2) for c in pair_counts.values())
    sum_a = sum(comb(c, 2) for c in a_counts.values())
    sum_b = sum(comb(c, 2) for c in b_counts.values())
    total = comb(len(keys), 2)
    expected = sum_a * sum_b / total
    maximum = (sum_a + sum_b) / 2
    if maximum == expected:
        return 1.0
    return (index - expected) / (maximum - expected)


def clone_representative(
    partition: dict[str, str], truth: dict[str, str], clone_cluster: str = "clone"
) -> str | None:
    """The representative chosen for the (largest bloc of the) clone cluster."""
    clone_members = [g for g, c in truth.items() if c == clone_cluster]
    if not clone_members:
        return None
    reps = Counter(partition[g] for g in clone_members if g in partition)
    return reps.most_common(1)[0][0] if reps else None
