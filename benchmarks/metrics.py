"""Result metrics for the scaling/bias benchmarks (stdlib only)."""

from __future__ import annotations

import json
import re
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


def newick_splits(text: str) -> tuple[frozenset[str], set[frozenset[str]]]:
    """Leaf set and non-trivial bipartitions of an unrooted newick tree.

    Each bipartition is canonicalized to the side that contains the
    alphabetically first leaf, so splits from two trees over the same leaf
    set compare directly. Support values and branch lengths are ignored.
    """
    text = text.strip().rstrip(";")
    pos = 0
    clades: list[frozenset[str]] = []

    def parse() -> frozenset[str]:
        nonlocal pos
        if text[pos] == "(":
            pos += 1
            children = [parse()]
            while text[pos] == ",":
                pos += 1
                children.append(parse())
            if text[pos] != ")":
                raise ValueError(f"malformed newick near offset {pos}")
            pos += 1
            match = re.match(r"[^,();]*", text[pos:])
            pos += match.end() if match else 0
            clade = frozenset().union(*children)
            clades.append(clade)
            return clade
        match = re.match(r"[^,();]+", text[pos:])
        if match is None:
            raise ValueError(f"malformed newick near offset {pos}")
        pos += match.end()
        return frozenset([match.group(0).split(":")[0]])

    try:
        leaves = parse()
    except IndexError:
        raise ValueError("malformed newick: unexpected end of input") from None
    first = min(leaves)
    splits = {
        clade if first in clade else leaves - clade
        for clade in clades
        if 1 < len(clade) < len(leaves) - 1
    }
    return leaves, splits


def robinson_foulds(newick_a: str, newick_b: str) -> dict[str, float]:
    """RF distance between two unrooted trees over the same leaf set."""
    leaves_a, splits_a = newick_splits(newick_a)
    leaves_b, splits_b = newick_splits(newick_b)
    if leaves_a != leaves_b:
        raise ValueError("trees have different leaf sets")
    rf = len(splits_a ^ splits_b)
    total = len(splits_a) + len(splits_b)
    return {
        "shared_splits": len(splits_a & splits_b),
        "rf_distance": rf,
        "normalized_rf": rf / total if total else 0.0,
    }


def clone_representative(
    partition: dict[str, str], truth: dict[str, str], clone_cluster: str = "clone"
) -> str | None:
    """The representative chosen for the (largest bloc of the) clone cluster."""
    clone_members = [g for g, c in truth.items() if c == clone_cluster]
    if not clone_members:
        return None
    reps = Counter(partition[g] for g in clone_members if g in partition)
    return reps.most_common(1)[0][0] if reps else None
