"""Small shared helpers for the live suite (imported by test modules)."""

from __future__ import annotations

import json
from pathlib import Path

from repgenr.core.contracts import CLUSTERS_TSV, read_clusters

Partition = set[frozenset[str]]


def truth_of(genomes_dir: Path) -> dict:
    return json.loads((genomes_dir / "truth.json").read_text(encoding="utf-8"))


def truth_partition(truth: dict) -> Partition:
    groups: dict[str, set[str]] = {}
    for name, cluster in truth["clusters"].items():
        groups.setdefault(cluster, set()).add(name)
    return {frozenset(g) for g in groups.values()}


def derep_partition(workdir: Path) -> Partition:
    clusters = read_clusters(workdir / "derep" / CLUSTERS_TSV)
    return {frozenset({rep, *members}) for rep, members in clusters.items()}


def representatives(workdir: Path) -> set[str]:
    return set(read_clusters(workdir / "derep" / CLUSTERS_TSV))


def newick_leaves(newick: str) -> set[str]:
    """Leaf labels of a Newick string (internal labels and lengths ignored)."""
    out: set[str] = set()
    token = ""
    leaf = False  # a token that starts right after "(" or "," is a leaf
    for ch in newick:
        if ch in "(),;":
            if leaf and token:
                out.add(token.split(":")[0].strip())
            token = ""
            leaf = ch in "(,"
        else:
            token += ch
    return {x for x in out if x}


def log_text(workdir: Path) -> str:
    return (workdir / "repgenr.log").read_text(encoding="utf-8")
