"""Unit tests for the sourmash greedy clustering (numpy-vectorized)."""

from __future__ import annotations

import numpy as np

from repgenr.dereplicators.base import STATUS_CONTAINED, STATUS_REPRESENTATIVE
from repgenr.dereplicators.sourmash import (
    _greedy_cluster,
    _parse_pairwise_csv,
    _sparse_greedy_cluster,
)


def _name_map(labels):
    return {label: label for label in labels}


def test_two_disjoint_pairs() -> None:
    labels = ["a", "b", "c", "d"]
    # a~b and c~d are similar; cross-pairs are not.
    sim = [
        [1.0, 0.999, 0.1, 0.1],
        [0.999, 1.0, 0.1, 0.1],
        [0.1, 0.1, 1.0, 0.999],
        [0.1, 0.1, 0.999, 1.0],
    ]
    clusters, status = _greedy_cluster(labels, sim, _name_map(labels), 0.99)
    # one representative per pair; members are the other element
    assert len(clusters) == 2
    assert all(len(m) == 1 for m in clusters.values())
    assert set(clusters) == {"a", "c"}  # first index of each pair wins on ties
    assert clusters["a"] == ["b"] and clusters["c"] == ["d"]
    assert status["a"] == STATUS_REPRESENTATIVE and status["b"] == STATUS_CONTAINED


def test_all_connected_collapses_to_one() -> None:
    labels = ["a", "b", "c"]
    sim = np.array([[1.0, 0.999, 0.999], [0.999, 1.0, 0.999], [0.999, 0.999, 1.0]])
    clusters, status = _greedy_cluster(labels, sim, _name_map(labels), 0.99)
    assert len(clusters) == 1
    rep = next(iter(clusters))
    assert sorted(clusters[rep]) == ["b", "c"]


def test_all_distinct_keeps_all() -> None:
    labels = ["a", "b", "c"]
    sim = np.eye(3)  # nothing above threshold off-diagonal
    clusters, status = _greedy_cluster(labels, sim, _name_map(labels), 0.99)
    assert len(clusters) == 3
    assert all(m == [] for m in clusters.values())
    assert all(s == STATUS_REPRESENTATIVE for s in status.values())


# --- sparse (branchwater pairwise) path ---------------------------------------

_PAIRWISE_HEADER = (
    "query_name,query_md5,match_name,match_md5,containment,max_containment,"
    "jaccard,average_containment,intersect_hashes,ksize,scaled,moltype,cosine\n"
)


def _pairwise_csv(tmp_path, rows):
    # rows: (query, match, jaccard) -- containment is set to jaccard so the file
    # looks like a real (already containment-thresholded) branchwater edge list.
    path = tmp_path / "pairwise.csv"
    lines = [_PAIRWISE_HEADER]
    for q, m, j in rows:
        lines.append(f"{q},md5q,{m},md5m,{j},{j},{j},{j},10,31,1000,DNA,\n")
    path.write_text("".join(lines))
    return path


def test_parse_pairwise_csv_thresholds_and_symmetrizes(tmp_path) -> None:
    known = {"a", "b", "c", "d"}
    csv_path = _pairwise_csv(
        tmp_path,
        [
            ("a", "b", 0.999),   # kept
            ("a", "a", 1.0),     # self-edge dropped
            ("c", "d", 0.5),     # below threshold dropped
            ("a", "z", 0.999),   # unknown node dropped
        ],
    )
    neighbors = _parse_pairwise_csv(csv_path, 0.99, known)
    assert neighbors == {"a": {"b"}, "b": {"a"}}  # symmetric, single kept edge


def test_sparse_two_disjoint_pairs() -> None:
    labels = ["a", "b", "c", "d"]
    neighbors = {"a": {"b"}, "b": {"a"}, "c": {"d"}, "d": {"c"}}
    clusters, status = _sparse_greedy_cluster(labels, neighbors, _name_map(labels))
    assert set(clusters) == {"a", "c"}  # first of each pair wins on ties
    assert clusters["a"] == ["b"] and clusters["c"] == ["d"]
    assert status["a"] == STATUS_REPRESENTATIVE and status["b"] == STATUS_CONTAINED


def test_sparse_loner_stays_singleton() -> None:
    labels = ["a", "b", "c"]
    neighbors = {"a": {"b"}, "b": {"a"}}  # c has no edges
    clusters, status = _sparse_greedy_cluster(labels, neighbors, _name_map(labels))
    assert set(clusters) == {"a", "c"}
    assert clusters["c"] == []
    assert status["c"] == STATUS_REPRESENTATIVE


def test_sparse_matches_dense_on_same_graph() -> None:
    # Build one similarity matrix; derive both a dense matrix and the equivalent
    # sparse adjacency from it, and confirm both clusterers agree.
    labels = ["a", "b", "c", "d", "e"]
    edges = {("a", "b"), ("b", "c"), ("d", "e")}  # a-b-c chain + d-e pair
    thr = 0.99
    n = len(labels)
    idx = {lab: i for i, lab in enumerate(labels)}
    sim = np.eye(n)
    neighbors: dict[str, set[str]] = {}
    for u, v in edges:
        sim[idx[u], idx[v]] = sim[idx[v], idx[u]] = 0.999
        neighbors.setdefault(u, set()).add(v)
        neighbors.setdefault(v, set()).add(u)

    dense, _ = _greedy_cluster(labels, sim, _name_map(labels), thr)
    sparse, _ = _sparse_greedy_cluster(labels, neighbors, _name_map(labels))
    assert {k: sorted(v) for k, v in dense.items()} == {
        k: sorted(v) for k, v in sparse.items()
    }
