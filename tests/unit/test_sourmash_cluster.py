"""Unit tests for the sourmash greedy clustering (numpy-vectorized)."""

from __future__ import annotations

import numpy as np

from repgenr.dereplicators.base import STATUS_CONTAINED, STATUS_REPRESENTATIVE
from repgenr.dereplicators.sourmash import _greedy_cluster


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
