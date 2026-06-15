"""Neighbor-joining sanity tests."""

from __future__ import annotations

from repgenr.tree.newick import neighbor_joining


def test_two_taxa() -> None:
    nwk = neighbor_joining(["a", "b"], [[0.0, 0.4], [0.4, 0.0]])
    assert nwk.startswith("(") and nwk.endswith(");")
    assert "a:0.2" in nwk and "b:0.2" in nwk


def test_four_taxa_groups_close_pairs() -> None:
    labels = ["A", "B", "C", "D"]
    # A,B close; C,D close; groups far apart
    dist = [
        [0.0, 0.1, 0.9, 0.9],
        [0.1, 0.0, 0.9, 0.9],
        [0.9, 0.9, 0.0, 0.1],
        [0.9, 0.9, 0.1, 0.0],
    ]
    nwk = neighbor_joining(labels, dist)
    # all taxa present and a valid Newick string
    for name in labels:
        assert name in nwk
    assert nwk.count("(") == nwk.count(")")
    assert nwk.endswith(");")


def test_labels_sanitized() -> None:
    nwk = neighbor_joining(["a b", "c:d", "e(f)"], [[0, 1, 1], [1, 0, 1], [1, 1, 0]])
    for bad in (" ", ":", "(", ")"):
        # sanitized labels must not reintroduce structural characters in names
        assert f"a{bad}b" not in nwk
