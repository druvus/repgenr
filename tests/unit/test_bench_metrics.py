"""Tests for benchmarks.metrics (scoring benchmark cells against truth.json)."""

from __future__ import annotations

import json

import pytest

from benchmarks.metrics import (
    adjusted_rand_index,
    clone_representative,
    load_truth,
    newick_splits,
    partition_from_clusters_tsv,
    robinson_foulds,
)


def test_load_truth_round_trip(tmp_path):
    clusters = {"a.fasta": "c1", "b.fasta": "c1", "c.fasta": "clone"}
    (tmp_path / "truth.json").write_text(
        json.dumps({"clusters": clusters, "seed": 1}), encoding="utf-8"
    )
    assert load_truth(tmp_path) == clusters


def test_partition_from_clusters_tsv(tmp_path):
    from repgenr.core.contracts import write_clusters

    tsv = tmp_path / "clusters.tsv"
    write_clusters(tsv, {"rep1.fasta": ["m1.fasta", "m2.fasta"], "rep2.fasta": []})
    partition = partition_from_clusters_tsv(tsv)
    assert partition["rep1.fasta"] == "rep1.fasta"
    assert partition["m1.fasta"] == "rep1.fasta"
    assert partition["m2.fasta"] == "rep1.fasta"
    assert partition["rep2.fasta"] == "rep2.fasta"


def test_ari_identical_labelings_is_one():
    a = {"g1": "x", "g2": "x", "g3": "y", "g4": "y"}
    b = {"g1": "p", "g2": "p", "g3": "q", "g4": "q"}  # same partition, renamed
    assert adjusted_rand_index(a, b) == 1.0


def test_ari_disagreement_is_below_one():
    a = {"g1": "x", "g2": "x", "g3": "y", "g4": "y"}
    b = {"g1": "x", "g2": "y", "g3": "y", "g4": "x"}
    assert adjusted_rand_index(a, b) < 0.5


def test_ari_uses_only_common_keys():
    a = {"g1": "x", "g2": "x", "g3": "y", "g4": "y", "extra": "z"}
    b = {"g1": "p", "g2": "p", "g3": "q", "g4": "q"}
    assert adjusted_rand_index(a, b) == 1.0


def test_ari_known_small_value():
    # Two clusterings of 6 items with one item moved between clusters.
    a = {f"g{i}": ("A" if i < 3 else "B") for i in range(6)}
    b = {f"g{i}": ("A" if i < 2 else "B") for i in range(6)}
    value = adjusted_rand_index(a, b)
    assert 0.0 < value < 1.0
    assert abs(value - adjusted_rand_index(b, a)) < 1e-12  # symmetric


def test_clone_representative_majority_bloc():
    truth = {"a": "clone", "b": "clone", "c": "clone", "d": "bg"}
    partition = {"a": "rep1", "b": "rep1", "c": "rep2", "d": "rep3"}
    assert clone_representative(partition, truth) == "rep1"


def test_clone_representative_none_without_clone_cluster():
    truth = {"a": "c1", "b": "c2"}
    partition = {"a": "a", "b": "b"}
    assert clone_representative(partition, truth) is None


def test_newick_splits_leaves_and_bipartitions():
    # unrooted 5-leaf tree with two non-trivial splits: {a,b} and {d,e}
    tree = "((a:0.1,b:0.2)0.9:0.05,c:0.3,(d:0.1,e:0.1)0.8:0.02);"
    leaves, splits = newick_splits(tree)
    assert leaves == frozenset("abcde")
    # splits are canonicalized to the side containing the first leaf "a"
    assert splits == {frozenset("ab"), frozenset("abc")}


def test_newick_splits_ignores_support_and_lengths():
    with_annotations = "((a:0.1,b:0.2)0.95:0.01,(c:0.3,d:0.4)0.5:0.02,e:0.5);"
    bare = "((a,b),(c,d),e);"
    assert newick_splits(with_annotations) == newick_splits(bare)


def test_robinson_foulds_identical_topologies():
    a = "((a:0.1,b:0.2):0.1,c:0.3,(d:0.1,e:0.2):0.1);"
    b = "((b,a),c,(e,d));"  # same topology, rotated
    result = robinson_foulds(a, b)
    assert result["rf_distance"] == 0
    assert result["normalized_rf"] == 0.0
    assert result["shared_splits"] == 2


def test_robinson_foulds_disjoint_topologies():
    a = "((a,b),(c,d),(e,f));"
    b = "((a,c),(b,e),(d,f));"
    result = robinson_foulds(a, b)
    assert result["shared_splits"] == 0
    assert result["normalized_rf"] == 1.0


def test_robinson_foulds_rejects_different_leaf_sets():
    with pytest.raises(ValueError, match="leaf sets"):
        robinson_foulds("((a,b),c,d);", "((a,b),c,e);")


def test_newick_splits_rejects_malformed_input():
    with pytest.raises(ValueError, match="malformed"):
        newick_splits("((a,b,c);")
