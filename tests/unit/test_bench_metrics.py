"""Tests for benchmarks.metrics (scoring benchmark cells against truth.json)."""

from __future__ import annotations

import json

from benchmarks.metrics import (
    adjusted_rand_index,
    clone_representative,
    load_truth,
    partition_from_clusters_tsv,
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
