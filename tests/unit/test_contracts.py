"""Contract round-trip tests."""

from __future__ import annotations

from pathlib import Path

from repgenr.core.contracts import (
    read_clusters,
    write_clusters,
    write_genome_status,
    write_genomes_map,
    write_tree2tax,
)


def test_clusters_round_trip(tmp_path: Path) -> None:
    clusters = {
        "rep_a.fasta": ["m1.fasta", "m2.fasta"],
        "rep_b.fasta": [],
    }
    path = tmp_path / "clusters.tsv"
    write_clusters(path, clusters)

    # representative always lists itself in the file
    lines = path.read_text().splitlines()
    assert lines[0] == "representative\tmember"
    assert "rep_a.fasta\trep_a.fasta" in lines
    assert "rep_b.fasta\trep_b.fasta" in lines

    # reading back excludes the self-edge
    assert read_clusters(path) == clusters


def test_genome_status_sorted(tmp_path: Path) -> None:
    path = tmp_path / "genome_status.tsv"
    write_genome_status(path, {"z.fasta": "contained", "a.fasta": "representative"})
    rows = [line.split("\t") for line in path.read_text().splitlines()[1:]]
    assert rows == [["a.fasta", "representative"], ["z.fasta", "contained"]]


def test_tree2tax_dedupes(tmp_path: Path) -> None:
    path = tmp_path / "tree2tax.tsv"
    write_tree2tax(path, [("leaf", "n1"), ("leaf", "n1"), ("n1", "root")])
    body = path.read_text().splitlines()[1:]
    assert body == ["leaf\tn1", "n1\troot"]


def test_genomes_map(tmp_path: Path) -> None:
    path = tmp_path / "genomes_map.tsv"
    write_genomes_map(path, [("GCA_000001", "Fam_gen_sp_GCA_000001")])
    assert path.read_text().strip() == "GCA_000001\tFam_gen_sp_GCA_000001"
