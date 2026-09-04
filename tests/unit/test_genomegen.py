"""Synthetic genome-set generator for the scaling/bias benchmarks."""

from __future__ import annotations

import json
from pathlib import Path

from benchmarks.genomegen import generate_set, pairwise_identity


def _read_seq(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    )


def test_deterministic_output(tmp_path: Path) -> None:
    a = tmp_path / "a"
    b = tmp_path / "b"
    generate_set(a, scenario="balanced", n=6, seed=42, genome_length=5000)
    generate_set(b, scenario="balanced", n=6, seed=42, genome_length=5000)
    files_a = sorted(p.name for p in a.glob("*.fasta"))
    files_b = sorted(p.name for p in b.glob("*.fasta"))
    assert files_a == files_b and len(files_a) == 6
    for name in files_a:
        assert _read_seq(a / name) == _read_seq(b / name)


def test_truth_json_written(tmp_path: Path) -> None:
    truth = generate_set(
        tmp_path / "s", scenario="clonal", n=10, seed=1, genome_length=4000, clone_fraction=0.5
    )
    on_disk = json.loads((tmp_path / "s" / "truth.json").read_text(encoding="utf-8"))
    assert on_disk == truth
    assert on_disk["scenario"] == "clonal"
    assert len(on_disk["clusters"]) == 10
    clone = on_disk["clone_cluster"]
    clone_members = [g for g, c in on_disk["clusters"].items() if c == clone]
    assert len(clone_members) == 5  # clone_fraction=0.5 of n=10


def test_intra_cluster_identity_higher_than_inter(tmp_path: Path) -> None:
    out = tmp_path / "s"
    truth = generate_set(out, scenario="balanced", n=8, seed=3, genome_length=20000)
    by_cluster: dict[str, list[str]] = {}
    for name, cluster in truth["clusters"].items():
        by_cluster.setdefault(cluster, []).append(name)
    clusters = [members for members in by_cluster.values() if len(members) >= 2]
    assert len(clusters) >= 2, "balanced scenario must produce multiple clusters"

    a1, a2 = clusters[0][0], clusters[0][1]
    b1 = clusters[1][0]
    seq = {n: _read_seq(out / n) for n in (a1, a2, b1)}
    intra = pairwise_identity(seq[a1], seq[a2])
    inter = pairwise_identity(seq[a1], seq[b1])
    assert intra > 0.99, f"intra-cluster identity too low: {intra}"
    assert inter < intra - 0.01, f"clusters not separated: intra={intra} inter={inter}"


def test_clone_block_is_near_identical(tmp_path: Path) -> None:
    out = tmp_path / "s"
    truth = generate_set(
        out, scenario="clonal", n=10, seed=7, genome_length=50000, clone_fraction=0.4
    )
    clone = truth["clone_cluster"]
    members = [g for g, c in truth["clusters"].items() if c == clone]
    s1, s2 = _read_seq(out / members[0]), _read_seq(out / members[1])
    assert pairwise_identity(s1, s2) >= 0.9995


def test_accession_order_controls_sort_position(tmp_path: Path) -> None:
    """clustered: the clone block occupies the alphabetically-first filenames;
    random: it does not (for a seeded set where that is checkable)."""
    clustered = generate_set(
        tmp_path / "c",
        scenario="clonal",
        n=10,
        seed=5,
        genome_length=2000,
        clone_fraction=0.4,
        order="clustered",
    )
    names = sorted(clustered["clusters"])
    clone = clustered["clone_cluster"]
    first_four = names[:4]
    assert all(clustered["clusters"][n] == clone for n in first_four)

    rand = generate_set(
        tmp_path / "r",
        scenario="clonal",
        n=10,
        seed=5,
        genome_length=2000,
        clone_fraction=0.4,
        order="random",
    )
    names_r = sorted(rand["clusters"])
    assert not all(rand["clusters"][n] == rand["clone_cluster"] for n in names_r[:4])


def test_canonical_filenames(tmp_path: Path) -> None:
    from repgenr.core.contracts import parse_genome_filename

    generate_set(tmp_path / "s", scenario="mixed", n=5, seed=2, genome_length=2000)
    for p in (tmp_path / "s").glob("*.fasta"):
        family, genus, species, accession = parse_genome_filename(p.name)
        assert family and genus and species and accession


def test_pairwise_identity_helper() -> None:
    assert pairwise_identity("ACGT", "ACGT") == 1.0
    assert pairwise_identity("ACGT", "ACGA") == 0.75
