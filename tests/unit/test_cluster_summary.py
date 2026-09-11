"""Per-cluster summary: one row per representative over clusters.tsv + quality."""

from __future__ import annotations

from pathlib import Path

from repgenr.core.contracts import read_cluster_summary, write_cluster_summary
from repgenr.stages.cluster_summary import summarise_clusters

REP = "Fam_Gen_tularensis_GCF_1.1.fasta"
M1 = "Fam_Gen_tularensis_GCF_2.1.fasta"
M2 = "Fam_Gen_holarctica_GCF_3.1.fasta"
SOLO = "Fam_Gen_novicida_GCF_4.1.fasta"


def test_summary_counts_members_and_species() -> None:
    rows = summarise_clusters({REP: [M1, M2], SOLO: []}, {})
    by_rep = {r.representative: r for r in rows}
    assert by_rep[REP].n_members == 2
    assert by_rep[REP].n_species == 2
    assert by_rep[REP].species == "tularensis,holarctica"
    assert by_rep[SOLO].n_members == 0
    assert by_rep[SOLO].n_species == 1


def test_summary_sorted_by_size_then_name() -> None:
    rows = summarise_clusters({"Zed_Gen_x_GCF_9.1.fasta": [], SOLO: [], REP: [M1, M2]}, {})
    assert [r.representative for r in rows] == [REP, SOLO, "Zed_Gen_x_GCF_9.1.fasta"]


def test_summary_quality_columns_blank_without_quality() -> None:
    (row,) = summarise_clusters({REP: [M1]}, {})
    assert row.rep_completeness is None
    assert row.rep_contamination is None
    assert row.member_max_completeness is None
    assert row.member_min_contamination is None
    assert row.best_member == ""


def test_summary_best_member_is_representative_when_keeper_is_best() -> None:
    quality = {REP: (99.0, 0.1), M1: (95.0, 1.0), M2: (98.0, 0.5)}
    (row,) = summarise_clusters({REP: [M1, M2]}, quality)
    assert row.rep_completeness == 99.0
    assert row.rep_contamination == 0.1
    assert row.member_max_completeness == 98.0
    assert row.member_min_contamination == 0.5
    assert row.best_member == REP


def test_summary_flags_better_member_than_keeper() -> None:
    quality = {REP: (90.0, 3.0), M1: (99.0, 0.2)}
    (row,) = summarise_clusters({REP: [M1]}, quality)
    assert row.best_member == M1


def test_summary_unscored_members_are_ignored_in_extremes() -> None:
    quality = {REP: (99.0, 0.1), M1: (95.0, 1.0)}
    (row,) = summarise_clusters({REP: [M1, M2]}, quality)
    assert row.member_max_completeness == 95.0
    assert row.member_min_contamination == 1.0
    assert row.n_members == 2


def test_summary_round_trip(tmp_path: Path) -> None:
    quality = {REP: (99.0, 0.1), M1: (95.0, 1.0)}
    rows = summarise_clusters({REP: [M1, M2], SOLO: []}, quality)
    path = tmp_path / "cluster_summary.tsv"
    write_cluster_summary(path, rows)
    header = path.read_text().splitlines()[0].split("\t")
    assert header == [
        "representative",
        "n_members",
        "n_species",
        "species",
        "rep_completeness",
        "rep_contamination",
        "member_max_completeness",
        "member_min_contamination",
        "best_member",
    ]
    assert read_cluster_summary(path) == rows
