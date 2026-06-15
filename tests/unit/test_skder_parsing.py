"""Offline tests for skDER output parsing + edge-based membership."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.dereplicators.skder import _parse_skder_output

_LOG = logging.getLogger("test")


def _make_skder_out(tmp_path: Path, edges: str) -> Path:
    out = tmp_path / "skder_out"
    reps = out / "Dereplicated_Representative_Genomes"
    reps.mkdir(parents=True)
    for name in ("repA.fasta", "repB.fasta"):
        (reps / name).write_text(">x\nACGT\n")
    (out / "Skani_Triangle_Edge_Output.txt").write_text(edges)
    return out


def test_membership_from_edges(tmp_path: Path) -> None:
    # memX is closest (highest ANI) to repA; memY only similar to repB
    edges = (
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "/g/repA.fasta\t/g/memX.fasta\t99.5\t90\t88\trA\tmX\n"
        "/g/repB.fasta\t/g/memX.fasta\t97.0\t80\t80\trB\tmX\n"
        "/g/memY.fasta\t/g/repB.fasta\t99.2\t85\t85\tmY\trB\n"
    )
    out = _make_skder_out(tmp_path, edges)
    genomes = [Path(f"/g/{n}.fasta") for n in ("repA", "repB", "memX", "memY")]

    result = _parse_skder_output(out, genomes, ani_cutoff=99.0, af_cutoff=50.0, logger=_LOG)

    assert {p.name for p in result.representatives} == {"repA.fasta", "repB.fasta"}
    # memX assigned to repA (99.5 > 97.0); memY assigned to repB
    assert "memX.fasta" in result.clusters["repA.fasta"]
    assert "memY.fasta" in result.clusters["repB.fasta"]
    assert result.genome_status["memX.fasta"] == "contained"
    assert result.genome_status["repA.fasta"] == "representative"


def test_edges_below_cutoff_not_clustered(tmp_path: Path) -> None:
    edges = (
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "/g/repA.fasta\t/g/memX.fasta\t90.0\t90\t90\trA\tmX\n"  # ANI below 99 cutoff
    )
    out = _make_skder_out(tmp_path, edges)
    genomes = [Path(f"/g/{n}.fasta") for n in ("repA", "repB", "memX")]
    result = _parse_skder_output(out, genomes, ani_cutoff=99.0, af_cutoff=50.0, logger=_LOG)
    # memX not assigned to any cluster, but still marked contained
    assert result.clusters["repA.fasta"] == []
    assert result.genome_status["memX.fasta"] == "contained"
