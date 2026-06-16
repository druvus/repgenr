"""Tests for MAF -> reference-anchored MSA-FASTA."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.converters.maf_to_fasta import maf_to_fasta

_MAF = """##maf version=1
a score=0
s ref.chr1 0 8 + 8 ACGTACGT
s qry.chr1 0 8 + 8 ACGAACGT

a score=0
s ref.chr1 8 4 + 12 GGGG
s qry.chr1 8 4 + 12 GGGA
"""


def test_maf_projection(tmp_path: Path) -> None:
    maf = tmp_path / "a.maf"
    maf.write_text(_MAF)
    out = tmp_path / "out.fasta"
    maf_to_fasta(maf, "ref", out)
    records = _read_fasta(out)
    # reference first, species keyed on genome part of name
    assert list(records)[0] == "ref"
    assert records["ref"] == "ACGTACGTGGGG"
    assert records["qry"] == "ACGAACGTGGGA"
    assert len(records["ref"]) == len(records["qry"])


def test_maf_reference_gap_dropped(tmp_path: Path) -> None:
    maf = tmp_path / "a.maf"
    maf.write_text(
        "a\n"
        "s ref.c 0 6 + 6 AC--GT\n"
        "s qry.c 0 6 + 6 ACTTGT\n"
    )
    out = tmp_path / "out.fasta"
    maf_to_fasta(maf, "ref", out)
    records = _read_fasta(out)
    assert records["ref"] == "ACGT"
    assert records["qry"] == "ACGT"


def test_name_map_reference_with_version_dot(tmp_path: Path) -> None:
    # SibeliaZ-style: MAF source names are raw contig IDs; name_map groups them
    # to genome labels that contain version dots (e.g. ``..._GCF_1.1``). The
    # reference is passed as such a dotted label. Regression for a bug where the
    # reference key got version-stripped and matched no row -> empty MSA.
    maf = tmp_path / "a.maf"
    maf.write_text(
        "a\n"
        "s c1 0 8 + 8 ACGTACGT\n"
        "s d1 0 8 + 8 ACGAACGT\n"
    )
    out = tmp_path / "out.fasta"
    name_map = {"c1": "G_GCF_1.1", "d1": "H_GCF_2.1"}
    maf_to_fasta(maf, "G_GCF_1.1", out, name_map=name_map)
    records = _read_fasta(out)
    assert list(records)[0] == "G_GCF_1.1"
    assert records["G_GCF_1.1"] == "ACGTACGT"
    assert records["H_GCF_2.1"] == "ACGAACGT"


def test_exclude_drops_pseudo_genome(tmp_path: Path) -> None:
    # Minigraph-Cactus adds a _MINIGRAPH_ backbone; exclude must drop it as a taxon.
    maf = tmp_path / "a.maf"
    maf.write_text(
        "a\n"
        "s ref.chr1 0 4 + 4 ACGT\n"
        "s qry.chr1 0 4 + 4 ACGA\n"
        "s _MINIGRAPH_.chr1 0 4 + 4 ACGT\n"
    )
    out = tmp_path / "out.fasta"
    maf_to_fasta(maf, "ref", out, exclude={"_MINIGRAPH_"})
    records = _read_fasta(out)
    assert set(records) == {"ref", "qry"}
    assert "_MINIGRAPH_" not in records


def test_empty_projection_raises(tmp_path: Path) -> None:
    # Reference label that matches no row must fail loudly, not write an empty MSA.
    maf = tmp_path / "a.maf"
    maf.write_text("a\ns c1 0 4 + 4 ACGT\ns d1 0 4 + 4 ACGA\n")
    out = tmp_path / "out.fasta"
    with pytest.raises(ValueError, match="zero-length"):
        maf_to_fasta(maf, "absent", out, name_map={"c1": "G", "d1": "H"})


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            records[name] = ""
        elif name is not None:
            records[name] += line.strip()
    return records
