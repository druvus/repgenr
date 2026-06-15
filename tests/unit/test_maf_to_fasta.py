"""Tests for MAF -> reference-anchored MSA-FASTA."""

from __future__ import annotations

from pathlib import Path

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
