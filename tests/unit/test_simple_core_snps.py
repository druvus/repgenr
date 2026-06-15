"""Tests for the simple SNP typer's core-SNP reduction."""

from __future__ import annotations

from pathlib import Path

from repgenr.snptypers.simple import _write_core_snps


def test_core_snp_reduction(tmp_path: Path) -> None:
    consensuses = {
        "ref": "ACGTACGT",
        "s1": "ACGAACGT",  # differs at col 3
        "s2": "ACGTACGA",  # differs at col 7
    }
    core = tmp_path / "core.fasta"
    matrix = tmp_path / "dist.tsv"
    n = _write_core_snps(consensuses, core, matrix)
    assert n == 2  # columns 3 and 7 are variable

    records = _read_fasta(core)
    assert records["ref"] == "TT"   # ref bases at the two variable columns
    assert records["s1"] == "AT"
    assert records["s2"] == "TA"

    # distance matrix: ref vs s1 = 1, ref vs s2 = 1, s1 vs s2 = 2
    lines = matrix.read_text().splitlines()
    assert lines[0].split("\t")[1:] == ["ref", "s1", "s2"]


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
