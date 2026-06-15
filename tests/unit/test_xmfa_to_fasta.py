"""Tests for the Python 3 XMFA -> FASTA port."""

from __future__ import annotations

from pathlib import Path

from repgenr.converters.xmfa_to_fasta import reverse_complement, xmfa_to_fasta

_XMFA = """#FormatVersion Mauve1
#Sequence1File\tgenomeA.fasta
#Sequence1Format\tFastA
#Sequence2File\tgenomeB.fasta
#Sequence2Format\tFastA
> 1:1-10 + genomeA.fasta
ACGTACGTAC
> 2:1-10 + genomeB.fasta
ACGTACGTAC
=
"""


def test_reverse_complement() -> None:
    assert reverse_complement(b"ACGT") == b"ACGT"
    assert reverse_complement(b"AAAA") == b"TTTT"
    assert reverse_complement(b"ACGTN".replace(b"N", b"")) == b"ACGT"


def test_basic_conversion(tmp_path: Path) -> None:
    xmfa = tmp_path / "aln.xmfa"
    xmfa.write_text(_XMFA)
    out = tmp_path / "out.fasta"
    xmfa_to_fasta(xmfa, "genomeA.fasta", 0, out)

    records = _read_fasta(out)
    # reference is written first
    assert list(records)[0] == "genomeA.fasta"
    assert records["genomeA.fasta"] == "ACGTACGTAC"
    assert records["genomeB.fasta"] == "ACGTACGTAC"
    # all sequences share the reference length
    assert len({len(s) for s in records.values()}) == 1


def test_reference_gap_removed(tmp_path: Path) -> None:
    # A gap in the reference column must be dropped from all sequences.
    xmfa = tmp_path / "aln.xmfa"
    xmfa.write_text(
        "#Sequence1File\tref.fasta\n"
        "#Sequence1Format\tFastA\n"
        "#Sequence2File\tqry.fasta\n"
        "#Sequence2Format\tFastA\n"
        "> 1:1-8 + ref.fasta\n"
        "ACGT--ACGT\n"
        "> 2:1-10 + qry.fasta\n"
        "ACGTTTACGT\n"
        "=\n"
    )
    out = tmp_path / "out.fasta"
    xmfa_to_fasta(xmfa, "ref.fasta", 0, out)
    records = _read_fasta(out)
    # reference has no gap characters once its own gaps are removed
    assert "-" not in records["ref.fasta"]
    assert records["ref.fasta"] == "ACGTACGT"
    # query loses the columns aligned to reference gaps
    assert records["qry.fasta"] == "ACGTACGT"


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
