"""Converters: reverse_complement, XMFA->FASTA happy path, GFA missing-tool."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.converters import gfa_to_fasta as gfa_mod
from repgenr.converters.gfa_to_fasta import gfa_to_fasta
from repgenr.converters.xmfa_to_fasta import reverse_complement, xmfa_to_fasta
from repgenr.core.errors import MissingBinaryError

_LOG = logging.getLogger("test")


def test_reverse_complement() -> None:
    assert reverse_complement(b"AAATTTGC") == b"GCAAATTT"
    assert reverse_complement(reverse_complement(b"ACGTRYMK")) == b"ACGTRYMK"
    assert reverse_complement(b"ACGT") == b"ACGT"  # palindrome


def _read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            out[name] = ""
        elif name is not None:
            out[name] += line.strip()
    return out


def test_xmfa_to_fasta_projects_two_sequences(tmp_path: Path) -> None:
    xmfa = tmp_path / "aln.xmfa"
    xmfa.write_text(
        "#FormatVersion Mauve1\n"
        "#Sequence1File\tref.fa\n#Sequence1Format\tFastA\n"
        "#Sequence2File\tqry.fa\n#Sequence2Format\tFastA\n"
        "> 1:1-8 + ref.fa\nACGTACGT\n"
        "> 2:1-8 + qry.fa\nACGTACGT\n"
        "=\n"
    )
    out = xmfa_to_fasta(xmfa, "ref.fa", 0, tmp_path / "msa.fasta")
    recs = _read_fasta(out)
    assert set(recs) == {"ref.fa", "qry.fa"}
    assert recs["ref.fa"] == "ACGTACGT"
    assert recs["qry.fa"] == "ACGTACGT"


def test_xmfa_reverse_strand_block(tmp_path: Path) -> None:
    # query aligned on the minus strand -> stored reverse-complemented onto ref coords
    xmfa = tmp_path / "rc.xmfa"
    xmfa.write_text(
        "#Sequence1File\tref.fa\n#Sequence2File\tqry.fa\n"
        "> 1:1-4 + ref.fa\nACGT\n"
        "> 2:1-4 - qry.fa\nACGT\n"
        "=\n"
    )
    out = xmfa_to_fasta(xmfa, "ref.fa", 0, tmp_path / "rc.fasta")
    recs = _read_fasta(out)
    assert recs["ref.fa"] == "ACGT"
    assert recs["qry.fa"] == "ACGT"  # ACGT is its own reverse complement


def test_gfa_requires_odgi(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(gfa_mod.shutil, "which", lambda _n: None)
    with pytest.raises(MissingBinaryError, match="odgi"):
        gfa_to_fasta(tmp_path / "g.gfa", tmp_path / "out.fasta", _LOG)
