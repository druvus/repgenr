"""Unit tests for the NCBI Virus selection helpers (no network, no mashtree)."""

from __future__ import annotations

import logging

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from repgenr.stages.vgenome import VgenomeParams
from repgenr.viral.ncbi_virus import VirusRecord
from repgenr.viral.selection import (
    _isolate_token,
    _length_range_records,
    _record_matches,
    _write_isolate_groups,
)

_LOG = logging.getLogger("test")


def _rec(acc, species, length, *, genus="Mastadenovirus", isolate="", organism=None):
    return VirusRecord(
        accession=acc, taxid="1", organism=organism or species, family="Adenoviridae",
        genus=genus, species=species, length=length, completeness="COMPLETE",
        segment="ANONYMOUS", isolate=isolate,
    )


def test_record_matches_levels() -> None:
    r = _rec("a", "Human mastadenovirus C", 300, organism="Human mastadenovirus C strain X")
    assert _record_matches(r, {"genus": ["mastadenovirus"]})
    assert _record_matches(r, {"species": ["human mastadenovirus c"]})
    assert not _record_matches(r, {"genus": ["lentivirus"]})
    # serotype matches a substring of the organism name
    assert _record_matches(r, {"serotype": ["strain x"]})
    # all requested levels must hold
    assert not _record_matches(r, {"genus": ["mastadenovirus"], "species": ["nope"]})


def test_length_range_records() -> None:
    recs = [_rec("a", "sp1", 300), _rec("b", "sp1", 310), _rec("c", "sp2", 900)]
    params = VgenomeParams(length_range="250-350")
    assert _length_range_records(recs, params, _LOG) == (250, 350)
    # without an explicit range: midpoint from per-species medians +/- deviation
    params = VgenomeParams(length_deviation=10, length_method="median_of_medians")
    lo, hi = _length_range_records(recs, params, _LOG)
    assert lo < hi and lo > 0


def test_isolate_token_sanitises() -> None:
    assert _isolate_token("A/duck/2019") == "iso-Aduck2019"
    assert _isolate_token("") == "iso-NA"


def test_write_isolate_groups(tmp_path, monkeypatch) -> None:
    # two segments of one isolate + one standalone record
    recs = [
        _rec("seg1", "Influenza A", 1000, isolate="A/duck/2019"),
        _rec("seg2", "Influenza A", 800, isolate="A/duck/2019"),
        _rec("solo", "Influenza A", 1300, isolate=""),
    ]
    seqs = {
        r.accession: SeqRecord(Seq("ACGT" * 10), id=r.accession, description=f"{r.accession} d")
        for r in recs
    }
    rows = _write_isolate_groups(tmp_path, recs, seqs, _LOG)
    # one grouped isolate genome + one singleton = 2 selection rows
    assert len(rows) == 2
    written = sorted(p.name for p in tmp_path.iterdir())
    assert any("iso-Aduck2019" in name for name in written)
    assert any(name.endswith("solo.fasta") for name in written)
    # the grouped genome concatenates both segments (longest first)
    grouped = next(p for p in tmp_path.iterdir() if "iso-Aduck2019" in p.name)
    assert "(2 segments)" in grouped.read_text()
