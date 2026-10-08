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


def _rec(
    acc, species, length, *, genus="Mastadenovirus", isolate="", organism=None, segment="ANONYMOUS"
):
    return VirusRecord(
        accession=acc,
        taxid="1",
        organism=organism or species,
        family="Adenoviridae",
        genus=genus,
        species=species,
        length=length,
        completeness="COMPLETE",
        segment=segment,
        isolate=isolate,
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
        _rec("seg1", "Influenza A", 1000, isolate="A/duck/2019", segment="4"),
        _rec("seg2", "Influenza A", 800, isolate="A/duck/2019", segment="6"),
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


# --- deep audit 2026-10 (entry-net) -------------------------------------------


def _seqs(records) -> dict:
    return {
        r.accession: SeqRecord(Seq("A" * r.length), id=r.accession, description=r.accession)
        for r in records
    }


def _group_rows(tmp_path, records):
    out = tmp_path / "g"
    out.mkdir()
    segments: dict[str, list[str]] = {}
    rows = _write_isolate_groups(out, records, _seqs(records), _LOG, segments)
    return out, rows, segments


def test_isolate_keeps_one_record_per_segment(tmp_path) -> None:
    """Live (Lassa 'Josiah'): three L and three S records of one isolate were
    concatenated into one 21 kb genome."""
    recs = [
        _rec("AY628202.1", "Lassa", 7285, isolate="Josiah", segment="L"),
        _rec("HQ688674.1", "Lassa", 7285, isolate="Josiah", segment="L"),
        _rec("NC_004297.1", "Lassa", 7279, isolate="Josiah", segment="L"),
        _rec("AY628203.1", "Lassa", 3402, isolate="Josiah", segment="S"),
        _rec("NC_004296.1", "Lassa", 3402, isolate="Josiah", segment="S"),
        _rec("HQ688672.1", "Lassa", 3401, isolate="Josiah", segment="S"),
    ]
    out, rows, segments = _group_rows(tmp_path, recs)
    assert [r.accession for r in rows] == ["iso-Josiah"]
    assert segments == {"iso-Josiah": ["AY628202.1", "AY628203.1"]}
    (genome,) = out.iterdir()
    seq = "".join(line for line in genome.read_text().splitlines() if not line.startswith(">"))
    assert len(seq) == 7285 + 3402


def test_isolate_name_shared_across_species_is_not_concatenated(tmp_path) -> None:
    """Live (Mammarenavirus): 'Acar 3080' named a Lassa S and a Mobala L."""
    recs = [
        _rec("AY628208.1", "Lassa", 3427, isolate="Acar 3080", segment="S"),
        _rec("DQ328876.1", "Mobala", 7325, isolate="Acar 3080", segment="L"),
        _rec("DQ328877.1", "Mobala", 3400, isolate="Acar 3080", segment="S"),
    ]
    _out, rows, segments = _group_rows(tmp_path, recs)
    assert sorted(r.accession for r in rows) == ["AY628208.1", "iso-Acar-3080"]
    assert segments == {"iso-Acar-3080": ["DQ328876.1", "DQ328877.1"]}


def test_isolate_tokens_that_collide_are_made_unique(tmp_path) -> None:
    """'Candid #1' and 'Candid-1' both sanitise to iso-Candid-1."""
    recs = [
        _rec("AY746353.1", "Junin", 3413, isolate="Candid-1", segment="S"),
        _rec("AY746354.1", "Junin", 7114, isolate="Candid-1", segment="L"),
        _rec("AY819707.2", "Junin", 7114, isolate="Candid #1", segment="L"),
        _rec("AY819708.1", "Junin", 3413, isolate="Candid #1", segment="S"),
    ]
    out, rows, _segments = _group_rows(tmp_path, recs)
    accessions = sorted(r.accession for r in rows)
    assert accessions == ["iso-Candid-1-AY746353.1", "iso-Candid-1-AY819707.2"]
    assert len(list(out.iterdir())) == 2


def test_unlabelled_record_of_a_segmented_isolate_stays_single(tmp_path) -> None:
    recs = [
        _rec("S1", "Lassa", 3400, isolate="X", segment="S"),
        _rec("L1", "Lassa", 7200, isolate="X", segment="L"),
        _rec("U1", "Lassa", 3000, isolate="X", segment="ANONYMOUS"),
    ]
    _out, rows, segments = _group_rows(tmp_path, recs)
    assert sorted(r.accession for r in rows) == ["U1", "iso-X"]
    assert segments == {"iso-X": ["L1", "S1"]}


def test_segment_label_variants_count_as_one_segment(tmp_path) -> None:
    """Live (Orthohantavirus): 'M', 'M; medium' and 'middle' are one segment,
    so an isolate keeps one M record and its grouped genome has three segments."""
    recs = [
        _rec("S.1", "Hantaan", 1700, isolate="76-118", segment="S; small"),
        _rec("M.1", "Hantaan", 3600, isolate="76-118", segment="M"),
        _rec("M.2", "Hantaan", 3610, isolate="76-118", segment="M; medium"),
        _rec("M.3", "Hantaan", 3500, isolate="76-118", segment="middle"),
        _rec("L.1", "Hantaan", 6500, isolate="76-118", segment="large"),
    ]
    out, rows, segments = _group_rows(tmp_path, recs)
    assert len(rows) == 1
    assert segments == {"iso-76-118": ["L.1", "M.2", "S.1"]}
    (genome,) = out.iterdir()
    assert "(3 segments)" in genome.read_text()


def test_numbered_rna_segments_are_not_merged(tmp_path) -> None:
    recs = [
        _rec("R1.1", "Tospo", 8900, isolate="X", segment="RNA 1"),
        _rec("R2.1", "Tospo", 4800, isolate="X", segment="RNA 2"),
        _rec("R3.1", "Tospo", 2900, isolate="X", segment="RNA 3"),
    ]
    _out, _rows, segments = _group_rows(tmp_path, recs)
    assert segments == {"iso-X": ["R1.1", "R2.1", "R3.1"]}


def test_unknown_segment_label_is_unlabelled(tmp_path) -> None:
    recs = [
        _rec("S1", "Phlebo", 1700, isolate="Y", segment="S"),
        _rec("U1", "Phlebo", 6400, isolate="Y", segment="Unknown"),
    ]
    _out, rows, segments = _group_rows(tmp_path, recs)
    assert segments == {} and len(rows) == 2
