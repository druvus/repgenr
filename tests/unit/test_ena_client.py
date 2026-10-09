"""ENA Portal client: query building, taxon resolution, paging, row normalisation."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from repgenr.core import ena, http
from repgenr.core.errors import UserInputError, WorkdirError

DATA = Path(__file__).parent / "data"


def _load(name: str):
    return json.loads((DATA / name).read_text(encoding="utf-8"))


def test_taxon_query_restricts_to_wgs_genomic() -> None:
    q = ena.taxon_query("2097")
    assert "tax_tree(2097)" in q
    assert 'library_strategy="WGS"' in q and 'library_source="GENOMIC"' in q


def test_accession_query_maps_each_accession_type_to_its_field() -> None:
    q = ena.accession_query(["SRR1", "ERR2", "PRJNA3", "SAMN4", "SRS5", "PRJEB6"])
    assert 'run_accession="SRR1"' in q and 'run_accession="ERR2"' in q
    assert 'study_accession="PRJNA3"' in q and 'study_accession="PRJEB6"' in q
    assert 'sample_accession="SAMN4"' in q
    assert 'secondary_sample_accession="SRS5"' in q
    assert " OR " in q


def test_accession_query_rejects_unknown_forms() -> None:
    with pytest.raises(UserInputError, match="GCF_000001"):
        ena.accession_query(["GCF_000001"])


def test_resolve_taxon_uses_any_name_and_returns_the_hit(monkeypatch) -> None:
    seen = {}

    def fake_get_json(url, params=None):
        seen["url"] = url
        return _load("ena_taxon_any_name.json")

    monkeypatch.setattr(http, "get_json", fake_get_json)
    hit = ena.resolve_taxon("Mycoplasma genitalium")
    assert seen["url"].endswith("/any-name/Mycoplasma%20genitalium")
    assert (hit.taxid, hit.scientific_name, hit.rank) == (
        "2097",
        "Mycoplasmoides genitalium",
        "species",
    )


def test_resolve_taxon_reports_ambiguity_and_absence(monkeypatch) -> None:
    monkeypatch.setattr(
        http,
        "get_json",
        lambda url, params=None: [
            {"taxId": "1", "scientificName": "A thing", "rank": "genus"},
            {"taxId": "2", "scientificName": "A thing", "rank": "species"},
        ],
    )
    with pytest.raises(UserInputError, match="taxid 1.*taxid 2|1 .* 2"):
        ena.resolve_taxon("A thing")
    monkeypatch.setattr(http, "get_json", lambda url, params=None: [])
    with pytest.raises(UserInputError, match="Nothing"):
        ena.resolve_taxon("Nonexistentus")


def test_search_runs_asks_for_every_record_in_one_request(monkeypatch) -> None:
    """The portal rejects offset (400 'Unsupported param offset'); limit=0
    returns every match, beyond the 10000 records a page held."""
    calls = []

    def fake_iter_lines(url, params=None):
        calls.append(dict(params))
        fields = params["fields"].split(",")
        yield "\t".join(fields)
        yield from ("\t".join([f"R{i}", *[""] * (len(fields) - 1)]) for i in range(4))

    monkeypatch.setattr(http, "iter_lines", fake_iter_lines)
    rows = ena.search_runs("tax_tree(1)")
    assert [r["run_accession"] for r in rows] == ["R0", "R1", "R2", "R3"]
    assert len(calls) == 1 and calls[0]["limit"] == 0 and "offset" not in calls[0]
    assert calls[0]["result"] == "read_run" and calls[0]["format"] == "tsv"
    assert calls[0]["fields"] == ",".join(ena.RUN_FIELDS)


def test_search_runs_requests_the_named_fields(monkeypatch) -> None:
    calls = []

    def fake_iter_lines(url, params=None):
        calls.append(params["fields"])
        return iter(["run_accession\ttax_id", "R1\t9"])

    monkeypatch.setattr(http, "iter_lines", fake_iter_lines)
    assert ena.search_runs("q", fields=("run_accession", "tax_id")) == [
        {"run_accession": "R1", "tax_id": "9"}
    ]
    assert calls == ["run_accession,tax_id"]


def test_tsv_records_have_the_shape_of_the_json_records() -> None:
    """The TSV reader returns the dicts the JSON format did: every value a
    string, an empty string where the portal has none, so to_read_rows and
    the census read them unchanged."""
    json_records = _load("ena_read_run_mixed.json")
    fields = list(json_records[0].keys())
    lines = ["\t".join(fields)]
    lines += ["\t".join(str(r.get(f, "")) for f in fields) for r in json_records]
    tsv_records = ena.read_tsv_records(iter(["", *lines, ""]))
    assert tsv_records == [{f: str(r.get(f, "")) for f in fields} for r in json_records]
    assert ena.to_read_rows(tsv_records) == ena.to_read_rows(json_records)


def test_tsv_reader_rejects_a_record_cut_short() -> None:
    assert ena.read_tsv_records(iter([])) == []
    assert ena.read_tsv_records(iter(["run_accession\ttax_id"])) == []
    with pytest.raises(WorkdirError, match="1 columns where the header has 2"):
        ena.read_tsv_records(iter(["run_accession\ttax_id", "R1\t9", "R2"]))


def test_tsv_reader_refuses_a_header_other_than_the_fields() -> None:
    """An HTML page served with status 200 is not read as records."""
    page = iter(["<!DOCTYPE html>", "<html><body>Service unavailable</body></html>"])
    with pytest.raises(WorkdirError, match="Unexpected answer.*DOCTYPE"):
        ena.read_tsv_records(page, fields=("run_accession", "tax_id"))
    with pytest.raises(WorkdirError, match="Unexpected answer"):
        ena.read_tsv_records(iter(["tax_id\trun_accession"]), fields=("run_accession", "tax_id"))
    assert ena.read_tsv_records(
        iter(["run_accession\ttax_id", "R1\t9"]), fields=("run_accession", "tax_id")
    ) == [{"run_accession": "R1", "tax_id": "9"}]


def test_search_runs_refuses_an_html_answer(monkeypatch) -> None:
    monkeypatch.setattr(
        http, "iter_lines", lambda url, params=None: iter(["<html>", "<p>maintenance</p>"])
    )
    with pytest.raises(WorkdirError, match="Unexpected answer"):
        ena.search_runs("q")


def test_count_runs_reads_the_count_endpoint(monkeypatch) -> None:
    calls = []

    def fake_get_json(url, params=None):
        calls.append((url, dict(params)))
        return {"count": "217249"}

    monkeypatch.setattr(http, "get_json", fake_get_json)
    assert ena.count_runs("tax_tree(1279)") == 217249
    url, params = calls[0]
    assert url == ena.COUNT_URL and url.endswith("/portal/api/count")
    assert params == {"result": "read_run", "query": "tax_tree(1279)", "format": "json"}
    monkeypatch.setattr(http, "get_json", lambda url, params=None: {"message": "no"})
    with pytest.raises(WorkdirError, match="Unexpected answer"):
        ena.count_runs("q")


def _warnings(caplog) -> list[str]:
    return [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]


def test_announce_search_warns_only_above_the_threshold(monkeypatch, caplog) -> None:
    log = logging.getLogger("repgenr.test")
    monkeypatch.setattr(ena, "count_runs", lambda query: ena.LARGE_SEARCH)
    with caplog.at_level(logging.INFO):
        assert ena.announce_search("q", "Staphylococcus", log) == ena.LARGE_SEARCH
    assert _warnings(caplog) == []
    assert any(f"ENA counts {ena.LARGE_SEARCH} " in r.getMessage() for r in caplog.records)

    caplog.clear()
    monkeypatch.setattr(ena, "count_runs", lambda query: 878821)
    with caplog.at_level(logging.INFO):
        assert ena.announce_search("q", "Salmonella", log) == 878821
    (warning,) = _warnings(caplog)
    assert "Salmonella has 878821 ENA runs (more than 100000)" in warning
    assert "GB of memory" in warning
    # Fewer fields, a smaller expected size.
    caplog.clear()
    with caplog.at_level(logging.INFO):
        ena.announce_search("q", "Salmonella", log, fields=("run_accession",))
    (small,) = _warnings(caplog)
    assert small != warning


def test_announce_search_continues_when_the_count_fails(monkeypatch, caplog) -> None:
    def down(query):
        raise WorkdirError("HTTP request failed: count")

    monkeypatch.setattr(ena, "count_runs", down)
    with caplog.at_level(logging.INFO):
        assert ena.announce_search("q", "Staphylococcus", logging.getLogger("t")) is None
    assert "Could not count" in _warnings(caplog)[0]


def test_to_read_rows_normalises_the_portal_records() -> None:
    rows = ena.to_read_rows(_load("ena_read_run_mixed.json"))
    by_run = {r.run_accession: r for r in rows}
    ont = by_run["SRR28800588"]
    assert (ont.platform, ont.layout, ont.taxid) == ("OXFORD_NANOPORE", "SINGLE", "2097")
    assert ont.organism == "Mycoplasmoides genitalium"
    assert ont.bases == 90992597 and len(ont.fastq_urls) == 1
    assert ont.fastq_urls[0].startswith("https://ftp.sra.ebi.ac.uk/")
    assert ont.fastq_bytes == (85028322,) and len(ont.fastq_md5[0]) == 32
    no_mirror = by_run["ERR2237850"]
    assert no_mirror.fastq_urls == () and no_mirror.fastq_bytes == ()
    assert no_mirror.platform == "PACBIO_SMRT"
    assert ont.library_selection == "RANDOM"
    taxon_rows = ena.to_read_rows(_load("ena_read_run_taxon.json"))
    paired = taxon_rows[0]
    assert len(paired.fastq_urls) == 2 and paired.layout == "PAIRED"
    assert {r.run_accession for r in taxon_rows if r.library_selection == "MDA"} == {"ERR17019821"}
    assert "library_selection" in ena.RUN_FIELDS
