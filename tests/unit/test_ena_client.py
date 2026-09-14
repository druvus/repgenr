"""ENA Portal client: query building, taxon resolution, paging, row normalisation."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from repgenr.core import ena, http
from repgenr.core.errors import UserInputError

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


def test_search_runs_pages_until_a_short_page(monkeypatch) -> None:
    calls = []

    def fake_get_json(url, params=None):
        calls.append(dict(params))
        offset = int(params["offset"])
        if offset == 0:
            return [{"run_accession": f"R{i}"} for i in range(3)]
        return [{"run_accession": "R3"}]

    monkeypatch.setattr(http, "get_json", fake_get_json)
    rows = ena.search_runs("tax_tree(1)", page=3)
    assert [r["run_accession"] for r in rows] == ["R0", "R1", "R2", "R3"]
    assert [c["offset"] for c in calls] == [0, 3]
    assert calls[0]["result"] == "read_run" and calls[0]["format"] == "json"


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
    paired = ena.to_read_rows(_load("ena_read_run_taxon.json"))[0]
    assert len(paired.fastq_urls) == 2 and paired.layout == "PAIRED"
