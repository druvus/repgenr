"""Offline tests for the Entrez taxonomy parser (no network).

The HTTP call is replaced with a fixture XML payload so the lineage parsing
(get_taxon_data_from_entrez / _parse_taxon) is exercised without touching NCBI.
"""

from __future__ import annotations

import logging

from repgenr.viral import entrez
from repgenr.viral.entrez import _parse_taxon, get_taxon_data_from_entrez

_LOG = logging.getLogger("test")

# One Taxon with a LineageEx, in the line-indentation shape _group_taxa expects:
# the opening <Taxon> is unindented, everything inside is indented.
_XML = """<?xml version="1.0"?>
<TaxaSet>
<Taxon>
    <TaxId>10535</TaxId>
    <ScientificName>Human mastadenovirus C</ScientificName>
    <Rank>species</Rank>
    <LineageEx>
        <Taxon>
            <TaxId>10508</TaxId>
            <ScientificName>Adenoviridae</ScientificName>
            <Rank>family</Rank>
        </Taxon>
        <Taxon>
            <TaxId>10509</TaxId>
            <ScientificName>Mastadenovirus</ScientificName>
            <Rank>genus</Rank>
        </Taxon>
    </LineageEx>
</Taxon>
</TaxaSet>
"""


def test_get_taxon_data_parses_lineage(monkeypatch) -> None:
    monkeypatch.setattr(entrez, "_send_query", lambda ids: _XML.split("\n"))
    monkeypatch.setattr(entrez, "sleep", lambda *_a, **_k: None)  # no rate-limit waits

    data, missing, _alts = get_taxon_data_from_entrez(["10535"], _LOG)

    assert not missing
    taxdata = data["10535"]["taxdata"]
    assert data["10535"]["name"] == "Human mastadenovirus C"
    assert taxdata["genus"]["name"] == "Mastadenovirus" and taxdata["genus"]["taxid"] == "10509"
    assert taxdata["family"]["name"] == "Adenoviridae" and taxdata["family"]["taxid"] == "10508"
    assert taxdata["species"]["name"] == "Human mastadenovirus C"


def test_get_taxon_data_marks_missing(monkeypatch) -> None:
    # Empty response -> the requested taxid is reported missing and gets an
    # all-None lineage placeholder rather than raising.
    monkeypatch.setattr(entrez, "_send_query", lambda ids: [""])
    monkeypatch.setattr(entrez, "sleep", lambda *_a, **_k: None)
    data, missing, _alts = get_taxon_data_from_entrez(["55555"], _LOG)
    assert missing == {"55555"}
    assert data["55555"]["taxdata"]["genus"]["name"] is None


def test_parse_taxon_ignores_unexpected_taxid() -> None:
    chunk = (
        "    <TaxId>999</TaxId>"
        "    <ScientificName>Ghost</ScientificName>"
        "    <Rank>species</Rank>"
    )
    # 999 was not requested, so the chunk is ignored (returns None).
    assert _parse_taxon(chunk, ["10535"], ["10535"], {}, _LOG) is None
