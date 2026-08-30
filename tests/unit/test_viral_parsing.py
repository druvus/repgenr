"""Offline tests for viral parsing helpers (no network)."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.viral._common import select_outgroup_from_matrix
from repgenr.viral.bvbrc import _matches, _read_ncbi
from repgenr.viral.entrez import TAXNAMES_ORDERED, _iter_taxa


def test_iter_taxa_returns_top_level_taxon_elements() -> None:
    xml = (
        "<TaxaSet>"
        "<Taxon><TaxId>11620</TaxId><ScientificName>Lassa</ScientificName></Taxon>"
        "<Taxon><TaxId>10535</TaxId><ScientificName>Adeno</ScientificName></Taxon>"
        "</TaxaSet>"
    )
    taxa = _iter_taxa(xml)
    assert [t.findtext("TaxId") for t in taxa] == ["11620", "10535"]


def test_read_ncbi_and_matches(tmp_path: Path) -> None:
    n = len(TAXNAMES_ORDERED)
    # build a metadata_ncbi.tsv row: taxid, name, num_with_tag, <n names>, <n taxids>
    names = ["" for _ in range(n)]
    taxids = ["" for _ in range(n)]
    gi = TAXNAMES_ORDERED.index("genus")
    si = TAXNAMES_ORDERED.index("species")
    names[gi], taxids[gi] = "Mastadenovirus", "10509"
    names[si], taxids[si] = "Human mastadenovirus C", "129951"

    taxid_cols = [f"{x}_taxid" for x in TAXNAMES_ORDERED]
    header = ["taxid", "name", "num_with_tag", *TAXNAMES_ORDERED, *taxid_cols]
    row = ["10535", "Adeno", "5", *names, *taxids]
    path = tmp_path / "metadata_ncbi.tsv"
    path.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n")

    ncbi = _read_ncbi(path)
    assert "10535" in ncbi
    # match by name (case-insensitive) and by taxid
    assert _matches(ncbi, "10535", "genus", ["mastadenovirus"])
    assert _matches(ncbi, "10535", "genus", ["10509"])
    assert not _matches(ncbi, "10535", "genus", ["francisella"])
    assert _matches(ncbi, "10535", "species", ["human mastadenovirus c"])


def test_select_outgroup_prefers_distant_candidate(tmp_path: Path) -> None:
    # header: row-id, then S/O sequence columns; rows mirror columns
    matrix = tmp_path / "dist.tsv"
    matrix.write_text(
        "\tS_a\tS_b\tO_x\tO_y\n"
        "S_a\t0\t0.02\t0.30\t0.50\n"
        "S_b\t0.02\t0\t0.31\t0.51\n"
        "O_x\t0.30\t0.31\t0\t0.40\n"
        "O_y\t0.50\t0.51\t0.40\t0\n"
    )
    chosen = select_outgroup_from_matrix(matrix, logger=_NullLogger())
    assert chosen in ("O_x", "O_y")


class _NullLogger:
    def info(self, *a, **k):
        pass

    def warning(self, *a, **k):
        pass


def test_parse_custom_filter_splits_at_first_colon() -> None:
    from repgenr.viral._common import parse_custom_filter

    assert parse_custom_filter("host:homo sapiens") == ("host", "homo sapiens")
    # values may themselves contain colons (e.g. strain names)
    assert parse_custom_filter("strain:a/wuhan:2019") == ("strain", "a/wuhan:2019")


@pytest.mark.parametrize("bad", ["nocolon", ":value", "key:", ":"])
def test_parse_custom_filter_rejects_malformed(bad: str) -> None:
    from repgenr.core.errors import UserInputError
    from repgenr.viral._common import parse_custom_filter

    with pytest.raises(UserInputError, match="target-custom"):
        parse_custom_filter(bad)


def test_bvbrc_custom_filter_value_with_colon(tmp_path: Path) -> None:
    # A colon inside the value crashed the old split(":") unpacking.
    import logging

    from repgenr.viral.bvbrc import _Record, _select_by_taxonomy

    records = [_Record(name="r1", bvbrc_id="b1", taxid="11", description="d", length=100)]
    ncbi = {"11": [{"taxlevelname": "strain", "taxname": "a/wuhan:2019", "taxid": "11"}]}
    selected, headers = _select_by_taxonomy(
        records, ncbi, {"custom": ["strain:a/wuhan:2019"]}, logging.getLogger("test")
    )
    assert "11" in selected and selected["11"] == {"b1": 100}
