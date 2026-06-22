"""Offline tests for viral parsing helpers (no network)."""

from __future__ import annotations

from pathlib import Path

from repgenr.viral._common import select_outgroup_from_matrix
from repgenr.viral.bvbrc import _matches, _read_ncbi
from repgenr.viral.entrez import TAXNAMES_ORDERED, _group_taxa


def test_group_taxa_splits_on_top_level_taxon() -> None:
    lines = [
        "<Taxon>",
        "  <TaxId>11620</TaxId>",
        "  <ScientificName>Lassa</ScientificName>",
        "</Taxon>",
        "<Taxon>",
        "  <TaxId>10535</TaxId>",
        "  <ScientificName>Adeno</ScientificName>",
        "</Taxon>",
    ]
    groups = _group_taxa(lines)
    assert len(groups) == 2
    assert "11620" in groups[0]
    assert "10535" in groups[1]


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
