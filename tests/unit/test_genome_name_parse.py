"""Canonical genome-name round-trip: one parser for bacterial and viral names."""

from __future__ import annotations

from repgenr.core.contracts import (
    accession_from_filename,
    genome_filename,
    parse_genome_filename,
)


def test_bacterial_round_trip() -> None:
    name = genome_filename("Francisellaceae", "Francisella", "tularensis", "GCF_000017785.1")
    assert name == "Francisellaceae_Francisella_tularensis_GCF_000017785.1.fasta"
    assert parse_genome_filename(name) == (
        "Francisellaceae", "Francisella", "tularensis", "GCF_000017785.1",
    )
    assert accession_from_filename(name) == "GCF_000017785.1"


def test_viral_accession_with_underscore() -> None:
    name = genome_filename("Retroviridae", "Lentivirus", "HIV-1", "NC_001802.1")
    assert accession_from_filename(name) == "NC_001802.1"


def test_viral_accession_without_underscore() -> None:
    # GenBank viral accessions often have no underscore -- must still round-trip.
    name = genome_filename("Coronaviridae", "Betacoronavirus", "SARS-CoV-2", "MN908947.3")
    fam, gen, sp, acc = parse_genome_filename(name)
    assert (fam, gen, sp, acc) == ("Coronaviridae", "Betacoronavirus", "SARS-CoV-2", "MN908947.3")


def test_parses_leaf_stem_without_suffix() -> None:
    leaf = "Francisellaceae_Francisella_tularensis_GCF_000017785.1"
    assert accession_from_filename(leaf) == "GCF_000017785.1"


def test_non_canonical_name() -> None:
    fam, gen, sp, acc = parse_genome_filename("weird.fasta")
    assert (fam, gen, sp) == ("", "", "")
    assert acc == "weird"
