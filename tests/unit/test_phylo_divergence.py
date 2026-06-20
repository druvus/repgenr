"""Divergence-aware routing: taxonomic spread + aligner warnings, and key=value parsing."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.cli.main import _parse_key_values
from repgenr.core.errors import UserInputError
from repgenr.stages.phylo import _taxonomic_spread, _warn_divergence


def _g(names: list[str]) -> list[Path]:
    return [Path(n) for n in names]


def test_taxonomic_spread_strains() -> None:
    # one species, many strains -> 1 genus, 1 species
    g = _g([
        "Francisellaceae_Francisella_tularensis_GCA_1.fasta",
        "Francisellaceae_Francisella_tularensis_GCA_2.fasta",
    ])
    assert _taxonomic_spread(g) == (1, 1)


def test_taxonomic_spread_genus() -> None:
    g = _g([
        "Francisellaceae_Francisella_tularensis_GCA_1.fasta",
        "Francisellaceae_Francisella_philomiragia_GCA_2.fasta",
    ])
    assert _taxonomic_spread(g) == (1, 2)  # 1 genus, 2 species


def test_taxonomic_spread_family() -> None:
    g = _g([
        "Francisellaceae_Francisella_tularensis_GCA_1.fasta",
        "Francisellaceae_Caedibacter_halobius_GCA_2.fasta",
    ])
    assert _taxonomic_spread(g) == (2, 2)  # 2 genera


def test_warn_family_level(caplog) -> None:
    g = _g([
        "Fam_Francisella_tularensis_GCA_1.fasta",
        "Fam_Caedibacter_halobius_GCA_2.fasta",
    ])
    with caplog.at_level(logging.WARNING):
        _warn_divergence("sibeliaz", g, logging.getLogger("t"))
    assert any("family-level" in r.message for r in caplog.records)


def test_warn_cactus_intraspecific(caplog) -> None:
    g = _g([
        "Fam_Francisella_tularensis_GCA_1.fasta",
        "Fam_Francisella_philomiragia_GCA_2.fasta",
    ])
    with caplog.at_level(logging.WARNING):
        _warn_divergence("cactus", g, logging.getLogger("t"))
    assert any("Minigraph-Cactus" in r.message for r in caplog.records)


def test_no_warn_for_strains(caplog) -> None:
    g = _g([
        "Fam_Francisella_tularensis_GCA_1.fasta",
        "Fam_Francisella_tularensis_GCA_2.fasta",
    ])
    with caplog.at_level(logging.WARNING):
        _warn_divergence("sibeliaz", g, logging.getLogger("t"))
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]


def test_parse_key_values() -> None:
    assert _parse_key_values(["kmer=15", "seed_weight=11"], "--aligner-arg") == {
        "kmer": "15", "seed_weight": "11",
    }
    assert _parse_key_values([], "--aligner-arg") == {}
    with pytest.raises(UserInputError):
        _parse_key_values(["bad"], "--aligner-arg")
