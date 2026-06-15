"""Offline tests for the GTDB-API metadata helpers (no network)."""

from __future__ import annotations

import pytest

from repgenr.core.errors import UserInputError
from repgenr.stages.metadata import (
    MetadataParams,
    _normalize_api_tax,
    _target_taxon,
    _validate,
)


def test_target_taxon_strings() -> None:
    p = MetadataParams(dataset="rep", level="genus", source="api", target_genus="francisella")
    assert _target_taxon(p) == "g__Francisella"

    p = MetadataParams(dataset="rep", level="family", source="api", target_family="francisellaceae")
    assert _target_taxon(p) == "f__Francisellaceae"

    p = MetadataParams(
        dataset="rep", level="species", source="api",
        target_genus="francisella", target_species="tularensis",
    )
    assert _target_taxon(p) == "s__Francisella tularensis"


def test_normalize_api_tax() -> None:
    row = {
        "gtdbFamily": "f__Francisellaceae",
        "gtdbGenus": "g__Francisella",
        "gtdbSpecies": "s__Francisella tularensis",
    }
    tax = _normalize_api_tax(row)
    assert tax["family"] == "Francisellaceae"
    assert tax["genus"] == "Francisella"
    # genus removed from species, spaces dropped
    assert tax["species"] == "tularensis"


def test_validate_api_does_not_require_release() -> None:
    # api source: release/version optional
    _validate(
        MetadataParams(dataset="rep", level="genus", source="api", target_genus="francisella")
    )


def test_validate_tsv_requires_release() -> None:
    with pytest.raises(UserInputError):
        _validate(MetadataParams(dataset="rep", level="genus", source="tsv", target_genus="x"))
