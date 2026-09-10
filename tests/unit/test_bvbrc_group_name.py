"""BV-BRC group names resolve against the server listing, case-insensitively."""

from __future__ import annotations

from repgenr.stages.vmetadata import resolve_group_name

LISTING = ["Adenoviridae.fna", "Hepatitis_E_virus.fna", "Picornaviridae.fna", "README"]


def test_lower_case_target_matches_mixed_case_group() -> None:
    assert resolve_group_name("hepatitis_e_virus", LISTING) == "Hepatitis_E_virus"


def test_exact_name_still_matches() -> None:
    assert resolve_group_name("Adenoviridae", LISTING) == "Adenoviridae"


def test_unknown_group_is_none() -> None:
    assert resolve_group_name("hepatovirus", LISTING) is None
