"""Unit tests for CLI input validators."""

from __future__ import annotations

import pytest

from repgenr.cli.base import _require_choice, _require_unit_interval
from repgenr.core.errors import UserInputError


def test_require_choice_accepts_valid() -> None:
    _require_choice("skder", {"skder", "galah"}, "--tool")  # no raise


def test_require_choice_rejects_invalid() -> None:
    with pytest.raises(UserInputError, match=r"Invalid --tool 'bogus'.*galah, skder"):
        _require_choice("bogus", {"skder", "galah"}, "--tool")


@pytest.mark.parametrize("value", [0.5, 1.0, 0.99995, None])
def test_unit_interval_accepts(value) -> None:
    _require_unit_interval(value, "--secondary-ani")  # no raise


@pytest.mark.parametrize("value", [0.0, -0.1, 1.5, 2.0])
def test_unit_interval_rejects(value) -> None:
    with pytest.raises(UserInputError, match="must be in"):
        _require_unit_interval(value, "--secondary-ani")
