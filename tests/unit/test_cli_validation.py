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


def test_phylo_build_mask_needs_the_snptype_source(tmp_path) -> None:
    """--mask on the stateless step is refused under the aligner source (D-6)."""
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    genomes = tmp_path / "g"
    genomes.mkdir()
    result = CliRunner().invoke(
        app,
        [
            "phylo-build",
            "--genomes-dir",
            str(genomes),
            "-o",
            str(tmp_path / "o"),
            "--mask",
            "gubbins",
        ],
    )
    assert result.exit_code != 0
    assert "--msa-source snptype" in result.output + str(result.exception or "")
