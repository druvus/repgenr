"""STAGE_DELIVERABLES: the outputs whose absence makes a completed stage rerun."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from repgenr.cli import base as cli


def _layout(workdir: Path) -> SimpleNamespace:
    return SimpleNamespace(workdir=workdir, derep_dir=workdir / "derep")


def test_every_stage_with_inputs_declares_deliverables() -> None:
    """Each workdir stage lists the outputs whose absence forces a rerun."""
    assert set(cli.STAGE_INPUTS) <= set(cli.STAGE_DELIVERABLES)


def test_derep_unpack_checks_unpacked_dir_with_representatives(tmp_path: Path) -> None:
    ctx = _layout(tmp_path)
    params = SimpleNamespace(no_representant=False)
    assert cli.missing_deliverables(ctx, "derep_unpack", params) == [tmp_path / "derep/unpacked"]
    (tmp_path / "derep" / "unpacked").mkdir(parents=True)
    # An empty unpacked/ still counts as missing when representatives are unpacked.
    assert cli.missing_deliverables(ctx, "derep_unpack", params) == [tmp_path / "derep/unpacked"]


def test_derep_unpack_without_representant_accepts_empty_unpacked_dir(tmp_path: Path) -> None:
    """--no-representant with only singleton clusters leaves unpacked/ empty;
    that must not rerun the stage on every call."""
    (tmp_path / "derep" / "unpacked").mkdir(parents=True)
    params = SimpleNamespace(no_representant=True)
    assert cli.missing_deliverables(_layout(tmp_path), "derep_unpack", params) == []
