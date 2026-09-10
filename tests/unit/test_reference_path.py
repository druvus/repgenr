"""snptype --reference resolution must not traverse outside the genome dirs."""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from repgenr.core.errors import UserInputError
from repgenr.stages.snptype import _reference_path


def _ctx(tmp_path: Path) -> SimpleNamespace:
    reps = tmp_path / "derep" / "representatives"
    genomes = tmp_path / "genomes"
    reps.mkdir(parents=True)
    genomes.mkdir()
    (genomes / "g1.fasta").write_text(">x\nACGT\n")
    return SimpleNamespace(representatives_dir=reps, genomes_dir=genomes)


_LOG = logging.getLogger("test")


def test_reference_resolves_by_basename(tmp_path: Path) -> None:
    ctx = _ctx(tmp_path)
    assert _reference_path(ctx, "g1.fasta", [], _LOG) == ctx.genomes_dir / "g1.fasta"


def test_reference_with_separator_rejected(tmp_path: Path) -> None:
    ctx = _ctx(tmp_path)
    (tmp_path / "secret.txt").write_text("x")
    with pytest.raises(UserInputError, match="basename"):
        _reference_path(ctx, "../secret.txt", [], _LOG)
    with pytest.raises(UserInputError, match="basename"):
        _reference_path(ctx, str(tmp_path / "secret.txt"), [], _LOG)


def test_default_reference_is_reported(tmp_path, caplog) -> None:
    """Without --reference the stage says which genome it picked (live audit)."""
    from pathlib import Path

    from repgenr.core.context import WorkdirContext

    ctx = WorkdirContext(tmp_path, create=True)
    genomes = [Path("b.fasta"), Path("c.fasta")]
    with caplog.at_level(logging.WARNING):
        assert _reference_path(ctx, None, genomes, logging.getLogger("t")) == Path("b.fasta")
    assert "No --reference given" in caplog.text and "b.fasta" in caplog.text
