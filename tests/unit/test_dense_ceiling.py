"""sourmash dense path refuses an oversized N x N matrix (memory guard)."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import WorkdirError
from repgenr.dereplicators import sourmash
from repgenr.dereplicators.sourmash import SourmashDereplicator

_LOG = logging.getLogger("test")


def test_dense_refuses_above_ceiling(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(sourmash, "_DENSE_MAX_GENOMES", 2)
    # never reach the tool: the size guard fires first
    monkeypatch.setattr(
        sourmash, "run_tool", lambda *a, **k: pytest.fail("should not sketch")
    )
    genomes = [tmp_path / f"g{i}.fasta" for i in range(3)]
    with pytest.raises(WorkdirError, match="dense compare needs an N x N matrix"):
        SourmashDereplicator()._dense_dereplicate(genomes, tmp_path / "o", 31, 1000, 0.9, _LOG)
