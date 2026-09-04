"""A genuine branchwater failure must surface, not silently switch algorithms.

Availability is pre-probed by ``_branchwater_available``, so an error raised
inside the sparse path is a real tool failure (bad parameters, OOM, broken
image). Silently degrading to the dense back-end can change which
representatives are picked; the fallback is therefore opt-in.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import ToolExecutionError
from repgenr.dereplicators import sourmash
from repgenr.dereplicators.base import DerepParams
from repgenr.dereplicators.sourmash import SourmashDereplicator

_LOG = logging.getLogger("test")


def _genomes(tmp_path: Path) -> list[Path]:
    out = []
    for name in ("g1.fasta", "g2.fasta"):
        p = tmp_path / name
        p.write_text(">x\nACGT\n")
        out.append(p)
    return out


def _failing_sparse_fake(calls: list[str]):
    """manysketch fails; dense sketch/compare succeed if reached."""

    def run_tool(caps, cmd, **kwargs):
        parts = [str(c) for c in cmd]
        if "manysketch" in parts:
            calls.append("manysketch")
            raise ToolExecutionError(parts, 137, output="killed")
        if "sketch" in parts:
            calls.append("sketch")
            fofn = Path(parts[parts.index("--from-file") + 1])
            outdir = Path(parts[parts.index("--outdir") + 1])
            outdir.mkdir(parents=True, exist_ok=True)
            for p in fofn.read_text(encoding="utf-8").splitlines():
                (outdir / f"{Path(p).stem}.sig").write_text("sig")
        elif "compare" in parts:
            calls.append("compare")
            fofn = Path(parts[parts.index("--from-file") + 1])
            labels = [Path(s).stem for s in fofn.read_text().splitlines()]
            n = len(labels)
            body = "\n".join(",".join("1.0" for _ in range(n)) for _ in range(n))
            Path(parts[parts.index("--csv") + 1]).write_text(
                ",".join(labels) + "\n" + body + "\n", encoding="utf-8"
            )
        return 0

    return run_tool


def _dereplicate(tmp_path: Path, extra: dict | None = None):
    params = DerepParams(primary_ani=0.99, secondary_ani=0.99, threads=1, extra=dict(extra or {}))
    return SourmashDereplicator().dereplicate(_genomes(tmp_path), tmp_path / "out", params, _LOG)


def test_sparse_failure_propagates(tmp_path: Path, monkeypatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(sourmash, "_branchwater_available", lambda caps, log: True)
    monkeypatch.setattr(sourmash, "run_tool", _failing_sparse_fake(calls))

    with pytest.raises(ToolExecutionError):
        _dereplicate(tmp_path)
    assert "sketch" not in calls  # dense path was NOT silently attempted


def test_dense_fallback_via_extra(tmp_path: Path, monkeypatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(sourmash, "_branchwater_available", lambda caps, log: True)
    monkeypatch.setattr(sourmash, "run_tool", _failing_sparse_fake(calls))

    result = _dereplicate(tmp_path, extra={"dense_fallback": True})
    assert calls[0] == "manysketch" and "compare" in calls
    assert len(result.representatives) == 1


def test_dense_fallback_via_env(tmp_path: Path, monkeypatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(sourmash, "_branchwater_available", lambda caps, log: True)
    monkeypatch.setattr(sourmash, "run_tool", _failing_sparse_fake(calls))
    monkeypatch.setenv("REPGENR_SOURMASH_DENSE_FALLBACK", "1")

    result = _dereplicate(tmp_path)
    assert "compare" in calls
    assert len(result.representatives) == 1
