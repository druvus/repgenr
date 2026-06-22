"""sourmash sketch caching: reuse signatures across --target-reps iterations."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.dereplicators import sourmash
from repgenr.dereplicators.sourmash import SourmashDereplicator

_LOG = logging.getLogger("test")


def _genomes(tmp_path: Path) -> list[Path]:
    out = []
    for name in ("g1.fasta", "g2.fasta"):
        p = tmp_path / name
        p.write_text(">x\nACGT\n")
        out.append(p)
    return out


def _fake_run_tool(calls: list[str], genomes: list[Path]):
    def run_tool(caps, cmd, **kwargs):
        parts = [str(c) for c in cmd]
        if "sketch" in parts:
            calls.append("sketch")
            outdir = Path(parts[parts.index("--outdir") + 1])
            outdir.mkdir(parents=True, exist_ok=True)
            for g in genomes:
                (outdir / f"{g.stem}.sig").write_text("sig")
        elif "compare" in parts:
            calls.append("compare")
            csv_path = Path(parts[parts.index("--csv") + 1])
            csv_path.write_text("g1,g2\n1.0,1.0\n1.0,1.0\n")
        return 0

    return run_tool


def test_dense_sketches_when_cache_empty(tmp_path: Path, monkeypatch) -> None:
    genomes = _genomes(tmp_path)
    calls: list[str] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls, genomes))
    cache = tmp_path / "sketches"

    clusters, _ = SourmashDereplicator()._dense_dereplicate(
        genomes, tmp_path / "out0", 31, 1000, 0.9, _LOG, sketch_cache=cache
    )
    assert "sketch" in calls  # nothing cached yet -> sketched
    assert len(clusters) == 1  # identical genomes collapse to one representative


def test_dense_reuses_cached_signatures(tmp_path: Path, monkeypatch) -> None:
    genomes = _genomes(tmp_path)
    cache = tmp_path / "sketches"
    cache.mkdir()
    for g in genomes:  # pre-populate as a prior iteration would have
        (cache / f"{g.stem}.sig").write_text("sig")

    calls: list[str] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls, genomes))

    SourmashDereplicator()._dense_dereplicate(
        genomes, tmp_path / "out1", 31, 1000, 0.9, _LOG, sketch_cache=cache
    )
    assert "sketch" not in calls  # reused the cache, only compared
    assert "compare" in calls
