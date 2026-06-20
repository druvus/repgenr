"""Aligner parameter passthrough: divergence-tuning flags reach the tool argv."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.aligners import progressivemauve as pm
from repgenr.aligners import sibeliaz as sz
from repgenr.aligners.base import AlignParams

_LOG = logging.getLogger("test")


def _genomes(tmp_path: Path, n: int) -> list[Path]:
    out = []
    for i in range(n):
        p = tmp_path / f"g{i}.fasta"
        p.write_text(">x\nACGTACGT\n")
        out.append(p)
    return out


def test_sibeliaz_forwards_tuning_flags(monkeypatch, tmp_path: Path) -> None:
    captured: dict[str, list[str]] = {}

    def fake_run(caps, cmd, **kw):
        captured["cmd"] = [str(c) for c in cmd]
        raise RuntimeError("stop after capture")

    monkeypatch.setattr(sz, "run_tool", fake_run)
    monkeypatch.setattr(sz, "_sibeliaz_invocation", lambda out_dir, logger: ["sibeliaz"])

    aln = sz.SibeliazAligner()
    params = AlignParams(threads=2, extra={"kmer": 15, "bubble": 200})
    with pytest.raises(RuntimeError):
        aln.align(_genomes(tmp_path, 3), None, tmp_path / "out", params, _LOG)

    cmd = captured["cmd"]
    assert "-k" in cmd and cmd[cmd.index("-k") + 1] == "15"
    assert "-b" in cmd and cmd[cmd.index("-b") + 1] == "200"
    # default: no -k when unset
    captured.clear()
    with pytest.raises(RuntimeError):
        aln.align(_genomes(tmp_path, 3), None, tmp_path / "out2", AlignParams(), _LOG)
    assert "-k" not in captured["cmd"]


def test_progressivemauve_forwards_seed_weight(monkeypatch, tmp_path: Path) -> None:
    captured: dict[str, list[str]] = {}

    def fake_run(caps, cmd, **kw):
        captured["cmd"] = [str(c) for c in cmd]
        raise RuntimeError("stop after capture")

    monkeypatch.setattr(pm, "run_tool", fake_run)

    aln = pm.ProgressiveMauveAligner()
    # 3 genomes (mauve needs >=3); threads=1 -> sequential -> exception propagates
    params = AlignParams(threads=1, extra={"seed_weight": 11})
    with pytest.raises(RuntimeError):
        aln.align(_genomes(tmp_path, 3), None, tmp_path / "out", params, _LOG)

    cmd = captured["cmd"]
    assert "--seed-weight" in cmd and cmd[cmd.index("--seed-weight") + 1] == "11"
