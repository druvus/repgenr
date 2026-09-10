"""FastTree honours --bootstrap as its resampling count (live audit finding)."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.treebuilders import fasttree as mod
from repgenr.treebuilders.base import TreeParams

_LOG = logging.getLogger("test")


def _msa(tmp_path: Path) -> Path:
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nACGT\n>b\nACGA\n", encoding="utf-8")
    return msa


def _run(tmp_path: Path, monkeypatch, bootstrap: int) -> list[str]:
    calls: list[list[str]] = []

    def fake(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        Path(k["stdout_path"]).write_text("(a,b);\n", encoding="utf-8")

    monkeypatch.setattr(mod, "run_tool", fake)
    mod.FasttreeBuilder().build(
        _msa(tmp_path), tmp_path / "out", TreeParams(bootstrap=bootstrap), _LOG
    )
    (cmd,) = calls
    return cmd


def test_bootstrap_sets_boot(tmp_path: Path, monkeypatch) -> None:
    cmd = _run(tmp_path, monkeypatch, 500)
    assert cmd[cmd.index("-boot") + 1] == "500"


def test_no_bootstrap_keeps_the_default(tmp_path: Path, monkeypatch) -> None:
    assert "-boot" not in _run(tmp_path, monkeypatch, 0)
