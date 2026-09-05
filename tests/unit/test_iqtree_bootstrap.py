"""IQ-TREE's ultrafast bootstrap needs at least 1000 replicates; refuse fewer up front."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError
from repgenr.treebuilders import iqtree as mod
from repgenr.treebuilders.base import TreeParams

_LOG = logging.getLogger("test")


def _msa(tmp_path: Path) -> Path:
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nACGT\n>b\nACGA\n", encoding="utf-8")
    return msa


def test_bootstrap_below_1000_is_refused_before_running(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", lambda caps, cmd, **k: calls.append([str(c) for c in cmd]))
    with pytest.raises(UserInputError, match="1000"):
        mod.IqtreeBuilder().build(_msa(tmp_path), tmp_path / "out", TreeParams(bootstrap=100), _LOG)
    assert calls == []


def test_bootstrap_of_1000_is_passed_through(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []

    def fake(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        Path(cmd[cmd.index("-s") + 1] + ".treefile").write_text("(a,b);\n", encoding="utf-8")

    monkeypatch.setattr(mod, "run_tool", fake)
    mod.IqtreeBuilder().build(_msa(tmp_path), tmp_path / "out", TreeParams(bootstrap=1000), _LOG)
    (cmd,) = calls
    assert cmd[cmd.index("-B") + 1] == "1000"
