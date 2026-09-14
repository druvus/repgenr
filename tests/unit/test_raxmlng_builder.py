"""RAxML-NG builder: publish the support-annotated tree when --all produced one."""

from __future__ import annotations

import logging
from pathlib import Path

import repgenr.treebuilders.raxmlng as mod
from repgenr.treebuilders.base import TreeParams

_BEST = "((a,b),(c,d));"
_SUPPORT = "((a,b)97,(c,d)100);"


def _fake_run(outputs: dict[str, str]):
    def run_tool(caps, command, *, logger, **kwargs):
        cmd = [str(c) for c in command]
        prefix = cmd[cmd.index("--prefix") + 1]
        for suffix, text in outputs.items():
            Path(prefix + suffix).write_text(text, encoding="utf-8")
        return 0

    return run_tool


def test_publishes_support_tree_when_present(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        mod, "run_tool", _fake_run({".raxml.bestTree": _BEST, ".raxml.support": _SUPPORT})
    )
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nAC\n", encoding="utf-8")
    tree = mod.RaxmlNgBuilder().build(msa, tmp_path / "out", TreeParams(), logging.getLogger())
    assert tree.read_text(encoding="utf-8") == _SUPPORT


def test_falls_back_to_best_tree_without_support(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(mod, "run_tool", _fake_run({".raxml.bestTree": _BEST}))
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nAC\n", encoding="utf-8")
    tree = mod.RaxmlNgBuilder().build(msa, tmp_path / "out", TreeParams(), logging.getLogger())
    assert tree.read_text(encoding="utf-8") == _BEST
