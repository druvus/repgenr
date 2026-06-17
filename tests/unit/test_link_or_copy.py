"""Unit tests for link_or_copy (hardlink staging with copy fallback)."""

from __future__ import annotations

import os
from pathlib import Path

from repgenr.core.process import link_or_copy


def test_hardlinks_same_filesystem(tmp_path: Path) -> None:
    src = tmp_path / "g.fasta"
    src.write_text(">s\nACGT\n")
    dst = tmp_path / "out" / "g.fasta"
    dst.parent.mkdir()
    link_or_copy(src, dst)
    assert dst.read_text() == ">s\nACGT\n"
    # a hardlink shares the inode with the source
    assert os.stat(src).st_ino == os.stat(dst).st_ino


def test_falls_back_to_copy_when_link_fails(tmp_path: Path, monkeypatch) -> None:
    src = tmp_path / "g.fasta"
    src.write_text(">s\nACGT\n")
    dst = tmp_path / "g_copy.fasta"

    def boom(*_a, **_k):
        raise OSError("cross-device link not permitted")

    monkeypatch.setattr(os, "link", boom)
    link_or_copy(src, dst)
    assert dst.read_text() == ">s\nACGT\n"
    # a real copy -> distinct inode
    assert os.stat(src).st_ino != os.stat(dst).st_ino
