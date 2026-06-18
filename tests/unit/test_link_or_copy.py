"""Unit tests for link_or_copy (hardlink staging with copy fallback)."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

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


def test_stages_real_content_through_a_symlink(tmp_path: Path) -> None:
    # Regression: skDER emits representatives as symlinks; hardlinking the symlink
    # itself (os.link on macOS) yielded a 0-byte staged genome. link_or_copy must
    # resolve to the real target and stage its content.
    real = tmp_path / "genome.fasta"
    real.write_text(">x\nACGTACGT\n")
    link = tmp_path / "rep_link.fasta"
    os.symlink(real, link)

    dst = tmp_path / "staged.fasta"
    link_or_copy(link, dst)

    assert dst.is_file()
    assert dst.stat().st_size == real.stat().st_size
    assert dst.read_text() == ">x\nACGTACGT\n"


def test_step_contract_rejects_empty_representative(tmp_path: Path) -> None:
    from repgenr.core.errors import WorkdirError
    from repgenr.dereplicators.base import STATUS_REPRESENTATIVE, DerepResult
    from repgenr.stages.derep_steps import _write_step_contract

    empty = tmp_path / "empty.fasta"
    empty.touch()  # zero length
    result = DerepResult(
        representatives=[empty],
        clusters={empty.name: []},
        genome_status={empty.name: STATUS_REPRESENTATIVE},
    )
    with pytest.raises(WorkdirError, match="empty"):
        _write_step_contract(tmp_path / "out", result, [])
