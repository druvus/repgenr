"""dRep genome staging: hardlink plain FASTAs, decompress only .gz."""

from __future__ import annotations

import gzip
import os
from pathlib import Path

from repgenr.dereplicators.drep import _stage_genome


def test_stage_plain_fasta_is_hardlinked(tmp_path: Path) -> None:
    src = tmp_path / "g.fasta"
    src.write_text(">x\nACGT\n")
    dest = tmp_path / "staged"
    dest.mkdir()
    out = _stage_genome(src, dest)
    assert out.read_text() == ">x\nACGT\n"
    # hardlink shares the inode (no extra disk), not a separate copy
    assert os.path.samefile(out, src) or out.stat().st_ino == src.stat().st_ino


def test_stage_gz_is_decompressed(tmp_path: Path) -> None:
    src = tmp_path / "g.fasta.gz"
    with gzip.open(src, "wb") as fo:
        fo.write(b">x\nACGT\n")
    dest = tmp_path / "staged"
    dest.mkdir()
    out = _stage_genome(src, dest)
    assert out.name == "g.fasta"  # .gz stripped
    assert out.read_text() == ">x\nACGT\n"
