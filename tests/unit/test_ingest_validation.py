"""Download/ingest validation: safe unzip, prune, FASTA sanity, disk preflight."""

from __future__ import annotations

import logging
import zipfile
from collections import namedtuple
from pathlib import Path

import pytest

from repgenr.core import process
from repgenr.core.errors import WorkdirError
from repgenr.stages import genome

_LOG = logging.getLogger("test")
_Usage = namedtuple("_Usage", "total used free")


def test_unzip_ok(tmp_path: Path) -> None:
    z = tmp_path / "ok.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr("a.txt", "hello")
    process.unzip(z, tmp_path / "out")
    assert (tmp_path / "out" / "a.txt").read_text() == "hello"


def test_unzip_corrupt_raises(tmp_path: Path) -> None:
    bad = tmp_path / "bad.zip"
    bad.write_bytes(b"not a zip file")
    with pytest.raises(WorkdirError, match="Corrupt or truncated"):
        process.unzip(bad, tmp_path / "out")


def test_prune_only_removes_stale_fasta(tmp_path: Path) -> None:
    (tmp_path / "keep.fasta").write_text(">k\nAC\n")
    (tmp_path / "stale.fasta").write_text(">s\nAC\n")
    (tmp_path / "notes.txt").write_text("user file")
    (tmp_path / "subdir").mkdir()
    genome._prune(tmp_path, {"keep.fasta"}, _LOG)
    remaining = {p.name for p in tmp_path.iterdir()}
    assert remaining == {"keep.fasta", "notes.txt", "subdir"}  # stale fasta gone, rest kept


def test_assert_fasta(tmp_path: Path) -> None:
    good = tmp_path / "g.fasta"
    good.write_text(">x\nACGT\n")
    genome._assert_fasta(good)  # no raise
    bad = tmp_path / "b.fasta"
    bad.write_text("<html>error</html>")
    with pytest.raises(WorkdirError, match="not FASTA"):
        genome._assert_fasta(bad)


def test_check_disk_floor_raises(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(genome.shutil, "disk_usage", lambda p: _Usage(100, 99, 10_000_000))
    with pytest.raises(WorkdirError, match="refusing to download"):
        genome._check_disk(tmp_path, 100, _LOG)


def test_check_disk_tight_warns(tmp_path: Path, monkeypatch, caplog) -> None:
    # above the 1 GB floor but below the per-genome estimate -> warn, no raise
    monkeypatch.setattr(genome.shutil, "disk_usage", lambda p: _Usage(0, 0, 2_000_000_000))
    with caplog.at_level(logging.WARNING, logger="test"):
        genome._check_disk(tmp_path, 1000, _LOG)
    assert any("Low disk" in r.message for r in caplog.records)


def test_check_disk_ample_ok(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(genome.shutil, "disk_usage", lambda p: _Usage(0, 0, 500_000_000_000))
    genome._check_disk(tmp_path, 1000, _LOG)  # no raise, no warning
