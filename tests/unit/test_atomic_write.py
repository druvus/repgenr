"""Atomic deliverable writes: a failed write never destroys the previous file."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.contracts import atomic_path, atomic_replace, write_clusters


def test_atomic_replace_success(tmp_path: Path) -> None:
    target = tmp_path / "out.tsv"
    with atomic_replace(target) as fo:
        fo.write("new content\n")
    assert target.read_text(encoding="utf-8") == "new content\n"
    assert list(tmp_path.iterdir()) == [target]  # no temp leftovers


def test_atomic_replace_failure_preserves_previous(tmp_path: Path) -> None:
    target = tmp_path / "out.tsv"
    target.write_text("previous good content\n", encoding="utf-8")
    with pytest.raises(RuntimeError):
        with atomic_replace(target) as fo:
            fo.write("half-writ")
            raise RuntimeError("crash mid-write")
    assert target.read_text(encoding="utf-8") == "previous good content\n"
    assert list(tmp_path.iterdir()) == [target]  # temp cleaned up


def test_atomic_path_success_and_failure(tmp_path: Path) -> None:
    target = tmp_path / "tree.nwk"
    target.write_text("(old);\n", encoding="utf-8")

    with atomic_path(target) as tmp:
        tmp.write_text("(new);\n", encoding="utf-8")
    assert target.read_text(encoding="utf-8") == "(new);\n"

    with pytest.raises(RuntimeError):
        with atomic_path(target) as tmp:
            tmp.write_text("(broken", encoding="utf-8")
            raise RuntimeError("tool failed")
    assert target.read_text(encoding="utf-8") == "(new);\n"
    assert list(tmp_path.iterdir()) == [target]


def test_atomic_path_failure_without_tmp_written(tmp_path: Path) -> None:
    target = tmp_path / "tree.nwk"
    with pytest.raises(RuntimeError):
        with atomic_path(target):
            raise RuntimeError("failed before writing anything")
    assert not target.exists()
    assert list(tmp_path.iterdir()) == []


def test_contract_writer_failure_preserves_previous_pair(tmp_path: Path, monkeypatch) -> None:
    """A crash inside a contract writer leaves the previous file intact."""
    clusters = tmp_path / "clusters.tsv"
    write_clusters(clusters, {"rep.fasta": ["m1.fasta"]})
    before = clusters.read_text(encoding="utf-8")

    class Boom:
        def writerow(self, row):
            raise OSError("disk full")

    import repgenr.core.contracts as contracts

    monkeypatch.setattr(contracts.csv, "writer", lambda fo, **kw: Boom())
    with pytest.raises(OSError):
        write_clusters(clusters, {"other.fasta": []})
    assert clusters.read_text(encoding="utf-8") == before


def test_remove_tree_tolerates_entries_that_vanish(tmp_path, monkeypatch) -> None:
    """A file removed by the OS mid-walk (an AppleDouble twin on an external
    volume) is not an error; other failures still raise (live audit)."""
    from repgenr.core.process import _ignore_vanished, remove_tree

    d = tmp_path / "scratch"
    (d / "sub").mkdir(parents=True)
    (d / "sub" / "a.txt").write_text("x", encoding="utf-8")
    remove_tree(d)
    assert not d.exists()
    remove_tree(d)  # already gone: a no-op
    _ignore_vanished(None, "gone", FileNotFoundError("gone"))
    import pytest

    with pytest.raises(PermissionError):
        _ignore_vanished(None, "kept", PermissionError("kept"))
