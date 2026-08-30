"""Digest helpers for input-aware resume fingerprints (core.inputs)."""

from __future__ import annotations

import os
from pathlib import Path

from repgenr.core.inputs import (
    ABSENT,
    dir_stat_digest,
    file_digest,
    inputs_digest,
    manifest_digest,
    path_digest,
)
from repgenr.core.manifest import GenomeRecord, Manifest


def _write(path: Path, content: str = "x") -> None:
    path.write_text(content, encoding="utf-8")


# --- file_digest --------------------------------------------------------------


def test_file_digest_changes_with_content(tmp_path) -> None:
    f = tmp_path / "a.tsv"
    _write(f, "one")
    d1 = file_digest(f)
    _write(f, "two")
    assert file_digest(f) != d1


def test_file_digest_stable_for_same_content(tmp_path) -> None:
    f = tmp_path / "a.tsv"
    _write(f, "same")
    assert file_digest(f) == file_digest(f)


def test_file_digest_absent_sentinel(tmp_path) -> None:
    assert file_digest(tmp_path / "missing.tsv") == ABSENT


# --- dir_stat_digest ----------------------------------------------------------


def test_dir_digest_absent_sentinel(tmp_path) -> None:
    assert dir_stat_digest(tmp_path / "nope") == ABSENT


def test_dir_digest_changes_on_add_remove_rename(tmp_path) -> None:
    d = tmp_path / "genomes"
    d.mkdir()
    _write(d / "g1.fasta")
    base = dir_stat_digest(d)

    _write(d / "g2.fasta")
    added = dir_stat_digest(d)
    assert added != base

    (d / "g2.fasta").rename(d / "g2renamed.fasta")
    renamed = dir_stat_digest(d)
    assert renamed != added

    (d / "g2renamed.fasta").unlink()
    assert dir_stat_digest(d) == base


def test_dir_digest_changes_on_resize_and_mtime(tmp_path) -> None:
    d = tmp_path / "genomes"
    d.mkdir()
    f = d / "g1.fasta"
    _write(f, "aa")
    base = dir_stat_digest(d)

    _write(f, "aaaa")  # size change
    resized = dir_stat_digest(d)
    assert resized != base

    # same size, bumped mtime (a stage rewriting a file always bumps mtime_ns)
    os.utime(f, ns=(1_000_000_000, 2_000_000_000))
    assert dir_stat_digest(d) != resized


def test_dir_digest_ignores_dotfiles(tmp_path) -> None:
    d = tmp_path / "genomes"
    d.mkdir()
    _write(d / "g1.fasta")
    base = dir_stat_digest(d)
    _write(d / "._g1.fasta")
    _write(d / ".DS_Store")
    assert dir_stat_digest(d) == base


# --- path_digest / inputs_digest ----------------------------------------------


def test_path_digest_dispatches(tmp_path) -> None:
    d = tmp_path / "adir"
    d.mkdir()
    f = tmp_path / "afile"
    _write(f)
    assert path_digest(d) == dir_stat_digest(d)
    assert path_digest(f) == file_digest(f)
    assert path_digest(tmp_path / "missing") == ABSENT


def test_inputs_digest_keys_are_workdir_relative(tmp_path) -> None:
    (tmp_path / "genomes").mkdir()
    _write(tmp_path / "selection.tsv")
    digests = inputs_digest(tmp_path, [tmp_path / "genomes", tmp_path / "selection.tsv"])
    assert set(digests) == {"genomes", "selection.tsv"}
    assert all(isinstance(v, str) and v for v in digests.values())


# --- manifest_digest ----------------------------------------------------------


def test_manifest_digest_order_independent_and_content_sensitive(tmp_path) -> None:
    m1 = Manifest(tmp_path / "a.sqlite")
    m1.upsert_many([
        GenomeRecord(accession="GCF_2", species="b"),
        GenomeRecord(accession="GCF_1", species="a"),
    ])
    m2 = Manifest(tmp_path / "b.sqlite")
    m2.upsert_many([
        GenomeRecord(accession="GCF_1", species="a"),
        GenomeRecord(accession="GCF_2", species="b"),
    ])
    assert manifest_digest(m1) == manifest_digest(m2)

    m2.upsert_many([GenomeRecord(accession="GCF_3", species="c")])
    assert manifest_digest(m1) != manifest_digest(m2)
