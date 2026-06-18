"""Unit tests for the SQLite manifest, incl. batched writes."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.errors import WorkdirError
from repgenr.core.manifest import (
    BUSY_TIMEOUT_MS,
    SCHEMA_VERSION,
    GenomeRecord,
    Manifest,
)


def _rec(i: int) -> GenomeRecord:
    return GenomeRecord(accession=f"GCF_{i:06d}.1", filename=f"f_{i}.fasta", source="gtdb")


def test_schema_version_and_busy_timeout(tmp_path: Path) -> None:
    m = Manifest(tmp_path / "manifest.sqlite")
    assert int(m._conn.execute("PRAGMA user_version").fetchone()[0]) == SCHEMA_VERSION
    assert int(m._conn.execute("PRAGMA busy_timeout").fetchone()[0]) == BUSY_TIMEOUT_MS
    m.close()


def test_adopts_pre_versioning_db(tmp_path: Path) -> None:
    # A pre-versioning DB reports user_version=0; opening it adopts v1 and keeps data.
    p = tmp_path / "manifest.sqlite"
    m = Manifest(p)
    m.upsert_many([_rec(1)])
    m._conn.execute("PRAGMA user_version=0")  # simulate an old, unversioned DB
    m._conn.commit()
    m.close()
    m2 = Manifest(p)
    assert int(m2._conn.execute("PRAGMA user_version").fetchone()[0]) == SCHEMA_VERSION
    assert m2.count() == 1
    m2.close()


def test_rejects_newer_schema(tmp_path: Path) -> None:
    p = tmp_path / "manifest.sqlite"
    m = Manifest(p)
    m._conn.execute(f"PRAGMA user_version={SCHEMA_VERSION + 1}")
    m._conn.commit()
    m.close()
    with pytest.raises(WorkdirError, match="newer than this"):
        Manifest(p)


def test_wal_mode_enabled(tmp_path: Path) -> None:
    m = Manifest(tmp_path / "manifest.sqlite")
    mode = m._conn.execute("PRAGMA journal_mode").fetchone()[0]
    assert mode.lower() == "wal"
    m.close()


def test_upsert_many_batches_and_persists(tmp_path: Path) -> None:
    m = Manifest(tmp_path / "manifest.sqlite")
    m.upsert_many([_rec(i) for i in range(50)])
    assert m.count() == 50
    # upsert is an upsert: re-inserting updates, does not duplicate
    again = _rec(0)
    again.species = "updated"
    m.upsert_many([again])
    assert m.count() == 50
    rows = {g.accession: g for g in m.all_genomes()}
    assert rows["GCF_000000.1"].species == "updated"
    m.close()


def test_set_derep_status_many(tmp_path: Path) -> None:
    m = Manifest(tmp_path / "manifest.sqlite")
    m.upsert_many([_rec(i) for i in range(5)])
    m.set_derep_status_many([
        ("GCF_000000.1", "representative", None),
        ("GCF_000001.1", "contained", "GCF_000000.1"),
        ("GCF_999999.1", "contained", "GCF_000000.1"),  # absent -> no-op, not an error
    ])
    reps = {g.accession for g in m.representatives()}
    assert reps == {"GCF_000000.1"}
    by_acc = {g.accession: g for g in m.all_genomes()}
    assert by_acc["GCF_000001.1"].derep_status == "contained"
    assert by_acc["GCF_000001.1"].representative == "GCF_000000.1"
    m.close()
