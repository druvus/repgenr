"""Unit tests for the SQLite manifest, incl. batched writes."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from repgenr.core import manifest as manifest_module
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


def test_readonly_rejects_newer_schema_without_leaking_the_connection(
    tmp_path: Path, monkeypatch
) -> None:
    """A version rejection in the read-only branch must still close the handle."""
    p = tmp_path / "manifest.sqlite"
    m = Manifest(p)
    m._conn.execute(f"PRAGMA user_version={SCHEMA_VERSION + 1}")
    m._conn.commit()
    m.close()

    opened: list[sqlite3.Connection] = []
    real_connect = sqlite3.connect

    def tracking_connect(*args, **kwargs):  # noqa: ANN002, ANN003, ANN202
        conn = real_connect(*args, **kwargs)
        opened.append(conn)
        return conn

    monkeypatch.setattr(manifest_module.sqlite3, "connect", tracking_connect)
    with pytest.raises(WorkdirError, match="newer than this"):
        Manifest.open_readonly(p)

    assert opened, "the read-only branch never opened a connection"
    for conn in opened:
        with pytest.raises(sqlite3.ProgrammingError):
            conn.execute("SELECT 1")  # a closed connection refuses further use


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


def test_replace_genomes_deletes_deselected(tmp_path: Path) -> None:
    m = Manifest(tmp_path / "m.sqlite")
    m.upsert_many([
        GenomeRecord(accession="A1", species="a"),
        GenomeRecord(accession="A2", species="b"),
        GenomeRecord(accession="OG", species="og", is_outgroup=True),
    ])
    # re-selection keeps A1, drops A2, new outgroup
    m.replace_genomes([
        GenomeRecord(accession="A1", species="a"),
        GenomeRecord(accession="OG2", species="og2", is_outgroup=True),
    ])
    accs = {g.accession for g in m.all_genomes(include_outgroup=True)}
    assert accs == {"A1", "OG2"}


def test_set_derep_status_many_clears_stale_rows(tmp_path: Path) -> None:
    m = Manifest(tmp_path / "m.sqlite")
    m.upsert_many([
        GenomeRecord(accession="A1"),
        GenomeRecord(accession="A2"),
        GenomeRecord(accession="A3"),
    ])
    m.set_derep_status_many([
        ("A1", "representative", None),
        ("A2", "contained", "A1"),
        ("A3", "representative", None),
    ])
    # a narrower re-dereplication: A3 no longer appears in the result at all
    m.set_derep_status_many([
        ("A1", "representative", None),
        ("A2", "contained", "A1"),
    ])
    by_acc = {g.accession: g for g in m.all_genomes(include_outgroup=True)}
    assert by_acc["A3"].derep_status is None, "stale representative status cleared"
    assert by_acc["A3"].representative is None
    assert by_acc["A1"].derep_status == "representative"


def test_manifest_digest_changes_after_reconciliation(tmp_path: Path) -> None:
    from repgenr.core.inputs import manifest_digest

    m = Manifest(tmp_path / "m.sqlite")
    m.upsert_many([GenomeRecord(accession="A1"), GenomeRecord(accession="A2")])
    before = manifest_digest(m)
    m.replace_genomes([GenomeRecord(accession="A1")])
    assert manifest_digest(m) != before
