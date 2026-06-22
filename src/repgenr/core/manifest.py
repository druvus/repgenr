"""SQLite genome manifest.

Replaces the ``str(dict)`` / ``pickle`` state blobs and the repeated
``os.listdir`` scans that do not scale to thousands of genomes. The manifest is
the single source of truth for which genomes are selected, where their files
live, their taxonomy, and their dereplication status.
"""

from __future__ import annotations

import os
import sqlite3
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

from .errors import WorkdirError

MANIFEST_FILENAME = "manifest.sqlite"
SCHEMA_VERSION = 1  # bump + add a migration step when the table layout changes
BUSY_TIMEOUT_MS = 30000  # wait up to 30s for a competing writer before erroring

_UPSERT_SQL = """
    INSERT INTO genomes (accession, filename, source, family, genus,
                         species, is_outgroup, derep_status, representative)
    VALUES (:accession, :filename, :source, :family, :genus,
            :species, :is_outgroup, :derep_status, :representative)
    ON CONFLICT(accession) DO UPDATE SET
        filename=excluded.filename,
        source=excluded.source,
        family=excluded.family,
        genus=excluded.genus,
        species=excluded.species,
        is_outgroup=excluded.is_outgroup,
        derep_status=excluded.derep_status,
        representative=excluded.representative
"""

_SET_DEREP_SQL = "UPDATE genomes SET derep_status=?, representative=? WHERE accession=?"

_SCHEMA = """
CREATE TABLE IF NOT EXISTS genomes (
    accession   TEXT PRIMARY KEY,
    filename    TEXT,
    source      TEXT,                  -- gtdb | bvbrc | ncbi
    family      TEXT,
    genus       TEXT,
    species     TEXT,
    is_outgroup INTEGER DEFAULT 0,
    derep_status TEXT,                 -- representative | contained | fail_qc | NULL
    representative TEXT                 -- accession of the representative, if contained
);
CREATE INDEX IF NOT EXISTS idx_genomes_species ON genomes(species);
CREATE INDEX IF NOT EXISTS idx_genomes_derep ON genomes(derep_status);
"""


@dataclass
class GenomeRecord:
    accession: str
    filename: str | None = None
    source: str | None = None
    family: str | None = None
    genus: str | None = None
    species: str | None = None
    is_outgroup: bool = False
    derep_status: str | None = None
    representative: str | None = None


class Manifest:
    """Thin SQLite wrapper for the genome inventory.

    Not thread-safe: the single connection (``check_same_thread`` default) must be
    used from the thread that opened it. Stages write the manifest on the main
    thread; do not call it from inside a ``parallel_map`` worker. WAL + the busy
    timeout cover concurrency across *processes* (Nextflow scatter), not threads.
    """

    def __init__(self, path: str | os.PathLike[str]):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(self.path)
        self._conn.row_factory = sqlite3.Row
        # WAL + synchronous=NORMAL: commits no longer pay a full fsync each, which
        # is the dominant cost for many small writes. The manifest is a workdir
        # artifact (regenerable from the stages), so the NORMAL durability
        # trade-off -- a power loss can lose only the last transaction -- is fine.
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA synchronous=NORMAL")
        # Concurrent writers (parallel stages, two invocations on one workdir)
        # wait for the lock instead of failing immediately with "database is
        # locked".
        self._conn.execute(f"PRAGMA busy_timeout={BUSY_TIMEOUT_MS}")
        self._conn.executescript(_SCHEMA)
        self._migrate()
        self._conn.commit()

    def _migrate(self) -> None:
        """Apply schema migrations keyed on ``PRAGMA user_version``.

        Existing pre-versioning databases report user_version=0; their layout
        already matches v1, so they are adopted as v1. Future schema changes add
        a numbered migration step and bump SCHEMA_VERSION.
        """
        version = int(self._conn.execute("PRAGMA user_version").fetchone()[0])
        if version > SCHEMA_VERSION:
            raise WorkdirError(
                f"Manifest at {self.path} has schema version {version}, newer than this "
                f"RepGenR supports ({SCHEMA_VERSION}). Upgrade RepGenR or use a new workdir."
            )
        # (no v0->v1 data change: the CREATE IF NOT EXISTS schema is v1)
        # Future: while version < SCHEMA_VERSION: apply step; version += 1
        if version != SCHEMA_VERSION:
            self._conn.execute(f"PRAGMA user_version={SCHEMA_VERSION}")

    @classmethod
    def open(cls, workdir: str | os.PathLike[str]) -> Manifest:
        return cls(Path(workdir) / MANIFEST_FILENAME)

    def close(self) -> None:
        self._conn.close()

    def __enter__(self) -> Manifest:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    @contextmanager
    def transaction(self) -> Iterator[sqlite3.Connection]:
        try:
            yield self._conn
            self._conn.commit()
        except Exception:
            self._conn.rollback()
            raise

    def upsert(self, record: GenomeRecord) -> None:
        with self.transaction() as conn:
            conn.execute(_UPSERT_SQL, _record_params(record))

    def upsert_many(self, records: list[GenomeRecord]) -> None:
        # One transaction (one commit/fsync) for the whole batch -- committing per
        # record is ~0.5 ms each, i.e. seconds of fsync overhead at 1000s genomes.
        with self.transaction() as conn:
            conn.executemany(_UPSERT_SQL, [_record_params(r) for r in records])

    def set_derep_status(
        self, accession: str, status: str, representative: str | None = None
    ) -> None:
        with self.transaction() as conn:
            conn.execute(_SET_DEREP_SQL, (status, representative, accession))

    def set_derep_status_many(
        self, updates: Sequence[tuple[str, str, str | None]]
    ) -> None:
        """Batch derep-status updates (accession, status, representative) in one
        transaction; avoids one commit per genome on large sets."""
        rows = [(status, rep, accession) for accession, status, rep in updates]
        with self.transaction() as conn:
            conn.executemany(_SET_DEREP_SQL, rows)

    def count(self) -> int:
        cur = self._conn.execute("SELECT COUNT(*) AS n FROM genomes WHERE is_outgroup=0")
        return int(cur.fetchone()["n"])

    def representatives(self) -> list[GenomeRecord]:
        cur = self._conn.execute(
            "SELECT * FROM genomes WHERE derep_status='representative'"
        )
        return [_row_to_record(row) for row in cur.fetchall()]

    def all_genomes(self, include_outgroup: bool = False) -> list[GenomeRecord]:
        query = "SELECT * FROM genomes"
        if not include_outgroup:
            query += " WHERE is_outgroup=0"
        cur = self._conn.execute(query)
        return [_row_to_record(row) for row in cur.fetchall()]


def _record_params(record: GenomeRecord) -> dict:
    return {
        "accession": record.accession,
        "filename": record.filename,
        "source": record.source,
        "family": record.family,
        "genus": record.genus,
        "species": record.species,
        "is_outgroup": int(record.is_outgroup),
        "derep_status": record.derep_status,
        "representative": record.representative,
    }


def _row_to_record(row: sqlite3.Row) -> GenomeRecord:
    return GenomeRecord(
        accession=row["accession"],
        filename=row["filename"],
        source=row["source"],
        family=row["family"],
        genus=row["genus"],
        species=row["species"],
        is_outgroup=bool(row["is_outgroup"]),
        derep_status=row["derep_status"],
        representative=row["representative"],
    )
