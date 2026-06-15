"""SQLite genome manifest.

Replaces the ``str(dict)`` / ``pickle`` state blobs and the repeated
``os.listdir`` scans that do not scale to thousands of genomes. The manifest is
the single source of truth for which genomes are selected, where their files
live, their taxonomy, and their dereplication status.
"""

from __future__ import annotations

import os
import sqlite3
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

MANIFEST_FILENAME = "manifest.sqlite"

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
    """Thin SQLite wrapper for the genome inventory."""

    def __init__(self, path: str | os.PathLike[str]):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(self.path)
        self._conn.row_factory = sqlite3.Row
        self._conn.executescript(_SCHEMA)
        self._conn.commit()

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
            conn.execute(
                """
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
                """,
                {
                    "accession": record.accession,
                    "filename": record.filename,
                    "source": record.source,
                    "family": record.family,
                    "genus": record.genus,
                    "species": record.species,
                    "is_outgroup": int(record.is_outgroup),
                    "derep_status": record.derep_status,
                    "representative": record.representative,
                },
            )

    def upsert_many(self, records: list[GenomeRecord]) -> None:
        for record in records:
            self.upsert(record)

    def set_derep_status(
        self, accession: str, status: str, representative: str | None = None
    ) -> None:
        with self.transaction() as conn:
            conn.execute(
                "UPDATE genomes SET derep_status=?, representative=? WHERE accession=?",
                (status, representative, accession),
            )

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
