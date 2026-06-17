"""Micro-benchmark for manifest (SQLite) write paths at scale.

Measures the per-commit cost of bulk upserts and derep-status updates on a real
on-disk SQLite file (so fsync cost is realistic). Run before/after batching.

    python scripts/bench_manifest.py
"""

from __future__ import annotations

import tempfile
import time
from pathlib import Path

from repgenr.core.manifest import GenomeRecord, Manifest


def _records(n: int) -> list[GenomeRecord]:
    return [
        GenomeRecord(
            accession=f"GCF_{i:09d}.1",
            filename=f"Fam_gen_sp_GCF_{i:09d}.1.fasta",
            source="gtdb",
            family="Fam",
            genus="gen",
            species="sp",
        )
        for i in range(n)
    ]


def bench_upsert_many(n: int) -> float:
    records = _records(n)
    with tempfile.TemporaryDirectory() as d:
        m = Manifest(Path(d) / "manifest.sqlite")
        start = time.perf_counter()
        m.upsert_many(records)
        elapsed = time.perf_counter() - start
        m.close()
    print(f"  upsert_many        n={n:>5}: {elapsed:7.3f}s")
    return elapsed


def bench_set_derep_status(n: int) -> float:
    records = _records(n)
    with tempfile.TemporaryDirectory() as d:
        m = Manifest(Path(d) / "manifest.sqlite")
        m.upsert_many(records)
        start = time.perf_counter()
        for i in range(n):
            m.set_derep_status(f"GCF_{i:09d}.1", "contained", "GCF_000000000.1")
        elapsed = time.perf_counter() - start
        m.close()
    print(f"  set_derep_status   n={n:>5}: {elapsed:7.3f}s")
    return elapsed


def bench_set_derep_status_many(n: int) -> float:
    records = _records(n)
    updates = [(f"GCF_{i:09d}.1", "contained", "GCF_000000000.1") for i in range(n)]
    with tempfile.TemporaryDirectory() as d:
        m = Manifest(Path(d) / "manifest.sqlite")
        m.upsert_many(records)
        start = time.perf_counter()
        m.set_derep_status_many(updates)
        elapsed = time.perf_counter() - start
        m.close()
    print(f"  set_derep_status_many n={n:>5}: {elapsed:7.3f}s")
    return elapsed


if __name__ == "__main__":
    print("manifest writes (on-disk SQLite, WAL + synchronous=NORMAL):")
    for n in (1000, 2000, 5000):
        bench_upsert_many(n)
    for n in (1000, 2000, 5000):
        bench_set_derep_status_many(n)
    print("single-call path (still one commit per call, but cheap under WAL):")
    for n in (1000, 5000):
        bench_set_derep_status(n)
