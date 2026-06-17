"""Unit tests for the SQLite manifest, incl. batched writes."""

from __future__ import annotations

from pathlib import Path

from repgenr.core.manifest import GenomeRecord, Manifest


def _rec(i: int) -> GenomeRecord:
    return GenomeRecord(accession=f"GCF_{i:06d}.1", filename=f"f_{i}.fasta", source="gtdb")


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
