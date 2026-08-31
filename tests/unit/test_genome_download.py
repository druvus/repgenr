"""Unit test for the genome-download batching."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.stages import genome

_LOG = logging.getLogger("test")


def test_download_splits_into_fixed_size_batches(monkeypatch, tmp_path: Path) -> None:
    seen: list[list[str]] = []
    monkeypatch.setattr(genome, "_DOWNLOAD_BATCH_SIZE", 5)
    # hermetic: the real free-disk guard must not fail this test on a full host
    monkeypatch.setattr(genome, "_check_disk", lambda scratch, n, logger: None)
    def fake_batch(batch, filenames, dest_dir, scratch_dir, logger, keep_files, bi):
        seen.append(batch)
        return []  # no missing accessions

    monkeypatch.setattr(genome, "_download_one_batch", fake_batch)
    accs = [f"GCF_{i:03d}.1" for i in range(12)]
    genome.download_accessions(
        accs, {}, tmp_path / "dest", tmp_path / "scratch", _LOG, keep_files=False
    )
    # 12 accessions / batch 5 -> [5, 5, 2], covering every accession exactly once
    assert [len(b) for b in seen] == [5, 5, 2]
    assert [a for b in seen for a in b] == accs
