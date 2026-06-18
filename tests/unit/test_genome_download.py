"""Unit test for the genome-download batching."""

from __future__ import annotations

import logging

from repgenr.stages import genome

_LOG = logging.getLogger("test")


def test_download_splits_into_fixed_size_batches(monkeypatch) -> None:
    seen: list[list[str]] = []
    monkeypatch.setattr(genome, "_DOWNLOAD_BATCH_SIZE", 5)
    monkeypatch.setattr(
        genome, "_download_one_batch",
        lambda ctx, batch, filenames, logger, keep_files, bi: seen.append(batch),
    )
    accs = [f"GCF_{i:03d}.1" for i in range(12)]
    genome._download_batch(
        ctx=object(), accessions=accs, filenames={}, logger=_LOG, keep_files=False
    )
    # 12 accessions / batch 5 -> [5, 5, 2], covering every accession exactly once
    assert [len(b) for b in seen] == [5, 5, 2]
    assert [a for b in seen for a in b] == accs
