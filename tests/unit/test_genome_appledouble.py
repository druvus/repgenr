"""Batch extraction ignores AppleDouble "._x.fna" side files (live audit finding).

On macOS volumes without native extended attributes every extracted file gets
a "._" twin; the twin is not FASTA, so it used to be reported as a discarded
download for the same accession on every run.
"""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.stages import genome


def test_appledouble_twins_are_skipped(monkeypatch, tmp_path: Path, caplog) -> None:
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    dest = tmp_path / "genomes"
    dest.mkdir()

    def fake_run_cmd(argv, **kwargs):
        if argv[1] == "download":
            Path(argv[-1]).write_bytes(b"")
        else:  # rehydrate: lay out the extracted tree the way datasets does
            extract = Path(argv[-1])
            acc = extract / "ncbi_dataset" / "data" / "GCF_000001.1"
            acc.mkdir(parents=True)
            (acc / "GCF_000001.1_genomic.fna").write_text(">x\nACGT\n", encoding="utf-8")
            (acc / "._GCF_000001.1_genomic.fna").write_bytes(b"\x00\x05\x16\x07")

    monkeypatch.setattr(genome, "_run_cmd", fake_run_cmd)
    monkeypatch.setattr(genome.process, "unzip", lambda zip_path, extract: extract.mkdir())

    with caplog.at_level(logging.WARNING):
        missing = genome._download_one_batch(
            ["GCF_000001.1"],
            {"GCF_000001.1": "Fam_Gen_sp_GCF_000001.1.fasta"},
            dest,
            scratch,
            logging.getLogger("test"),
            False,
            0,
        )
    assert missing == []
    assert (dest / "Fam_Gen_sp_GCF_000001.1.fasta").read_text(encoding="utf-8") == ">x\nACGT\n"
    assert "Discarding non-FASTA" not in caplog.text
