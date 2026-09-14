"""Contig filtering and assembly summary shared by every assembler."""

from __future__ import annotations

from pathlib import Path

from repgenr.assemblers.contigs import filter_contigs


def test_filter_drops_short_contigs_renames_headers_and_summarises(tmp_path: Path) -> None:
    src = tmp_path / "raw.fa"
    src.write_text(
        ">NODE_1 cov=3\nACGT" * 1 + "A" * 996 + "\n>NODE_2\nACGTAC\n>NODE_3 x\n" + "G" * 600 + "\n",
        encoding="utf-8",
    )
    dst = tmp_path / "SRR1.fasta"
    stats = filter_contigs(src, dst, min_length=500, prefix="SRR1")
    text = dst.read_text(encoding="utf-8")
    assert text.startswith(">SRR1_contig1\n") and ">SRR1_contig2\n" in text and "NODE_2" not in text
    assert (stats.n_contigs, stats.total_length, stats.largest_contig) == (2, 1600, 1000)
    assert stats.n50 == 1000  # the 1000 bp contig alone covers half the total


def test_n50_on_equal_contigs(tmp_path: Path) -> None:
    src = tmp_path / "raw.fa"
    src.write_text(">a\n" + "A" * 700 + "\n>b\n" + "C" * 700 + "\n", encoding="utf-8")
    stats = filter_contigs(src, tmp_path / "out.fa", min_length=500, prefix="x")
    assert (stats.n_contigs, stats.n50) == (2, 700)
