"""Contig filtering and assembly summary shared by every assembler."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..core.contracts import atomic_replace


@dataclass(frozen=True)
class ContigStats:
    n_contigs: int
    total_length: int
    n50: int
    largest_contig: int


def filter_contigs(src: Path, dst: Path, *, min_length: int, prefix: str) -> ContigStats:
    """Copy ``src`` to ``dst`` keeping contigs of at least ``min_length`` bases,
    renamed ``<prefix>_contig<n>`` in the order kept; return the summary."""
    kept: list[str] = []
    for _header, seq in _records(src):
        if len(seq) >= min_length:
            kept.append(seq)
    lengths = sorted((len(s) for s in kept), reverse=True)
    with atomic_replace(dst) as fo:
        for i, seq in enumerate(kept, start=1):
            fo.write(f">{prefix}_contig{i}\n")
            for j in range(0, len(seq), 80):
                fo.write(seq[j : j + 80] + "\n")
    return ContigStats(
        n_contigs=len(lengths),
        total_length=sum(lengths),
        n50=_n50(lengths),
        largest_contig=lengths[0] if lengths else 0,
    )


def _records(path: Path):
    header: str | None = None
    parts: list[str] = []
    with open(path, encoding="utf-8") as fo:
        for line in fo:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:], []
            elif header is not None:
                parts.append(line.strip())
    if header is not None:
        yield header, "".join(parts)


def _n50(lengths_desc: list[int]) -> int:
    half = sum(lengths_desc) / 2
    running = 0
    for length in lengths_desc:
        running += length
        if running >= half:
            return length
    return 0
