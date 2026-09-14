"""Simulate sequencing reads from a genome FASTA for tests (no external tools).

Short reads: paired 150 bp fragments with a few substitutions; long reads:
5-20 kb fragments with more substitutions. Enough for a real assembler to
recover a synthetic genome; not a model of any instrument.
"""

from __future__ import annotations

import gzip
import random
from pathlib import Path

_COMP = str.maketrans("ACGT", "TGCA")


def _genome(fasta: Path) -> str:
    return "".join(
        line.strip()
        for line in fasta.read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    )


def _mutate(seq: str, rate: float, rng: random.Random) -> str:
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice("ACGT".replace(out[i], "") or "A")
    return "".join(out)


def _write(path: Path, reads: list[tuple[str, str]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as fo:
        for name, seq in reads:
            fo.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def simulate_paired_reads(
    fasta: Path,
    out_dir: Path,
    *,
    coverage: float = 30.0,
    read_length: int = 150,
    fragment: int = 500,
    fragment_sd: int = 100,
    error_rate: float = 0.005,
    seed: int = 1,
) -> tuple[Path, Path]:
    """Write ``<stem>_1.fastq.gz`` and ``<stem>_2.fastq.gz`` under ``out_dir``.

    Fragment sizes are drawn around ``fragment`` so that, as in a real library,
    only some pairs overlap and can be merged by an assembler's read merger.
    """
    rng = random.Random(seed)
    genome = _genome(fasta)
    n_pairs = int(len(genome) * coverage / (2 * read_length))
    r1, r2 = [], []
    for i in range(n_pairs):
        size = max(read_length, int(rng.gauss(fragment, fragment_sd)))
        start = rng.randint(0, max(0, len(genome) - size))
        frag = genome[start : start + size]
        fwd = _mutate(frag[:read_length], error_rate, rng)
        rev = _mutate(frag[-read_length:].translate(_COMP)[::-1], error_rate, rng)
        r1.append((f"read{i}/1", fwd))
        r2.append((f"read{i}/2", rev))
    out_dir.mkdir(parents=True, exist_ok=True)
    p1, p2 = out_dir / f"{fasta.stem}_1.fastq.gz", out_dir / f"{fasta.stem}_2.fastq.gz"
    _write(p1, r1)
    _write(p2, r2)
    return p1, p2


def simulate_long_reads(
    fasta: Path, out_dir: Path, *, coverage: float = 30.0, error_rate: float = 0.03, seed: int = 1
) -> Path:
    """Write ``<stem>.fastq.gz`` of 5-20 kb reads under ``out_dir``."""
    rng = random.Random(seed)
    genome = _genome(fasta)
    reads, total = [], 0
    while total < len(genome) * coverage:
        length = rng.randint(5000, 20000)
        start = rng.randint(0, max(0, len(genome) - 1000))
        seq = genome[start : start + length]
        if rng.random() < 0.5:
            seq = seq.translate(_COMP)[::-1]
        reads.append((f"read{len(reads)}", _mutate(seq, error_rate, rng)))
        total += len(seq)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{fasta.stem}.fastq.gz"
    _write(path, reads)
    return path
