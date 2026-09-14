"""The assemble stage with a real assembler on simulated reads (skipped without it)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from readsim import simulate_paired_reads  # noqa: E402

from repgenr.core.context import WorkdirContext  # noqa: E402
from repgenr.core.contracts import (  # noqa: E402
    READS_TSV,
    ReadRow,
    read_assembly_stats,
    write_reads,
)
from repgenr.stages.assemble import AssembleParams, run  # noqa: E402


def _synthetic_genome(path: Path, length: int = 60_000, seed: int = 7) -> None:
    import random

    rng = random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(length))
    path.write_text(f">syn\n{seq}\n", encoding="utf-8")


@pytest.mark.requires_binary("skesa")
def test_skesa_recovers_a_synthetic_genome(workdir: Path, tmp_path: Path) -> None:
    genome = tmp_path / "syn.fasta"
    _synthetic_genome(genome)
    r1, r2 = simulate_paired_reads(genome, tmp_path / "reads", coverage=40)
    ctx = WorkdirContext(workdir, create=True)
    write_reads(
        workdir / READS_TSV,
        [
            ReadRow(
                "SRRSYN",
                "SAM1",
                "PRJ1",
                "Synthetic organism",
                "1",
                "ILLUMINA",
                "MiSeq",
                "PAIRED",
                60_000 * 40,
                16_000,
                "Synfam",
                "Syngen",
                "syn",
                (str(r1), str(r2)),
                (),
                (r1.stat().st_size, r2.stat().st_size),
            )
        ],
    )
    assert run(ctx, AssembleParams(assembler="skesa", threads=2, jobs=1)) == 1
    stats = read_assembly_stats(workdir / "assembly_stats.tsv")[0]
    assert stats.total_length > 50_000 and stats.n_contigs < 20
