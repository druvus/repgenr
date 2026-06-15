"""Minigraph-Cactus (cactus-pangenome) aligner.

Builds a pangenome from same-species genomes. We run ``cactus-pangenome`` with a
generated seqFile, then normalize its HAL output to a reference-anchored
MSA-FASTA via HAL -> MAF -> FASTA. Cactus is resource heavy and drives its own
Toil job store; the Nextflow process that wraps this stage requests a high
resource label.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..converters.hal_to_maf import hal_to_maf
from ..converters.maf_to_fasta import maf_to_fasta
from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from .base import Aligner, AlignParams, AlignResult


class CactusAligner(Aligner):
    capabilities = ToolCapabilities(
        name="cactus",
        container="quay.io/comparative-genomics-toolkit/cactus:v2.9.3",
        required_binaries=(
            BinarySpec("cactus-pangenome", version_args=("--version",)),
            BinarySpec("hal2maf", version_args=()),
        ),
        recommended_max_genomes=2000,
        threads_param=None,  # Toil manages its own parallelism
    )

    def align(
        self,
        genomes: Sequence[Path],
        reference: Path | None,
        out_dir: Path,
        params: AlignParams,
        logger: logging.Logger,
    ) -> AlignResult:
        genomes = list(genomes)
        if reference is None:
            reference = genomes[0]
        out_dir.mkdir(parents=True, exist_ok=True)

        seqfile = out_dir / "seqfile.txt"
        with open(seqfile, "w") as fo:
            for g in genomes:
                fo.write(f"{_sample_name(g)}\t{g.resolve()}\n")

        job_store = out_dir / "jobstore"
        results = out_dir / "cactus_out"
        ref_name = _sample_name(reference)
        run_tool(self.capabilities, 
            [
                "cactus-pangenome",
                job_store,
                seqfile,
                "--outDir", results,
                "--outName", "pangenome",
                "--reference", ref_name,
            ],
            logger=logger,
            log_prefix="cactus",
        )

        hal = _find_hal(results)
        if hal is None:
            raise WorkdirError(f"cactus-pangenome produced no HAL under {results}")

        maf = out_dir / "pangenome.maf"
        hal_to_maf(hal, ref_name, maf, logger, caps=self.capabilities)
        msa = out_dir / "msa.fasta"
        maf_to_fasta(maf, ref_name, msa)
        return AlignResult(msa_fasta=msa, native_format=hal)


def _sample_name(genome: Path) -> str:
    """Cactus sample names must avoid '.' (reserved for haplotype suffixes)."""
    return genome.stem.replace(".", "_")


def _find_hal(results: Path) -> Path | None:
    candidates = sorted(results.rglob("*.hal"))
    return candidates[0] if candidates else None
