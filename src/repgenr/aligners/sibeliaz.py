"""SibeliaZ aligner.

SibeliaZ builds locally collinear blocks and emits a MAF for moderately
divergent genomes (up to ~0.09 substitutions/site). The MAF is projected onto
reference coordinates by :mod:`repgenr.converters.maf_to_fasta`.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..converters.maf_to_fasta import maf_to_fasta
from ..core.binaries import BinarySpec
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from ..core.process import run
from .base import Aligner, AlignParams, AlignResult


class SibeliazAligner(Aligner):
    capabilities = ToolCapabilities(
        name="sibeliaz",
        required_binaries=(BinarySpec("sibeliaz", version_args=("-v",)),),
        recommended_max_genomes=2000,
        threads_param="-t",
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

        # Default mode runs the block alignment step, producing alignment.maf.
        # (The -n flag would emit only block coordinates, no MAF.)
        run(
            [
                "sibeliaz",
                "-t", str(params.threads),
                "-o", out_dir,
                *[str(g.resolve()) for g in genomes],
            ],
            logger=logger,
            log_prefix="sibeliaz",
        )
        maf = out_dir / "alignment.maf"
        if not maf.exists():
            # SibeliaZ writes blocks_coords.gff + alignment.maf; locate the MAF
            candidates = sorted(out_dir.rglob("*.maf"))
            if not candidates:
                raise WorkdirError("SibeliaZ did not produce a MAF file")
            maf = candidates[0]

        # SibeliaZ's MAF uses sequence/contig IDs (FASTA header first token), not
        # genome filenames; build the seqid -> genome-stem map for the converter.
        name_map = _build_seqid_map(genomes)
        msa = out_dir / "msa.fasta"
        maf_to_fasta(maf, reference.stem, msa, name_map=name_map)
        return AlignResult(msa_fasta=msa, native_format=maf)


def _build_seqid_map(genomes) -> dict[str, str]:
    """Map each FASTA sequence ID (header first token) to its genome stem."""
    name_map: dict[str, str] = {}
    for genome in genomes:
        stem = genome.stem
        with open(genome) as fo:
            for line in fo:
                if line.startswith(">"):
                    seqid = line[1:].split()[0]
                    name_map[seqid] = stem
    return name_map
