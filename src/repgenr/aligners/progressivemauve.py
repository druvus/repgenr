"""progressiveMauve aligner.

Reproduces the old ``phylo.py`` accurate path: align every genome to a single
reference with progressiveMauve, project each pairwise XMFA onto reference
coordinates (:mod:`repgenr.converters.xmfa_to_fasta`), then concatenate into one
reference-anchored MSA-FASTA. Aligning to a single reference keeps the work
linear in the number of genomes.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..converters.xmfa_to_fasta import xmfa_to_fasta
from ..core.binaries import BinarySpec
from ..core.errors import UserInputError
from ..core.plugins import ToolCapabilities
from ..core.process import run
from .base import Aligner, AlignParams, AlignResult


class ProgressiveMauveAligner(Aligner):
    capabilities = ToolCapabilities(
        name="progressivemauve",
        required_binaries=(BinarySpec("progressiveMauve", version_args=()),),
        recommended_max_genomes=500,
        threads_param=None,  # progressiveMauve is single-threaded per alignment
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
        if reference not in genomes:
            genomes = [reference, *genomes]
        if len(genomes) < 3:
            raise UserInputError(
                "progressiveMauve alignment needs at least 3 genomes (outgroup included)."
            )

        out_dir.mkdir(parents=True, exist_ok=True)
        xmfa_dir = out_dir / "xmfa"
        xmfa_dir.mkdir(exist_ok=True)
        ref_arg = str(reference.resolve())

        per_query_fastas: list[Path] = []
        for query in genomes:
            if query == reference:
                continue
            stem = query.stem
            xmfa = xmfa_dir / f"{stem}.xmfa"
            fa = xmfa_dir / f"{stem}.fa"
            run(
                ["progressiveMauve", "--output", xmfa, ref_arg, str(query.resolve())],
                logger=logger,
                log_prefix="progressivemauve",
            )
            xmfa_to_fasta(xmfa, ref_arg, 0, fa)
            per_query_fastas.append(fa)

        msa = out_dir / "msa.fasta"
        _concatenate(per_query_fastas, reference, msa)
        return AlignResult(msa_fasta=msa)


def _concatenate(per_query_fastas: list[Path], reference: Path, out_path: Path) -> None:
    """Write the reference row once, then each query row; leaf names are stems."""
    ref_stem = reference.stem
    written_ref = False
    with open(out_path, "w") as out:
        for fa in per_query_fastas:
            for name, seq in _read_fasta(fa):
                leaf = Path(name).stem
                if leaf == ref_stem:
                    if written_ref:
                        continue
                    written_ref = True
                out.write(f">{leaf}\n")
                for pos in range(0, len(seq), 80):
                    out.write(seq[pos : pos + 80] + "\n")


def _read_fasta(path: Path):
    name = None
    seq: list[str] = []
    with open(path) as fo:
        for line in fo:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)
