"""Snippy SNP typer (per-genome calling + snippy-core)."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from ..core.process import run
from .base import SnpParams, SnpResult, SnpTyper


class SnippyTyper(SnpTyper):
    capabilities = ToolCapabilities(
        name="snippy",
        required_binaries=(
            BinarySpec("snippy", version_args=("--version",)),
            BinarySpec("snippy-core", version_args=("--version",)),
        ),
        recommended_max_genomes=1000,
        threads_param="--cpus",
    )
    requires_reference = True

    def call(
        self,
        genomes: Sequence[Path],
        reference: Path | None,
        out_dir: Path,
        params: SnpParams,
        logger: logging.Logger,
    ) -> SnpResult:
        genomes = list(genomes)
        if reference is None:
            reference = genomes[0]
        out_dir.mkdir(parents=True, exist_ok=True)

        sample_dirs = []
        for genome in genomes:
            if genome.resolve() == reference.resolve():
                continue
            sdir = out_dir / genome.stem
            run(
                [
                    "snippy", "--cpus", str(params.threads),
                    "--outdir", sdir, "--ref", reference.resolve(),
                    "--ctgs", genome.resolve(), "--force",
                ],
                logger=logger,
                log_prefix="snippy",
            )
            sample_dirs.append(sdir)

        core_prefix = out_dir / "core"
        run(
            ["snippy-core", "--ref", reference.resolve(), "--prefix", core_prefix, *sample_dirs],
            logger=logger,
            log_prefix="snippy-core",
        )
        core_aln = Path(str(core_prefix) + ".aln")
        if not core_aln.exists():
            raise WorkdirError("snippy-core did not produce a core alignment (.aln)")
        core_fasta = out_dir / "core_snp.fasta"
        core_fasta.write_text(core_aln.read_text())
        return SnpResult(core_snp_fasta=core_fasta, masked=False)
