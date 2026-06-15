"""ParSNP SNP typer (Harvest suite)."""

from __future__ import annotations

import logging
import shutil
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from .base import SnpParams, SnpResult, SnpTyper


class ParsnpTyper(SnpTyper):
    capabilities = ToolCapabilities(
        name="parsnp",
        conda=("bioconda::parsnp", "bioconda::harvesttools"),
        required_binaries=(
            BinarySpec("parsnp", version_args=("--version",)),
            BinarySpec("harvesttools", version_args=("--version",)),
        ),
        recommended_max_genomes=2000,
        threads_param="-p",
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

        # ParSNP wants a directory of query genomes (excluding the reference).
        gdir = out_dir / "input_genomes"
        gdir.mkdir(exist_ok=True)
        for genome in genomes:
            if genome.resolve() == reference.resolve():
                continue
            shutil.copy2(genome, gdir / genome.name)

        results = out_dir / "parsnp_out"
        cmd: list[str | Path] = [
            "parsnp", "-r", reference.resolve(), "-d", gdir,
            "-o", results, "-p", str(params.threads),
        ]
        run_tool(self.capabilities, 
            cmd,
            logger=logger,
            log_prefix="parsnp",
        )
        ggr = results / "parsnp.ggr"
        if not ggr.exists():
            raise WorkdirError("ParSNP did not produce parsnp.ggr")

        core_fasta = out_dir / "core_snp.fasta"
        run_tool(self.capabilities, 
            ["harvesttools", "-i", ggr, "-S", core_fasta],
            logger=logger,
            log_prefix="harvesttools",
        )
        return SnpResult(core_snp_fasta=core_fasta, masked=False)
