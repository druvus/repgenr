"""kSNP4 SNP typer (alignment-free k-mer SNP discovery)."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from .base import SnpParams, SnpResult, SnpTyper


class Ksnp4Typer(SnpTyper):
    capabilities = ToolCapabilities(
        name="ksnp",
        conda=("bioconda::ksnp4",),
        required_binaries=(BinarySpec("kSNP4", version_args=()),),
        recommended_max_genomes=2000,
        threads_param="-CPU",
    )
    requires_reference = False  # kSNP is reference-free

    def call(
        self,
        genomes: Sequence[Path],
        reference: Path | None,
        out_dir: Path,
        params: SnpParams,
        logger: logging.Logger,
    ) -> SnpResult:
        genomes = list(genomes)
        out_dir.mkdir(parents=True, exist_ok=True)

        # kSNP wants a tab-separated "path<TAB>name" input list.
        in_list = out_dir / "in_list.txt"
        with open(in_list, "w") as fo:
            for genome in genomes:
                fo.write(f"{genome.resolve()}\t{genome.stem}\n")

        results = out_dir / "ksnp_out"
        kmer = str(params.extra.get("kmer", 21))
        run_tool(self.capabilities, 
            ["kSNP4", "-in", in_list, "-outdir", results, "-k", kmer, "-CPU", str(params.threads)],
            logger=logger,
            log_prefix="ksnp",
        )
        matrix = results / "core_SNPs_matrix.fasta"
        if not matrix.exists():
            matrix = results / "SNPs_all_matrix.fasta"
        if not matrix.exists():
            raise WorkdirError("kSNP4 did not produce a SNP matrix FASTA")
        core_fasta = out_dir / "core_snp.fasta"
        core_fasta.write_text(matrix.read_text())
        return SnpResult(core_snp_fasta=core_fasta, masked=False)
