"""Shovill assembler (Illumina paired-end; SPAdes with read correction and trimming)."""

from __future__ import annotations

import logging
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import UserInputError
from ..core.plugins import ToolCapabilities
from .base import AssembleParams, Assembler, AssemblyResult, ReadSet


class ShovillAssembler(Assembler):
    capabilities = ToolCapabilities(
        name="shovill",
        container="quay.io/biocontainers/shovill:1.4.2--hdfd78af_1",
        conda=("bioconda::shovill",),
        required_binaries=(BinarySpec("shovill", version_args=("--version",), min_version="1.1"),),
        default_params={"assembler": "spades"},
        accepted_extras=frozenset({"assembler", "depth", "minlen"}),
    )
    read_types = frozenset({"ILLUMINA"})
    layouts = frozenset({"PAIRED"})

    def assemble(
        self, reads: ReadSet, out_dir: Path, params: AssembleParams, logger: logging.Logger
    ) -> AssemblyResult:
        if len(reads.files) != 2:
            raise UserInputError(
                f"shovill needs paired-end reads; run {reads.run_accession} is "
                f"{reads.layout.lower()}. Use --assembler skesa for single-end Illumina."
            )
        out_dir.mkdir(parents=True, exist_ok=True)
        result_dir = out_dir / "shovill_out"
        defaults = self.capabilities.default_params
        cmd: list[str | Path] = [
            "shovill",
            "--R1",
            reads.files[0],
            "--R2",
            reads.files[1],
            "--outdir",
            result_dir,
            "--force",
            "--cpus",
            str(params.threads),
            "--ram",
            str(params.memory_gb),
            "--assembler",
            str(params.extra.get("assembler", defaults["assembler"])),
            "--minlen",
            str(params.extra.get("minlen", 200)),
        ]
        if "depth" in params.extra:
            cmd += ["--depth", str(params.extra["depth"])]
        run_tool(self.capabilities, cmd, logger=logger, log_prefix="shovill", cwd=out_dir)
        return AssemblyResult(contigs=result_dir / "contigs.fa")
