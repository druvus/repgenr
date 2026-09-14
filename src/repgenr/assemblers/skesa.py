"""SKESA assembler (Illumina; conservative, low memory)."""

from __future__ import annotations

import logging
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.plugins import ToolCapabilities
from .base import AssembleParams, Assembler, AssemblyResult, ReadSet


class SkesaAssembler(Assembler):
    capabilities = ToolCapabilities(
        name="skesa",
        container="quay.io/biocontainers/skesa:2.5.1--h077b44d_3",
        conda=("bioconda::skesa",),
        required_binaries=(BinarySpec("skesa", version_args=("--version",), min_version="2.4"),),
    )
    read_types = frozenset({"ILLUMINA"})

    def assemble(
        self, reads: ReadSet, out_dir: Path, params: AssembleParams, logger: logging.Logger
    ) -> AssemblyResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        contigs = out_dir / "contigs.fa"
        cmd: list[str | Path] = [
            "skesa",
            "--reads",
            ",".join(str(f) for f in reads.files),
            "--cores",
            str(params.threads),
            "--memory",
            str(params.memory_gb),
            "--contigs_out",
            contigs,
        ]
        run_tool(self.capabilities, cmd, logger=logger, log_prefix="skesa", cwd=out_dir)
        return AssemblyResult(contigs=contigs)
