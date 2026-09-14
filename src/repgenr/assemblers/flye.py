"""Flye assembler (long reads: Oxford Nanopore and PacBio)."""

from __future__ import annotations

import logging
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import UserInputError
from ..core.plugins import ToolCapabilities
from .base import AssembleParams, Assembler, AssemblyResult, ReadSet

_MODES = ("nano-raw", "nano-hq", "pacbio-raw", "pacbio-hifi")
# Instruments whose PacBio output is HiFi (CCS) rather than CLR.
_HIFI_INSTRUMENTS = ("sequel ii", "revio")


class FlyeAssembler(Assembler):
    capabilities = ToolCapabilities(
        name="flye",
        container="quay.io/biocontainers/flye:2.9.6--py312h734f728_1",
        conda=("bioconda::flye",),
        required_binaries=(BinarySpec("flye", version_args=("--version",), min_version="2.9"),),
        accepted_extras=frozenset({"mode", "genome_size"}),
        # Flye sizes its own memory; the hint is not a flag it takes.
        ignored_params=frozenset({"memory_gb"}),
    )
    read_types = frozenset({"OXFORD_NANOPORE", "PACBIO_SMRT"})
    layouts = frozenset({"SINGLE"})

    def assemble(
        self, reads: ReadSet, out_dir: Path, params: AssembleParams, logger: logging.Logger
    ) -> AssemblyResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        result_dir = out_dir / "flye_out"
        mode = str(params.extra.get("mode") or _default_mode(reads))
        if mode not in _MODES:
            raise UserInputError(f"--tool-arg mode={mode!r}: choose one of {', '.join(_MODES)}.")
        cmd: list[str | Path] = [
            "flye",
            f"--{mode}",
            *reads.files,
            "-o",
            result_dir,
            "--threads",
            str(params.threads),
        ]
        if "genome_size" in params.extra:
            cmd += ["--genome-size", str(params.extra["genome_size"])]
        run_tool(self.capabilities, cmd, logger=logger, log_prefix="flye", cwd=out_dir)
        return AssemblyResult(contigs=result_dir / "assembly.fasta")


def _default_mode(reads: ReadSet) -> str:
    """Flye's read mode from the platform and instrument: R10-era ONT as
    nano-hq, Sequel II and Revio output as HiFi, older PacBio as raw CLR."""
    if reads.platform == "OXFORD_NANOPORE":
        return "nano-hq"
    model = reads.instrument_model.lower()
    return "pacbio-hifi" if any(tag in model for tag in _HIFI_INSTRUMENTS) else "pacbio-raw"
