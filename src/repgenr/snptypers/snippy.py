"""Snippy SNP typer (per-genome calling + snippy-core)."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.contracts import strip_fasta_suffix
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from ..core.process import warn_argv_bytes
from .base import SnpParams, SnpResult, SnpTyper


class SnippyTyper(SnpTyper):
    capabilities = ToolCapabilities(
        name="snippy",
        container="quay.io/biocontainers/snippy:4.6.0--hdfd78af_6",
        conda=("bioconda::snippy",),
        required_binaries=(
            BinarySpec("snippy", version_args=("--version",), min_version="4.6"),
            BinarySpec("snippy-core", version_args=("--version",), min_version="4.6"),
        ),
        recommended_max_genomes=1000,
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
            sdir = out_dir / strip_fasta_suffix(genome.name)
            run_tool(
                self.capabilities,
                [
                    "snippy",
                    "--cpus",
                    str(params.threads),
                    "--outdir",
                    sdir,
                    "--ref",
                    reference.resolve(),
                    "--ctgs",
                    genome.resolve(),
                    "--force",
                ],
                logger=logger,
                log_prefix="snippy",
            )
            sample_dirs.append(sdir)

        core_prefix = out_dir / "core"
        core_cmd: list[str | Path] = [
            "snippy-core",
            "--ref",
            reference.resolve(),
            "--prefix",
            core_prefix,
            *sample_dirs,
        ]
        warn_argv_bytes("snippy-core", core_cmd, logger)
        run_tool(
            self.capabilities,
            core_cmd,
            logger=logger,
            log_prefix="snippy-core",
        )
        core_aln = Path(str(core_prefix) + ".aln")
        if not core_aln.exists():
            raise WorkdirError("snippy-core did not produce a core alignment (.aln)")
        core_fasta = out_dir / "core_snp.fasta"
        _copy_naming_reference(core_aln, core_fasta, strip_fasta_suffix(reference.name))

        full_aln = Path(str(core_prefix) + ".full.aln")
        full: Path | None = None
        if full_aln.exists():
            full = out_dir / "full_alignment.fasta"
            _copy_naming_reference(full_aln, full, strip_fasta_suffix(reference.name))

        return SnpResult(core_snp_fasta=core_fasta, masked=False, full_alignment=full)


def _copy_naming_reference(src: Path, dst: Path, reference_name: str) -> None:
    """Copy a snippy-core alignment, renaming its 'Reference' record.

    snippy-core labels the reference 'Reference'; every other typer names it
    by its genome, as the tree leaves and tree2tax expect. Streams line by
    line, since the whole-genome alignment is the largest file snippy-core
    writes.
    """
    with open(src, encoding="utf-8") as fi, open(dst, "w", encoding="utf-8") as fo:
        for line in fi:
            if line.rstrip("\r\n") == ">Reference":
                line = f">{reference_name}\n"
            fo.write(line)
