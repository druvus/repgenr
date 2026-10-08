"""ParSNP SNP typer (Harvest suite)."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.contracts import rename_fasta_records, strip_fasta_suffix
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from ..core.process import link_or_copy
from .base import SnpParams, SnpResult, SnpTyper


class ParsnpTyper(SnpTyper):
    capabilities = ToolCapabilities(
        name="parsnp",
        conda=("bioconda::parsnp", "bioconda::harvesttools"),
        required_binaries=(
            BinarySpec("parsnp", version_args=("--version",), min_version="2.0"),
            BinarySpec("harvesttools", version_args=("--version",), min_version="1.3"),
        ),
        recommended_max_genomes=2000,
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
        # Hardlinked where the file system allows, as the representatives are:
        # copies would hold a second set of every genome in scratch.
        gdir = out_dir / "input_genomes"
        gdir.mkdir(exist_ok=True)
        for genome in genomes:
            if genome.resolve() == reference.resolve():
                continue
            link_or_copy(genome, gdir / genome.name)

        results = out_dir / "parsnp_out"
        cmd: list[str | Path] = [
            "parsnp",
            "-r",
            reference.resolve(),
            "-d",
            gdir,
            "-o",
            results,
            "-p",
            str(params.threads),
        ]
        run_tool(
            self.capabilities,
            cmd,
            logger=logger,
            log_prefix="parsnp",
        )
        ggr = results / "parsnp.ggr"
        if not ggr.exists():
            raise WorkdirError("ParSNP did not produce parsnp.ggr")

        names = _record_names(genomes, reference)
        core_fasta = out_dir / "core_snp.fasta"
        raw_core = out_dir / "harvest_core_snp.fasta"
        run_tool(
            self.capabilities,
            ["harvesttools", "-i", ggr, "-S", raw_core],
            logger=logger,
            log_prefix="harvesttools",
        )
        rename_fasta_records(raw_core, core_fasta, names)

        full_fasta = out_dir / "full_alignment.fasta"
        raw_full = out_dir / "harvest_full_alignment.fasta"
        run_tool(
            self.capabilities,
            ["harvesttools", "-i", ggr, "-M", raw_full],
            logger=logger,
            log_prefix="harvesttools",
        )
        rename_fasta_records(raw_full, full_fasta, names)

        return SnpResult(core_snp_fasta=core_fasta, masked=False, full_alignment=full_fasta)


def _record_names(genomes: Sequence[Path], reference: Path) -> dict[str, str]:
    """harvesttools record name -> genome stem.

    harvesttools names a record by its file name, and the reference by its
    file name plus '.ref'. Every other typer names records by genome stem
    (the file name without its FASTA suffix and .gz),
    which the tree leaves, the outgroup lookup in tree2tax and the masker's
    outgroup exclusion all expect. The reference is passed resolved, so its
    target's name is mapped as well.
    """
    names = {genome.name: strip_fasta_suffix(genome.name) for genome in genomes}
    for name in {reference.name, reference.resolve().name}:
        names[f"{name}.ref"] = strip_fasta_suffix(reference.name)
    return names
