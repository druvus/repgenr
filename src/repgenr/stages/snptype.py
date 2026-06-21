"""SNP typing stage.

Selects a SNP typer, runs it against a reference (or reference-free),
optionally masks recombination with Gubbins, and writes the canonical SNP
outputs: ``snp/core_snp.fasta`` (+ optional VCF and SNP distance matrix). The
core-SNP alignment is both a standalone typing deliverable and an MSA source for
the phylo stage.

The compute is factored into :func:`snptype_core`, a stateless engine that takes
explicit input/output directories and never touches the run config or manifest.
The workdir-bound :func:`run` resolves paths from the context, calls the core and
records provenance; the data-channel phylo step reuses the core directly.
"""

from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import CORE_SNP_FASTA, list_fasta
from ..core.errors import UserInputError, WorkdirError
from ..snptypers.base import SnpParams, SnpResult
from ..snptypers.base import registry as snp_registry


@dataclass
class SnptypeParams:
    tool: str = "simple"
    threads: int = 16
    reference: str | None = None
    all_genomes: bool = False
    mask: str = "none"  # none | gubbins
    extra: dict = field(default_factory=dict)


def snptype_core(
    genomes: list[Path],
    reference: Path | None,
    snp_dir: Path,
    scratch: Path,
    params: SnptypeParams,
    logger: logging.Logger,
) -> tuple[SnpResult, dict[str, str]]:
    """Run a SNP typer over ``genomes`` into ``snp_dir`` (stateless; no config).

    ``reference`` is the already-resolved reference path (or None). For a
    reference-requiring typer a None reference falls back to ``genomes[0]``;
    reference-free typers ignore it. Returns the SNP result and tool versions.
    """
    if not genomes:
        raise WorkdirError("No genomes found for SNP typing. Run the genome (and derep) stages.")

    typer = snp_registry.create(params.tool)
    versions = typer.preflight()

    ref = None
    if typer.requires_reference:
        ref = reference if reference is not None else genomes[0]

    snp_dir.mkdir(parents=True, exist_ok=True)
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True, exist_ok=True)

    snp_params = SnpParams(
        threads=params.threads,
        reference=ref,
        mask=params.mask,
        extra=dict(params.extra),
    )
    logger.info("SNP typing %d genomes with %s", len(genomes), params.tool)
    result = typer.call(genomes, ref, scratch, snp_params, logger)

    core = snp_dir / CORE_SNP_FASTA
    masked = False
    if params.mask == "gubbins":
        from ..snptypers.gubbins import mask_recombination

        filtered = mask_recombination(result.core_snp_fasta, scratch / "gubbins", logger)
        shutil.copy2(filtered, core)
        masked = True
    elif params.mask not in ("none", ""):
        raise UserInputError(f"Unknown mask '{params.mask}' (none|gubbins)")
    else:
        shutil.copy2(result.core_snp_fasta, core)

    if result.vcf is not None:
        shutil.copy2(result.vcf, snp_dir / "variants.vcf")
    if result.snp_distance_matrix is not None:
        shutil.copy2(result.snp_distance_matrix, snp_dir / "snp_distance_matrix.tsv")

    return (
        SnpResult(
            core_snp_fasta=core,
            vcf=(snp_dir / "variants.vcf") if result.vcf else None,
            snp_distance_matrix=(snp_dir / "snp_distance_matrix.tsv")
            if result.snp_distance_matrix
            else None,
            masked=masked,
        ),
        versions,
    )


def run(ctx: WorkdirContext, params: SnptypeParams) -> SnpResult:
    logger = ctx.logger
    genomes = _genome_set(ctx, params.all_genomes)
    if not genomes:
        raise WorkdirError("No genomes found for SNP typing. Run the genome (and derep) stages.")

    reference = None
    if snp_registry.create(params.tool).requires_reference:
        reference = _reference_path(ctx, params.reference, genomes)

    result, versions = snptype_core(
        genomes, reference, ctx.snp_dir, ctx.scratch_dir / "snptype", params, logger
    )

    ctx.config.record_stage(
        "snptype",
        tool=params.tool,
        params={
            "all_genomes": params.all_genomes,
            "reference": reference.name if reference else None,
            "mask": params.mask,
        },
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("SNP typing complete: %s", result.core_snp_fasta)
    return result


def _genome_set(ctx: WorkdirContext, all_genomes: bool) -> list[Path]:
    source = ctx.genomes_dir if all_genomes else ctx.representatives_dir
    return list_fasta(source)


def _reference_path(ctx, reference_name, genomes) -> Path:
    if reference_name:
        for base in (ctx.representatives_dir, ctx.genomes_dir):
            cand = base / reference_name
            if cand.exists():
                return cand
        raise UserInputError(f"Reference genome not found: {reference_name}")
    return genomes[0]
