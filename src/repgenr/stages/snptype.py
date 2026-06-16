"""SNP typing stage.

Selects a SNP typer, runs it against a reference (or reference-free),
optionally masks recombination with Gubbins, and writes the canonical SNP
outputs: ``snp/core_snp.fasta`` (+ optional VCF and SNP distance matrix). The
core-SNP alignment is both a standalone typing deliverable and an MSA source for
the phylo stage.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import CORE_SNP_FASTA
from ..core.errors import UserInputError, WorkdirError
from ..snptypers.base import SnpParams, SnpResult
from ..snptypers.base import registry as snp_registry

_FASTA_SUFFIXES = (".fasta", ".fasta.gz", ".fa", ".fna", ".fas")


@dataclass
class SnptypeParams:
    tool: str = "simple"
    threads: int = 16
    reference: str | None = None
    all_genomes: bool = False
    mask: str = "none"  # none | gubbins
    extra: dict = field(default_factory=dict)


def run(ctx: WorkdirContext, params: SnptypeParams) -> SnpResult:
    logger = ctx.logger
    genomes = _genome_set(ctx, params.all_genomes)
    if not genomes:
        raise WorkdirError("No genomes found for SNP typing. Run the genome (and derep) stages.")

    typer = snp_registry.create(params.tool)
    versions = typer.preflight()

    reference = None
    if typer.requires_reference:
        reference = _reference_path(ctx, params.reference, genomes)

    ctx.snp_dir.mkdir(parents=True, exist_ok=True)
    scratch = ctx.scratch_dir / "snptype"
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True, exist_ok=True)

    snp_params = SnpParams(
        threads=params.threads,
        reference=reference,
        mask=params.mask,
        extra=dict(params.extra),
    )
    logger.info("SNP typing %d genomes with %s", len(genomes), params.tool)
    result = typer.call(genomes, reference, scratch, snp_params, logger)

    core = ctx.snp_dir / CORE_SNP_FASTA
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
        shutil.copy2(result.vcf, ctx.snp_dir / "variants.vcf")
    if result.snp_distance_matrix is not None:
        shutil.copy2(result.snp_distance_matrix, ctx.snp_dir / "snp_distance_matrix.tsv")

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
    logger.info("SNP typing complete: %s", core)

    return SnpResult(
        core_snp_fasta=core,
        vcf=(ctx.snp_dir / "variants.vcf") if result.vcf else None,
        snp_distance_matrix=(ctx.snp_dir / "snp_distance_matrix.tsv")
        if result.snp_distance_matrix
        else None,
        masked=masked,
    )


def _genome_set(ctx: WorkdirContext, all_genomes: bool) -> list[Path]:
    source = ctx.genomes_dir if all_genomes else ctx.representatives_dir
    if not source.exists():
        return []
    return sorted(
        p for p in source.iterdir()
        if not p.name.startswith(".") and any(p.name.endswith(s) for s in _FASTA_SUFFIXES)
    )


def _reference_path(ctx, reference_name, genomes) -> Path:
    if reference_name:
        for base in (ctx.representatives_dir, ctx.genomes_dir):
            cand = base / reference_name
            if cand.exists():
                return cand
        raise UserInputError(f"Reference genome not found: {reference_name}")
    return genomes[0]
