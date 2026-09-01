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
from ..core.contracts import CLUSTERS_TSV, CORE_SNP_FASTA, atomic_path, list_fasta
from ..core.errors import UserInputError, WorkdirError
from ..core.integrity import check_genome_completeness, check_representatives_consistency
from ..core.plugins import scale_warning, warn_unconsumed_extras
from ..snptypers.base import SnpParams, SnpResult
from ..snptypers.base import registry as snp_registry


@dataclass
class SnptypeParams:
    tool: str = "simple"
    threads: int = 16
    reference: str | None = None
    all_genomes: bool = False
    mask: str = "none"  # none | gubbins
    allow_incomplete: bool = False
    extra: dict = field(default_factory=dict)


def _check_inputs(ctx, all_genomes: bool, allow_incomplete: bool, logger) -> None:
    if all_genomes:
        check_genome_completeness(
            ctx.genomes_dir, ctx.workdir, logger=logger, allow_incomplete=allow_incomplete
        )
    else:
        check_representatives_consistency(
            ctx.representatives_dir, ctx.derep_dir / CLUSTERS_TSV,
            logger=logger, allow_incomplete=allow_incomplete,
        )


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

    warn = scale_warning(snp_registry, params.tool, len(genomes))
    if warn:
        limit, alts = warn
        logger.warning(
            "SNP typer '%s' is tuned for <=%d genomes but you have %d; consider: %s",
            params.tool, limit, len(genomes), ", ".join(alts) or "none",
        )
    typer = snp_registry.create(params.tool)
    warn_unconsumed_extras(typer.capabilities, params.extra, logger, family="SNP typer")
    versions = typer.preflight()

    ref = None
    if typer.requires_reference:
        if reference is not None:
            ref = reference
        else:
            ref = genomes[0]
            logger.warning(
                "No --reference given; SNP calling against the alphabetically "
                "first genome '%s'. Reference-private errors bias every SNP "
                "distance; pass --reference to choose deliberately.",
                ref.name,
            )

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
    if params.mask not in ("none", ""):
        from ..maskers.base import registry as masker_registry

        masker = masker_registry.create(params.mask)
        versions.update(masker.preflight())
        filtered = masker.mask(
            result.core_snp_fasta, scratch / params.mask, logger
        )
        with atomic_path(core) as tmp:
            shutil.copy2(filtered, tmp)
        masked = True
    else:
        with atomic_path(core) as tmp:
            shutil.copy2(result.core_snp_fasta, tmp)

    if result.vcf is not None:
        with atomic_path(snp_dir / "variants.vcf") as tmp:
            shutil.copy2(result.vcf, tmp)
    if result.snp_distance_matrix is not None:
        with atomic_path(snp_dir / "snp_distance_matrix.tsv") as tmp:
            shutil.copy2(result.snp_distance_matrix, tmp)

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
    _check_inputs(ctx, params.all_genomes, params.allow_incomplete, logger)
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
        # Resolved against the workdir genome dirs only: a path component would
        # let the lookup escape them.
        if Path(reference_name).name != reference_name:
            raise UserInputError(
                f"--reference must be a genome file basename, not a path: {reference_name}"
            )
        for base in (ctx.representatives_dir, ctx.genomes_dir):
            cand = base / reference_name
            if cand.exists():
                return cand
        raise UserInputError(f"Reference genome not found: {reference_name}")
    return genomes[0]
