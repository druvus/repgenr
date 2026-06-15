"""Phylogenetics stage: compose an MSA source with a tree builder.

Three orthogonal choices:
  * genome set    -- dereplicated representatives (default) or all genomes
  * MSA source    -- a whole-genome aligner OR a SNP typer's core-SNP alignment
                     (skipped entirely for alignment-free tree builders)
  * tree builder  -- iqtree / fasttree / raxmlng (MSA) or mashtree / sourmash
                     (alignment-free)

Outgroup rooting is handled here once, regardless of the tools chosen.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..aligners.base import AlignParams
from ..aligners.base import registry as aligner_registry
from ..core.context import WorkdirContext
from ..core.contracts import TREE_NWK
from ..core.errors import UserInputError, WorkdirError
from ..core.plugins import auto_select, scale_warning
from ..treebuilders.base import InputKind, TreeParams
from ..treebuilders.base import registry as treebuilder_registry

_FASTA_SUFFIXES = (".fasta", ".fasta.gz", ".fa", ".fna", ".fas")


@dataclass
class PhyloParams:
    treebuilder: str = "iqtree"
    msa_source: str = "aligner"  # aligner | snptype
    aligner: str = "progressivemauve"
    snptyper: str = "simple"
    all_genomes: bool = False
    no_outgroup: bool = False
    threads: int = 16
    bootstrap: int = 0
    reference: str | None = None
    extra: dict = field(default_factory=dict)


def run(ctx: WorkdirContext, params: PhyloParams) -> Path:
    logger = ctx.logger
    genomes = _genome_set(ctx, params.all_genomes)
    if not genomes:
        raise WorkdirError(
            "No genomes found for phylo. Run the genome (and derep) stages first."
        )

    outgroup_file, outgroup_leaf = _resolve_outgroup(ctx, params.no_outgroup, logger)

    treebuilder = params.treebuilder
    if treebuilder == "auto":
        treebuilder = auto_select(treebuilder_registry, len(genomes)) or "iqtree"
        logger.info("Auto-selected tree builder '%s' for %d genomes", treebuilder, len(genomes))
    else:
        warn = scale_warning(treebuilder_registry, treebuilder, len(genomes))
        if warn:
            limit, alts = warn
            logger.warning(
                "Tree builder '%s' is tuned for <=%d genomes but you have %d; consider: %s",
                treebuilder, limit, len(genomes), ", ".join(alts) or "none",
            )

    builder = treebuilder_registry.create(treebuilder)
    versions = builder.preflight()

    tree_params = TreeParams(
        threads=params.threads,
        outgroup=None if params.no_outgroup else outgroup_leaf,
        bootstrap=params.bootstrap,
        extra=dict(params.extra),
    )
    ctx.tree_dir.mkdir(parents=True, exist_ok=True)

    if builder.input_kind == InputKind.GENOMES:
        inputs = list(genomes)
        if outgroup_file is not None:
            inputs.append(outgroup_file)
        logger.info(
            "Building tree with %s (alignment-free) over %d genomes",
            treebuilder, len(inputs),
        )
        tree = builder.build(inputs, ctx.tree_dir, tree_params, logger)
        source_tool = None
    else:
        msa, source_tool, source_versions = _build_msa(ctx, params, genomes, outgroup_file, logger)
        versions = {**versions, **source_versions}
        logger.info("Building tree with %s from MSA %s", treebuilder, msa)
        tree = builder.build(msa, ctx.tree_dir, tree_params, logger)

    final = ctx.tree_dir / TREE_NWK
    if tree.resolve() != final.resolve():
        final.write_text(Path(tree).read_text())

    ctx.config.record_stage(
        "phylo",
        tool=treebuilder,
        params={
            "requested_treebuilder": params.treebuilder,
            "msa_source": params.msa_source if builder.input_kind == InputKind.MSA_FASTA else None,
            "aligner": params.aligner if params.msa_source == "aligner" else None,
            "snptyper": params.snptyper if params.msa_source == "snptype" else None,
            "all_genomes": params.all_genomes,
            "bootstrap": params.bootstrap,
            "outgroup": None if params.no_outgroup else outgroup_leaf,
        },
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Phylogenetic tree written to %s", final)
    return final


def _genome_set(ctx: WorkdirContext, all_genomes: bool) -> list[Path]:
    source = ctx.genomes_dir if all_genomes else ctx.representatives_dir
    if not source.exists():
        return []
    return sorted(
        p for p in source.iterdir()
        if not p.name.startswith(".") and any(p.name.endswith(s) for s in _FASTA_SUFFIXES)
    )


def _resolve_outgroup(
    ctx: WorkdirContext, no_outgroup: bool, logger
) -> tuple[Path | None, str | None]:
    if no_outgroup:
        return None, None
    acc_file = ctx.workdir / "outgroup_accession.txt"
    if not acc_file.exists() or not ctx.outgroup_dir.exists():
        logger.warning("No outgroup found; proceeding without one")
        return None, None
    accession = acc_file.read_text().strip()
    for f in ctx.outgroup_dir.iterdir():
        if accession in f.name:
            logger.info("Using %s as outgroup", f.name)
            return f, f.stem
    logger.warning("Outgroup accession %s not found in %s", accession, ctx.outgroup_dir)
    return None, None


def _build_msa(ctx, params, genomes, outgroup_file, logger):
    inputs = list(genomes)
    if outgroup_file is not None:
        inputs.append(outgroup_file)

    if params.msa_source == "aligner":
        aligner = aligner_registry.create(params.aligner)
        versions = aligner.preflight()
        reference = _reference_path(ctx, params.reference, genomes)
        align_params = AlignParams(threads=params.threads, reference=reference)
        result = aligner.align(inputs, reference, ctx.align_dir, align_params, logger)
        return result.msa_fasta, params.aligner, versions

    if params.msa_source == "snptype":
        # Reuse the snptype stage output as the MSA source.
        from .snptype import SnptypeParams
        from .snptype import run as snptype_run

        snp_params = SnptypeParams(
            tool=params.snptyper,
            threads=params.threads,
            reference=params.reference,
            all_genomes=params.all_genomes,
            mask=params.extra.get("mask", "none"),
        )
        snp_result = snptype_run(ctx, snp_params)
        record = ctx.config.stages.get("snptype")
        versions = record.tool_versions if record else {}
        return snp_result.core_snp_fasta, params.snptyper, versions

    raise UserInputError(f"Unknown msa-source '{params.msa_source}' (aligner|snptype)")


def _reference_path(ctx, reference_name, genomes) -> Path:
    if reference_name:
        cand = ctx.representatives_dir / reference_name
        if cand.exists():
            return cand
        cand = ctx.genomes_dir / reference_name
        if cand.exists():
            return cand
        raise UserInputError(f"Reference genome not found: {reference_name}")
    return genomes[0]
