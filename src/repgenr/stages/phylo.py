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

import logging
from collections.abc import Sequence
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..aligners.base import AlignParams
from ..aligners.base import registry as aligner_registry
from ..core.context import WorkdirContext
from ..core.contracts import TREE_NWK, parse_genome_filename
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


@dataclass
class PhyloDirs:
    """Output directories for the stateless tree-building core."""

    tree_dir: Path
    align_dir: Path
    snp_dir: Path
    scratch_dir: Path


@dataclass
class PhyloOutcome:
    tree: Path
    treebuilder: str
    versions: dict[str, str]
    outgroup_leaf: str | None


def build_tree(
    genomes: list[Path],
    outgroup_file: Path | None,
    outgroup_leaf: str | None,
    dirs: PhyloDirs,
    params: PhyloParams,
    logger: logging.Logger,
) -> PhyloOutcome:
    """Build a phylogeny from explicit inputs into ``dirs`` (stateless; no config).

    Selects (or auto-picks) the tree builder, derives the MSA from an aligner or
    a SNP typer when needed, roots by the outgroup and writes ``tree/tree.nwk``.
    The reference (for aligner/SNP sources) is resolved by basename against the
    genome set, so the core needs no working directory.
    """
    if not genomes:
        raise WorkdirError(
            "No genomes found for phylo. Run the genome (and derep) stages first."
        )

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
    dirs.tree_dir.mkdir(parents=True, exist_ok=True)

    if builder.input_kind == InputKind.GENOMES:
        inputs = list(genomes)
        if outgroup_file is not None:
            inputs.append(outgroup_file)
        logger.info(
            "Building tree with %s (alignment-free) over %d genomes",
            treebuilder, len(inputs),
        )
        tree = builder.build(inputs, dirs.tree_dir, tree_params, logger)
    else:
        msa, source_versions = _build_msa(genomes, outgroup_file, dirs, params, logger)
        versions = {**versions, **source_versions}
        logger.info("Building tree with %s from MSA %s", treebuilder, msa)
        tree = builder.build(msa, dirs.tree_dir, tree_params, logger)

    final = dirs.tree_dir / TREE_NWK
    if tree.resolve() != final.resolve():
        final.write_text(Path(tree).read_text())
    logger.info("Phylogenetic tree written to %s", final)
    return PhyloOutcome(
        tree=final, treebuilder=treebuilder, versions=versions, outgroup_leaf=outgroup_leaf
    )


@dataclass
class PhyloBuildParams:
    """Inputs for the stateless phylo step (explicit paths, no workdir)."""

    genomes_dir: Path
    out_dir: Path
    outgroup_dir: Path | None = None
    outgroup_accession: Path | None = None
    phylo: PhyloParams = field(default_factory=PhyloParams)


def phylo_build(params: PhyloBuildParams, logger: logging.Logger) -> Path:
    """Build a phylogeny from an explicit genomes directory (data-channel step)."""
    genomes = list_fasta(params.genomes_dir)
    if not genomes:
        raise WorkdirError(f"No genome FASTA files found in {params.genomes_dir}.")

    outgroup_file: Path | None = None
    outgroup_leaf: str | None = None
    if (
        not params.phylo.no_outgroup
        and params.outgroup_dir is not None
        and params.outgroup_accession is not None
    ):
        outgroup_file, outgroup_leaf = resolve_outgroup_files(
            params.outgroup_dir, params.outgroup_accession, logger
        )

    dirs = PhyloDirs(
        tree_dir=params.out_dir / "tree",
        align_dir=params.out_dir / "align",
        snp_dir=params.out_dir / "snp",
        scratch_dir=params.out_dir / "scratch",
    )
    outcome = build_tree(genomes, outgroup_file, outgroup_leaf, dirs, params.phylo, logger)
    return outcome.tree


def run(ctx: WorkdirContext, params: PhyloParams) -> Path:
    logger = ctx.logger
    genomes = _genome_set(ctx, params.all_genomes)
    if not genomes:
        raise WorkdirError(
            "No genomes found for phylo. Run the genome (and derep) stages first."
        )

    outgroup_file, outgroup_leaf = _resolve_outgroup(ctx, params.no_outgroup, logger)

    dirs = PhyloDirs(
        tree_dir=ctx.tree_dir,
        align_dir=ctx.align_dir,
        snp_dir=ctx.snp_dir,
        scratch_dir=ctx.scratch_dir,
    )
    outcome = build_tree(genomes, outgroup_file, outgroup_leaf, dirs, params, logger)

    is_msa = treebuilder_registry.create(outcome.treebuilder).input_kind == InputKind.MSA_FASTA
    ctx.config.record_stage(
        "phylo",
        tool=outcome.treebuilder,
        params={
            "requested_treebuilder": params.treebuilder,
            "msa_source": params.msa_source if is_msa else None,
            "aligner": params.aligner if params.msa_source == "aligner" else None,
            "snptyper": params.snptyper if params.msa_source == "snptype" else None,
            "all_genomes": params.all_genomes,
            "bootstrap": params.bootstrap,
            "outgroup": None if params.no_outgroup else outcome.outgroup_leaf,
        },
        tool_versions=outcome.versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return outcome.tree


def _genome_set(ctx: WorkdirContext, all_genomes: bool) -> list[Path]:
    source = ctx.genomes_dir if all_genomes else ctx.representatives_dir
    return list_fasta(source)


def list_fasta(source: Path) -> list[Path]:
    """Sorted FASTA files directly under ``source`` (empty if it does not exist)."""
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
    return resolve_outgroup_files(ctx.outgroup_dir, acc_file, logger)


def resolve_outgroup_files(
    outgroup_dir: Path, accession_file: Path, logger: logging.Logger
) -> tuple[Path | None, str | None]:
    """Resolve the outgroup genome file and leaf name from explicit paths."""
    if not accession_file.exists() or not outgroup_dir.exists():
        logger.warning("No outgroup found; proceeding without one")
        return None, None
    accession = accession_file.read_text().strip()
    if not accession:
        logger.warning("No outgroup accession recorded; proceeding without one")
        return None, None
    for f in sorted(outgroup_dir.iterdir()):
        if accession in f.name:
            logger.info("Using %s as outgroup", f.name)
            return f, f.stem
    logger.warning("Outgroup accession %s not found in %s", accession, outgroup_dir)
    return None, None


def _build_msa(
    genomes: list[Path],
    outgroup_file: Path | None,
    dirs: PhyloDirs,
    params: PhyloParams,
    logger: logging.Logger,
) -> tuple[Path, dict[str, str]]:
    inputs = list(genomes)
    if outgroup_file is not None:
        inputs.append(outgroup_file)

    if params.msa_source == "aligner":
        aligner = aligner_registry.create(params.aligner)
        versions = aligner.preflight()
        reference = _resolve_reference(params.reference, genomes, outgroup_file)
        _warn_divergence(params.aligner, inputs, logger)
        align_params = AlignParams(
            threads=params.threads, reference=reference, extra=dict(params.extra),
        )
        result = aligner.align(inputs, reference, dirs.align_dir, align_params, logger)
        return result.msa_fasta, versions

    if params.msa_source == "snptype":
        # Reuse the SNP typer's core-SNP alignment as the MSA source.
        from .snptype import SnptypeParams, snptype_core

        snp_params = SnptypeParams(
            tool=params.snptyper,
            threads=params.threads,
            reference=params.reference,
            all_genomes=params.all_genomes,
            mask=params.extra.get("mask", "none"),
        )
        snp_reference: Path | None = (
            _resolve_reference(params.reference, genomes, outgroup_file)
            if params.reference
            else None
        )
        snp_result, versions = snptype_core(
            genomes, snp_reference, dirs.snp_dir, dirs.scratch_dir / "snptype", snp_params, logger
        )
        return snp_result.core_snp_fasta, versions

    raise UserInputError(f"Unknown msa-source '{params.msa_source}' (aligner|snptype)")


def _resolve_reference(
    reference_name: str | None, genomes: Sequence[Path], outgroup_file: Path | None
) -> Path:
    """Resolve a named reference by basename against the genome set (else genomes[0])."""
    if reference_name:
        pool = list(genomes)
        if outgroup_file is not None:
            pool.append(outgroup_file)
        for p in pool:
            if p.name == reference_name:
                return p
        raise UserInputError(f"Reference genome not found: {reference_name}")
    return genomes[0]


def _taxonomic_spread(genomes: Sequence[Path]) -> tuple[int, int]:
    """Distinct (genera, species) among the inputs, read from the canonical
    ``Family_Genus_species_Accession.fasta`` filenames. Used to gauge divergence
    without the manifest, so it works in the shared-workdir and data-channel paths.
    """
    genera: set[str] = set()
    species: set[tuple[str, str]] = set()
    for g in genomes:
        _family, genus, sp, _acc = parse_genome_filename(Path(g).name)
        if genus:
            genera.add(genus.lower())
            species.add((genus.lower(), sp.lower()))
    return len(genera), len(species)


def _warn_divergence(aligner_name: str, genomes: Sequence[Path], logger) -> None:
    """Warn when a whole-genome aligner is run on a divergent (genus/family-level)
    set, where the shared collinear core shrinks and the alignment degrades.
    """
    n_genera, n_species = _taxonomic_spread(genomes)
    if aligner_name == "cactus" and n_species > 1:
        logger.warning(
            "Aligner 'cactus' (Minigraph-Cactus) targets same-species genomes, but the input "
            "spans %d species; it will likely drop divergent genomes from the graph. For "
            "genus/family-level data use an alignment-free tree builder (mashtree/sourmash).",
            n_species,
        )
    elif n_genera > 1:
        logger.warning(
            "Whole-genome aligner '%s' on a family-level set (%d genera): the shared collinear "
            "core shrinks sharply with divergence, so the alignment may be small or fragmentary. "
            "Consider an alignment-free tree builder (mashtree/sourmash), or loosen the aligner "
            "seeds (e.g. --aligner-arg kmer=15 for sibeliaz).",
            aligner_name, n_genera,
        )
    elif n_species > 1:
        logger.info(
            "Whole-genome aligner '%s' on a genus-level set (%d species): expect a reduced core "
            "alignment as divergence increases; alignment-free builders scale better.",
            aligner_name, n_species,
        )
