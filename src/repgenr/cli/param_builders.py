"""Single construction point for each pipeline stage's parameters.

Both the manual per-stage commands and the ``run`` orchestrator build stage
parameters through these functions, so the two entry points cannot drift:
identical CLI intent produces identical params objects and therefore identical
resume fingerprints. Defaults live on the params dataclasses -- a builder
argument left unset is simply not passed, so a default changed on the
dataclass takes effect everywhere at once. Shared validation
(tool/threshold checks) also lives here so both entry points reject the same
inputs the same way.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from .base import _require_choice, _require_unit_interval

if TYPE_CHECKING:
    from ..stages.dereplicate import DereplicateParams
    from ..stages.genome import GenomeParams
    from ..stages.ingest import IngestParams
    from ..stages.metadata import MetadataParams
    from ..stages.phylo import PhyloParams
    from ..stages.tree2tax import Tree2taxParams
    from ..stages.vgenome import VgenomeParams
    from ..stages.vmetadata import VmetadataParams

# Sentinel for "caller did not set this": the argument is not passed to the
# dataclass, whose own default then applies.
_UNSET: Any = object()

# Closed choice sets shared with `run`, which validates the same values under
# its own flag names (--metadata-source, --viral-source).
METADATA_DATASETS = frozenset({"all", "rep"})
METADATA_LEVELS = frozenset({"family", "genus", "species"})
METADATA_SOURCES = frozenset({"tsv", "api"})
VIRAL_SOURCES = frozenset({"ncbi_virus", "bvbrc"})
LENGTH_METHODS = frozenset({"median_of_medians", "mean"})
DEREP_STOCK_ACTIONS = frozenset({"list", "pack", "unpack", "delete"})


def _build(cls: type, **kwargs: Any):
    return cls(**{k: v for k, v in kwargs.items() if v is not _UNSET})


def require_mask(mask: str) -> None:
    """``--mask`` is ``none`` or a registered masker."""
    from ..maskers.base import registry as _mask_registry

    _require_choice(mask, {"none", *_mask_registry.names()}, "--mask")


def metadata_params(
    *,
    dataset: str,
    level: str,
    source: Any = _UNSET,
    release: Any = _UNSET,
    version: Any = _UNSET,
    target_family: Any = _UNSET,
    target_genus: Any = _UNSET,
    target_species: Any = _UNSET,
    outgroup_accession: Any = _UNSET,
    metadata_path: Any = _UNSET,
    nodownload: Any = _UNSET,
    limit: Any = _UNSET,
) -> MetadataParams:
    from ..stages.metadata import MetadataParams

    _require_choice(dataset, METADATA_DATASETS, "--dataset")
    _require_choice(level, METADATA_LEVELS, "--level")
    if source is not _UNSET:
        _require_choice(source, METADATA_SOURCES, "--source")
    return _build(
        MetadataParams,
        dataset=dataset,
        level=level,
        source=source,
        release=release,
        version=version,
        target_family=target_family,
        target_genus=target_genus,
        target_species=target_species,
        outgroup_accession=outgroup_accession,
        metadata_path=metadata_path,
        nodownload=nodownload,
        limit=limit,
    )


def genome_params(*, accession_list_only: Any = _UNSET, keep_files: Any = _UNSET) -> GenomeParams:
    from ..stages.genome import GenomeParams

    return _build(GenomeParams, accession_list_only=accession_list_only, keep_files=keep_files)


def ingest_params(
    *,
    genomes_dir: str,
    selection: Any = _UNSET,
    outgroup: Any = _UNSET,
    copy: Any = _UNSET,
) -> IngestParams:
    from ..stages.ingest import IngestParams

    return _build(
        IngestParams, genomes_dir=genomes_dir, selection=selection, outgroup=outgroup, copy=copy
    )


def vmetadata_params(
    *,
    target: Any = _UNSET,
    source: Any = _UNSET,
    filter: Any = _UNSET,
    list_targets: Any = _UNSET,
    host: Any = _UNSET,
    complete_only: Any = _UNSET,
    released_after: Any = _UNSET,
) -> VmetadataParams:
    from ..stages.vmetadata import VmetadataParams

    if source is not _UNSET:
        _require_choice(source, VIRAL_SOURCES, "--source")
    return _build(
        VmetadataParams,
        target=target,
        source=source,
        filter=filter,
        list_targets=list_targets,
        host=host,
        complete_only=complete_only,
        released_after=released_after,
    )


def vgenome_params(
    *,
    target_genus: Any = _UNSET,
    target_species: Any = _UNSET,
    target_serotype: Any = _UNSET,
    target_custom: Any = _UNSET,
    length_all: Any = _UNSET,
    length_deviation: Any = _UNSET,
    length_method: Any = _UNSET,
    length_range: Any = _UNSET,
    discard: Any = _UNSET,
    no_outgroup: Any = _UNSET,
    group_segments: Any = _UNSET,
    outgroup_candidates_taxid_min_genomes: Any = _UNSET,
    outgroup_treebuilder: Any = _UNSET,
    glance: Any = _UNSET,
    print_fasta_headers: Any = _UNSET,
    ignore_duplicates: Any = _UNSET,
    keep_files: Any = _UNSET,
) -> VgenomeParams:
    from ..stages.vgenome import VgenomeParams
    from ..viral._outgroup import distance_matrix_builders

    if length_method is not _UNSET:
        _require_choice(length_method, LENGTH_METHODS, "--length-method")
    if outgroup_treebuilder is not _UNSET:
        _require_choice(
            outgroup_treebuilder, set(distance_matrix_builders()), "--outgroup-treebuilder"
        )
    return _build(
        VgenomeParams,
        target_genus=target_genus,
        target_species=target_species,
        target_serotype=target_serotype,
        target_custom=target_custom,
        length_all=length_all,
        length_deviation=length_deviation,
        length_method=length_method,
        length_range=length_range,
        discard=discard,
        no_outgroup=no_outgroup,
        group_segments=group_segments,
        outgroup_candidates_taxid_min_genomes=outgroup_candidates_taxid_min_genomes,
        outgroup_treebuilder=outgroup_treebuilder,
        glance=glance,
        print_fasta_headers=print_fasta_headers,
        ignore_duplicates=ignore_duplicates,
        keep_files=keep_files,
    )


def dereplicate_params(
    *,
    tool: Any = _UNSET,
    primary_ani: Any = _UNSET,
    secondary_ani: Any = _UNSET,
    aligned_fraction: Any = _UNSET,
    threads: Any = _UNSET,
    process_size: Any = _UNSET,
    num_processes: Any = _UNSET,
    pre_primary_ani: Any = _UNSET,
    pre_secondary_ani: Any = _UNSET,
    reduce: Any = _UNSET,
    target_reps: Any = _UNSET,
    extra: Any = _UNSET,
    allow_incomplete: Any = _UNSET,
    keeper: Any = _UNSET,
) -> DereplicateParams:
    from ..core.errors import UserInputError
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.dereplicate import DereplicateParams

    if tool is not _UNSET:
        _require_choice(tool, {"auto", *_derep_registry.names()}, "--tool")
    if reduce is not _UNSET:
        _require_choice(reduce, {"none", "species", "genus"}, "--reduce")
    if keeper is not _UNSET:
        _require_choice(keeper, {"quality", "tool"}, "--keeper")
    if target_reps is not _UNSET and target_reps < 0:
        raise UserInputError(f"--target-reps must be >= 0, got {target_reps}.")
    for value, label in (
        (primary_ani, "--primary-ani"),
        (secondary_ani, "--secondary-ani"),
        (aligned_fraction, "--aligned-fraction"),
        (pre_primary_ani, "--pre-primary-ani"),
        (pre_secondary_ani, "--pre-secondary-ani"),
    ):
        if value is not _UNSET:
            _require_unit_interval(value, label)
    return _build(
        DereplicateParams,
        tool=tool,
        primary_ani=primary_ani,
        secondary_ani=secondary_ani,
        aligned_fraction=aligned_fraction,
        threads=threads,
        process_size=process_size,
        num_processes=num_processes,
        pre_primary_ani=pre_primary_ani,
        pre_secondary_ani=pre_secondary_ani,
        reduce=reduce,
        target_reps=target_reps,
        extra=extra,
        allow_incomplete=allow_incomplete,
        keeper=keeper,
    )


def phylo_params(
    *,
    treebuilder: Any = _UNSET,
    msa_source: Any = _UNSET,
    aligner: Any = _UNSET,
    snptyper: Any = _UNSET,
    all_genomes: Any = _UNSET,
    no_outgroup: Any = _UNSET,
    bootstrap: Any = _UNSET,
    reference: Any = _UNSET,
    threads: Any = _UNSET,
    extra: Any = _UNSET,
    allow_incomplete: Any = _UNSET,
) -> PhyloParams:
    from ..aligners.base import registry as _aln_registry
    from ..core.errors import UserInputError
    from ..snptypers.base import registry as _snp_registry
    from ..stages.phylo import PhyloParams
    from ..treebuilders.base import registry as _tb_registry

    if treebuilder is not _UNSET:
        _require_choice(treebuilder, {"auto", *_tb_registry.names()}, "--treebuilder")
    effective_source = msa_source if msa_source is not _UNSET else "aligner"
    _require_choice(effective_source, {"aligner", "snptype"}, "--msa-source")
    if effective_source == "aligner":
        if aligner is not _UNSET:
            _require_choice(aligner, set(_aln_registry.names()), "--aligner")
    elif snptyper is not _UNSET:
        _require_choice(snptyper, set(_snp_registry.names()), "--snptyper")
    # The phylo stage reads extra["mask"] only on the snptype path; a mask
    # requested with the aligner source would be dropped without a trace.
    mask = extra.get("mask") if extra is not _UNSET else None
    if mask is not None:
        require_mask(mask)
        if effective_source != "snptype":
            raise UserInputError("--mask applies only with --msa-source snptype.")
    return _build(
        PhyloParams,
        treebuilder=treebuilder,
        msa_source=msa_source,
        aligner=aligner,
        snptyper=snptyper,
        all_genomes=all_genomes,
        no_outgroup=no_outgroup,
        bootstrap=bootstrap,
        reference=reference,
        threads=threads,
        extra=extra,
        allow_incomplete=allow_incomplete,
    )


def tree2tax_params(
    *,
    node_basename: Any = _UNSET,
    root_name: Any = _UNSET,
    remove_outgroup: Any = _UNSET,
    include_dereplicated: Any = _UNSET,
    collapse_support: Any = _UNSET,
    collapse_length: Any = _UNSET,
) -> Tree2taxParams:
    from ..stages.tree2tax import Tree2taxParams

    return _build(
        Tree2taxParams,
        node_basename=node_basename,
        root_name=root_name,
        remove_outgroup=remove_outgroup,
        include_dereplicated=include_dereplicated,
        collapse_support=collapse_support,
        collapse_length=collapse_length,
    )
