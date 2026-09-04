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
    from ..stages.metadata import MetadataParams
    from ..stages.phylo import PhyloParams
    from ..stages.tree2tax import Tree2taxParams
    from ..stages.vgenome import VgenomeParams
    from ..stages.vmetadata import VmetadataParams

# Sentinel for "caller did not set this": the argument is not passed to the
# dataclass, whose own default then applies.
_UNSET: Any = object()


def _build(cls: type, **kwargs: Any):
    return cls(**{k: v for k, v in kwargs.items() if v is not _UNSET})


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
        _require_choice(source, {"ncbi_virus", "bvbrc"}, "--source")
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


def vgenome_params(**kwargs: Any) -> VgenomeParams:
    from ..stages.vgenome import VgenomeParams

    return _build(VgenomeParams, **kwargs)


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
    all_genomes: Any = _UNSET,
    include_dereplicated: Any = _UNSET,
) -> Tree2taxParams:
    from ..stages.tree2tax import Tree2taxParams

    return _build(
        Tree2taxParams,
        node_basename=node_basename,
        root_name=root_name,
        remove_outgroup=remove_outgroup,
        all_genomes=all_genomes,
        include_dereplicated=include_dereplicated,
    )
