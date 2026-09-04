"""Stateless data-channel steps: genome-fetch, dereplicate-chunk, dereplicate-merge.

These run as discrete Nextflow process steps (no shared workdir); they read
explicit inputs (a selection.tsv or a file-of-filenames) and write a result
directory, rather than going through the workdir-bound :func:`_run` harness.
"""

from __future__ import annotations

from pathlib import Path

import typer

from ..core.errors import UserInputError
from ..core.logging import configure_logging
from .base import (
    _RUN_STATE,
    DEFAULT_THREADS,
    _aligner_help,
    _derep_help,
    _parse_key_values,
    _read_path_fofn,
    _require_choice,
    _require_unit_interval,
    _snp_help,
    _tree_help,
    app,
    gated_extra,
    stage_errors,
)


@app.command(name="genome-fetch")
def genome_fetch_cmd(
    selection: Path = typer.Option(
        ..., "--selection", help="selection.tsv from the metadata stage."
    ),
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output dir for downloaded genomes."),
    keep_files: bool = typer.Option(False, "--keep-files", help="Keep download intermediates."),
    versions_out: Path | None = typer.Option(
        None, "--versions-out", help="Write resolved tool versions (YAML fragment) here."
    ),
) -> None:
    """Download genomes listed in a selection.tsv (stateless data-channel step)."""
    from ..stages.genome_steps import GenomeFetchParams, genome_fetch

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        genome_fetch(
            GenomeFetchParams(
                selection_tsv=selection,
                out_dir=out_dir,
                keep_files=keep_files,
                versions_out=versions_out,
            ),
            logger,
        )


@app.command(name="dereplicate-chunk")
def dereplicate_chunk_cmd(
    genomes_fofn: Path = typer.Option(
        ..., "--genomes-fofn", help="File of genome FASTA paths, one per line."
    ),
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output directory for the chunk result."),
    tool: str = typer.Option("skder", "--tool", help=_derep_help(auto=False)),
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani"),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani"),
    aligned_fraction: float = typer.Option(0.50, "-af", "--aligned-fraction"),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", min=1),
    virus: bool = typer.Option(False, "--virus", help="Pass virus-tuned parameters to the tool."),
    tool_arg: list[str] = typer.Option(
        [], "--tool-arg", help="Tool tuning as key=value (repeatable)."
    ),
    selection_tsv: Path | None = typer.Option(
        None,
        "--selection-tsv",
        help="selection.tsv with quality columns; enables quality-aware representatives.",
    ),
    keeper: str = typer.Option(
        "quality",
        "--keeper",
        help="Representative choice when --selection-tsv is given: "
        "quality (manifest completeness/contamination) or tool (adapter's own pick).",
    ),
    versions_out: Path | None = typer.Option(
        None, "--versions-out", help="Write resolved tool versions (YAML fragment) here."
    ),
) -> None:
    """Dereplicate one chunk of genomes (scatter step; writes a chunk result dir)."""
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.derep_steps import ChunkParams, dereplicate_chunk

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        # A concrete tool (not 'auto') so every scattered chunk and the merge agree.
        _require_choice(tool, set(_derep_registry.names()), "--tool")
        _require_unit_interval(primary_ani, "--primary-ani")
        _require_unit_interval(secondary_ani, "--secondary-ani")
        _require_unit_interval(aligned_fraction, "--aligned-fraction")
        _require_choice(keeper, {"quality", "tool"}, "--keeper")
        genomes = _read_path_fofn(genomes_fofn)
        dereplicate_chunk(
            ChunkParams(
                tool=tool,
                genomes=genomes,
                out_dir=out_dir,
                primary_ani=primary_ani,
                secondary_ani=secondary_ani,
                aligned_fraction=aligned_fraction,
                threads=threads,
                extra={
                    **_parse_key_values(tool_arg, "--tool-arg"),
                    **(gated_extra(_derep_registry, tool, "virus", True) if virus else {}),
                },
                selection_tsv=selection_tsv,
                keeper=keeper,
                versions_out=versions_out,
            ),
            logger,
        )


@app.command(name="phylo-build")
def phylo_build_cmd(
    genomes_dir: Path = typer.Option(
        ..., "--genomes-dir", help="Directory of genome FASTA files to build the tree from."
    ),
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output dir (writes tree/tree.nwk)."),
    outgroup_dir: Path | None = typer.Option(
        None, "--outgroup-dir", help="Directory holding the outgroup genome file(s)."
    ),
    outgroup_accession: Path | None = typer.Option(
        None, "--outgroup-accession", help="File naming the outgroup accession."
    ),
    treebuilder: str = typer.Option("iqtree", "--treebuilder", help=_tree_help()),
    msa_source: str = typer.Option("aligner", "--msa-source", help="aligner or snptype."),
    aligner: str = typer.Option("progressivemauve", "--aligner", help=_aligner_help()),
    snptyper: str = typer.Option("simple", "--snptyper", help=_snp_help()),
    no_outgroup: bool = typer.Option(False, "--no-outgroup", help="Do not root with an outgroup."),
    bootstrap: int = typer.Option(
        0, "-B", "--bootstrap", min=0, help="Bootstrap replicates (>=1000)."
    ),
    reference: str | None = typer.Option(None, "--reference", help="Reference genome filename."),
    aligner_arg: list[str] = typer.Option(
        [], "--aligner-arg", help="Aligner tuning as key=value (repeatable)."
    ),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", min=1),
    versions_out: Path | None = typer.Option(
        None, "--versions-out", help="Write resolved tool versions (YAML fragment) here."
    ),
) -> None:
    """Build a phylogeny from a genomes directory (stateless data-channel step)."""
    from ..aligners.base import registry as _aln_registry
    from ..snptypers.base import registry as _snp_registry
    from ..stages.phylo import PhyloBuildParams, PhyloParams, phylo_build
    from ..treebuilders.base import registry as _tb_registry

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        _require_choice(treebuilder, {"auto", *_tb_registry.names()}, "--treebuilder")
        _require_choice(msa_source, {"aligner", "snptype"}, "--msa-source")
        if msa_source == "aligner":
            _require_choice(aligner, set(_aln_registry.names()), "--aligner")
        else:
            _require_choice(snptyper, set(_snp_registry.names()), "--snptyper")

        phylo_params = PhyloParams(
            treebuilder=treebuilder,
            msa_source=msa_source,
            aligner=aligner,
            snptyper=snptyper,
            no_outgroup=no_outgroup,
            bootstrap=bootstrap,
            reference=reference,
            threads=threads,
            extra=_parse_key_values(aligner_arg, "--aligner-arg"),
        )
        phylo_build(
            PhyloBuildParams(
                genomes_dir=genomes_dir,
                out_dir=out_dir,
                outgroup_dir=outgroup_dir,
                outgroup_accession=outgroup_accession,
                phylo=phylo_params,
                versions_out=versions_out,
            ),
            logger,
        )


@app.command(name="tree2tax-relations")
def tree2tax_relations_cmd(
    tree: Path = typer.Option(..., "--tree", help="Rooted/unrooted tree in Newick (tree.nwk)."),
    out_dir: Path = typer.Option(
        ..., "-o", "--out", help="Output dir (writes tree2tax.tsv + genomes_map.tsv)."
    ),
    clusters: Path | None = typer.Option(
        None, "--clusters", help="derep clusters.tsv (for --include-dereplicated)."
    ),
    outgroup_dir: Path | None = typer.Option(
        None, "--outgroup-dir", help="Directory holding the outgroup genome file(s)."
    ),
    outgroup_accession: Path | None = typer.Option(
        None, "--outgroup-accession", help="File naming the outgroup accession."
    ),
    node_basename: str | None = typer.Option(None, "--node-basename", help="Prefix for nodes."),
    root_name: str = typer.Option("root", "-r", "--root-name", help="Name for the root node."),
    remove_outgroup: bool = typer.Option(False, "--remove-outgroup", help="Drop outgroup."),
    include_dereplicated: bool = typer.Option(
        False, "--include-dereplicated", help="List redundant genomes under their representative."
    ),
    versions_out: Path | None = typer.Option(
        None, "--versions-out", help="Write resolved tool versions (YAML fragment) here."
    ),
) -> None:
    """Emit FlexTaxD relations from a tree (stateless data-channel step)."""
    from ..stages.tree2tax import Tree2taxStepParams, tree2tax_relations

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        tree2tax_relations(
            Tree2taxStepParams(
                tree=tree,
                out_dir=out_dir,
                clusters=clusters,
                outgroup_dir=outgroup_dir,
                outgroup_accession=outgroup_accession,
                node_basename=node_basename,
                root_name=root_name,
                remove_outgroup=remove_outgroup,
                include_dereplicated=include_dereplicated,
                versions_out=versions_out,
            ),
            logger,
        )


@app.command(name="dereplicate-merge")
def dereplicate_merge_cmd(
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output dir for the merged result."),
    chunk_dir: list[Path] = typer.Option(
        [], "--chunk-dir", help="A chunk result directory (repeatable)."
    ),
    chunk_fofn: Path | None = typer.Option(
        None, "--chunk-fofn", help="File listing chunk result directories, one per line."
    ),
    tool: str = typer.Option("skder", "--tool", help=_derep_help(auto=False)),
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani"),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani"),
    aligned_fraction: float = typer.Option(0.50, "-af", "--aligned-fraction"),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", min=1),
    virus: bool = typer.Option(False, "--virus", help="Pass virus-tuned parameters to the tool."),
    tool_arg: list[str] = typer.Option(
        [], "--tool-arg", help="Tool tuning as key=value (repeatable)."
    ),
    selection_tsv: Path | None = typer.Option(
        None,
        "--selection-tsv",
        help="selection.tsv with quality columns; enables quality-aware representatives.",
    ),
    keeper: str = typer.Option(
        "quality",
        "--keeper",
        help="Representative choice when --selection-tsv is given: "
        "quality (manifest completeness/contamination) or tool (adapter's own pick).",
    ),
    versions_out: Path | None = typer.Option(
        None, "--versions-out", help="Write resolved tool versions (YAML fragment) here."
    ),
) -> None:
    """Dereplicate the union of chunk representatives (gather step)."""
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.derep_steps import MergeParams, dereplicate_merge

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        _require_choice(tool, set(_derep_registry.names()), "--tool")
        _require_unit_interval(primary_ani, "--primary-ani")
        _require_unit_interval(secondary_ani, "--secondary-ani")
        _require_unit_interval(aligned_fraction, "--aligned-fraction")
        _require_choice(keeper, {"quality", "tool"}, "--keeper")
        chunk_dirs = list(chunk_dir)
        if chunk_fofn is not None:
            chunk_dirs += _read_path_fofn(chunk_fofn)
        if not chunk_dirs:
            raise UserInputError("Provide at least one --chunk-dir or a --chunk-fofn.")
        dereplicate_merge(
            MergeParams(
                tool=tool,
                chunk_dirs=chunk_dirs,
                out_dir=out_dir,
                primary_ani=primary_ani,
                secondary_ani=secondary_ani,
                aligned_fraction=aligned_fraction,
                threads=threads,
                extra={
                    **_parse_key_values(tool_arg, "--tool-arg"),
                    **(gated_extra(_derep_registry, tool, "virus", True) if virus else {}),
                },
                selection_tsv=selection_tsv,
                keeper=keeper,
                versions_out=versions_out,
            ),
            logger,
        )
