"""Stateless data-channel steps: genome-fetch, dereplicate-chunk, dereplicate-merge,
assemble-run, genome-qc, reads-gather.

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
    HELP_ALIGNED_FRACTION,
    HELP_PRIMARY_ANI,
    HELP_SECONDARY_ANI,
    HELP_THREADS,
    _aligner_help,
    _assembler_help,
    _classifier_help,
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
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani", help=HELP_PRIMARY_ANI),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani", help=HELP_SECONDARY_ANI),
    aligned_fraction: float = typer.Option(
        0.50, "-af", "--aligned-fraction", help=HELP_ALIGNED_FRACTION
    ),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS),
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
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS),
    mask: str = typer.Option(
        "none", "--mask", help="Recombination masking for --msa-source snptype."
    ),
    msa_only: bool = typer.Option(
        False,
        "--msa-only",
        help="Build the alignment and stop, writing msa.fasta (for a separate tree step).",
    ),
    msa: Path | None = typer.Option(
        None, "--msa", help="Build the tree from this alignment instead of constructing one."
    ),
    versions_out: Path | None = typer.Option(
        None, "--versions-out", help="Write resolved tool versions (YAML fragment) here."
    ),
) -> None:
    """Build a phylogeny from a genomes directory (stateless data-channel step)."""
    from ..aligners.base import registry as _aln_registry
    from ..snptypers.base import registry as _snp_registry
    from ..stages.phylo import PhyloBuildParams, PhyloParams, phylo_build
    from ..treebuilders.base import registry as _tb_registry
    from .param_builders import require_mask

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        _require_choice(treebuilder, {"auto", *_tb_registry.names()}, "--treebuilder")
        _require_choice(msa_source, {"aligner", "snptype"}, "--msa-source")
        if msa_source == "aligner":
            _require_choice(aligner, set(_aln_registry.names()), "--aligner")
        else:
            _require_choice(snptyper, set(_snp_registry.names()), "--snptyper")
        require_mask(mask)
        if mask != "none" and msa_source != "snptype":
            raise UserInputError("--mask applies only with --msa-source snptype.")
        if msa_only and msa is not None:
            raise UserInputError("--msa-only builds an alignment; --msa consumes one. Pick one.")
        if msa is not None and not msa.is_file():
            raise UserInputError(f"Alignment not found: {msa}")

        phylo_params = PhyloParams(
            treebuilder=treebuilder,
            msa_source=msa_source,
            aligner=aligner,
            snptyper=snptyper,
            no_outgroup=no_outgroup,
            bootstrap=bootstrap,
            reference=reference,
            threads=threads,
            extra={
                **_parse_key_values(aligner_arg, "--aligner-arg"),
                **({"mask": mask} if mask != "none" else {}),
            },
        )
        phylo_build(
            PhyloBuildParams(
                genomes_dir=genomes_dir,
                out_dir=out_dir,
                outgroup_dir=outgroup_dir,
                outgroup_accession=outgroup_accession,
                phylo=phylo_params,
                versions_out=versions_out,
                msa_only=msa_only,
                msa=msa,
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
    root_name: str = typer.Option("root", "--root-name", help="Name for the root node."),
    remove_outgroup: bool = typer.Option(False, "--remove-outgroup", help="Drop outgroup."),
    include_dereplicated: bool = typer.Option(
        True,
        "--include-dereplicated/--no-include-dereplicated",
        help="List redundant genomes under their representative.",
    ),
    versions_out: Path | None = typer.Option(
        None, "--versions-out", help="Write resolved tool versions (YAML fragment) here."
    ),
    collapse_support: float | None = typer.Option(
        None,
        "--collapse-support",
        min=0.0,
        max=1.0,
        help="Merge nodes whose support is below this fraction into their parent.",
    ),
    collapse_length: float | None = typer.Option(
        None,
        "--collapse-length",
        min=0.0,
        help="Merge nodes whose branch is shorter than this length into their parent.",
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
                collapse_support=collapse_support,
                collapse_length=collapse_length,
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
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani", help=HELP_PRIMARY_ANI),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani", help=HELP_SECONDARY_ANI),
    aligned_fraction: float = typer.Option(
        0.50, "-af", "--aligned-fraction", help=HELP_ALIGNED_FRACTION
    ),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS),
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
    reduce: str = typer.Option(
        "none",
        "--reduce",
        help="Taxonomy-aware reduction after the merge: none, species, or genus "
        "(one representative per taxon; taxonomy from --selection-tsv or the filenames).",
    ),
    target_reps: int = typer.Option(
        0,
        "--target-reps",
        help="Target representative count: search --secondary-ani of the merge pass "
        "to land near it (0 = off; re-runs the merge per search step).",
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
        _require_choice(reduce, {"none", "species", "genus"}, "--reduce")
        if target_reps < 0:
            raise UserInputError("--target-reps must be 0 (off) or a positive count.")
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
                reduce=reduce,
                target_reps=target_reps,
                versions_out=versions_out,
            ),
            logger,
        )


# --- reads chain ------------------------------------------------------------------


@app.command(name="assemble-run")
def assemble_run_cmd(
    reads_tsv: Path = typer.Option(..., "--reads-tsv", help="reads.tsv from the reads stage."),
    run: str = typer.Option(
        ..., "--run", help="The run accession (a row of reads.tsv) to assemble."
    ),
    out_dir: Path = typer.Option(
        ..., "-o", "--out", help="Output dir: contigs.fasta and assembly.ok, or excused_runs.tsv."
    ),
    assembler: str = typer.Option("auto", "--assembler", help=_assembler_help()),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS),
    memory_gb: int = typer.Option(
        16,
        "--memory-gb",
        min=1,
        help="Memory hint for the assembly, in GB, for tools that cap RAM.",
    ),
    min_contig_length: int = typer.Option(
        500, "--min-contig-length", min=0, help="Drop contigs shorter than this many bases."
    ),
    keep_reads: bool = typer.Option(
        False, "--keep-reads", help="Keep the downloaded FASTQ files after assembling."
    ),
    keep_files: bool = typer.Option(
        False, "--keep-files", help="Keep the assembler scratch directory."
    ),
    tool_arg: list[str] = typer.Option(
        [], "--tool-arg", help="Assembler tuning as key=value (repeatable), e.g. mode=nano-raw."
    ),
    versions_out: Path | None = typer.Option(
        None, "--versions-out", help="Write resolved tool versions (YAML fragment) here."
    ),
) -> None:
    """Fetch and assemble one run of a reads.tsv (stateless data-channel step)."""
    from ..assemblers.base import registry as asm_registry
    from ..stages.assemble_steps import AssembleRunParams, assemble_run

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        _require_choice(assembler, {"auto", *asm_registry.names()}, "--assembler")
        assemble_run(
            AssembleRunParams(
                reads_tsv=reads_tsv,
                run=run,
                out_dir=out_dir,
                assembler=assembler,
                threads=threads,
                memory_gb=memory_gb,
                min_contig_length=min_contig_length,
                keep_reads=keep_reads,
                keep_files=keep_files,
                extra=_parse_key_values(tool_arg, "--tool-arg"),
                versions_out=versions_out,
            ),
            logger,
        )


@app.command(name="genome-qc")
def genome_qc_cmd(
    assemblies: Path = typer.Option(
        ..., "--assemblies", help="Directory of assemble-run output dirs, one per run."
    ),
    out_dir: Path = typer.Option(
        ..., "-o", "--out", help="Output dir for quality.tsv and classification.tsv."
    ),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS),
    checkm2_db: Path | None = typer.Option(
        None,
        "--checkm2-db",
        help="CheckM2 DIAMOND database; enables quality scoring (or set CHECKM2DB).",
    ),
    classifier: str = typer.Option("auto", "--classifier", help=_classifier_help()),
    gtdb_sketch: Path | None = typer.Option(
        None,
        "--gtdb-sketch",
        help="GTDB sourmash sketch database (.sig.zip); enables classification "
        "(or set REPGENR_GTDB_SKETCH).",
    ),
    gtdb_lineages: Path | None = typer.Option(
        None,
        "--gtdb-lineages",
        help="The lineages CSV published with the sketch (or set REPGENR_GTDB_LINEAGES).",
    ),
    tool_arg: list[str] = typer.Option(
        [], "--tool-arg", help="Classifier tuning as key=value (repeatable)."
    ),
    versions_out: Path | None = typer.Option(
        None, "--versions-out", help="Write resolved tool versions (YAML fragment) here."
    ),
) -> None:
    """Score (CheckM2) and classify a batch of assemblies (stateless data-channel step)."""
    from ..classifiers.base import registry as cls_registry
    from ..stages.assemble_steps import GenomeQcParams, genome_qc

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        _require_choice(classifier, {"auto", "none", *cls_registry.names()}, "--classifier")
        genome_qc(
            GenomeQcParams(
                assemblies_dir=assemblies,
                out_dir=out_dir,
                threads=threads,
                checkm2_db=None if checkm2_db is None else str(checkm2_db),
                classifier=classifier,
                gtdb_sketch=None if gtdb_sketch is None else str(gtdb_sketch),
                gtdb_lineages=None if gtdb_lineages is None else str(gtdb_lineages),
                extra=_parse_key_values(tool_arg, "--tool-arg"),
                versions_out=versions_out,
            ),
            logger,
        )


@app.command(name="reads-gather")
def reads_gather_cmd(
    reads_tsv: Path = typer.Option(..., "--reads-tsv", help="reads.tsv from the reads stage."),
    assemblies: Path = typer.Option(
        ..., "--assemblies", help="Directory of assemble-run output dirs, one per run."
    ),
    out_dir: Path = typer.Option(
        ..., "-o", "--out", help="Output dir for genomes/, selection.tsv and the stats tables."
    ),
    qc: Path | None = typer.Option(
        None, "--qc", help="genome-qc output dir (quality.tsv, classification.tsv), if it ran."
    ),
    min_completeness: float = typer.Option(
        50.0, "--min-completeness", min=0.0, max=100.0, help="CheckM2 completeness floor."
    ),
    max_contamination: float = typer.Option(
        10.0, "--max-contamination", min=0.0, max=100.0, help="CheckM2 contamination ceiling."
    ),
) -> None:
    """Write the genome contract from per-run assemblies (stateless data-channel step)."""
    from ..stages.assemble_steps import ReadsGatherParams, reads_gather

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        reads_gather(
            ReadsGatherParams(
                reads_tsv=reads_tsv,
                assemblies_dir=assemblies,
                out_dir=out_dir,
                qc_dir=qc,
                min_completeness=min_completeness,
                max_contamination=max_contamination,
            ),
            logger,
        )
