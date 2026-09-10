"""One-shot pipeline orchestrator: ``repgenr run``.

Chains the canonical stages so a user need not invoke five commands by hand. It
forwards the common options and relies on stage defaults for the rest; for full
per-stage control use the individual commands. Each stage goes through the same
:func:`_run` harness, so the resume guard applies -- re-running ``run`` skips
stages already completed with the same parameters.
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
    HELP_NO_OUTGROUP,
    HELP_OUTGROUP_ACCESSION,
    HELP_PRIMARY_ANI,
    HELP_SECONDARY_ANI,
    HELP_TARGET_FAMILY,
    HELP_TARGET_GENUS,
    HELP_TARGET_SPECIES,
    HELP_THREADS,
    PIPELINE_BACTERIAL,
    PIPELINE_VIRAL,
    _aligner_help,
    _derep_help,
    _parse_key_values,
    _require_choice,
    _run,
    _snp_help,
    _tree_help,
    app,
    gated_extra,
    stage_errors,
)
from .cmd_viral import _validate_released_after


def _virus_extra(derep_tool: str, viral: bool) -> dict:
    """Virus-tuned extras, injected only when the tool actually reads them.

    Passing extra["virus"] to a tool that ignores it would silently change the
    resume fingerprint for nothing (and imply tuning that never happened).
    """
    if not viral:
        return {}
    from ..dereplicators.base import registry as _derep_registry

    return gated_extra(_derep_registry, derep_tool, "virus", True)


def _preflight_tools(
    derep_tool: str,
    treebuilder: str,
    msa_source: str,
    aligner: str,
    snptyper: str,
) -> None:
    """Check every external tool the chain will need before the first stage runs.

    Each stage preflights its own adapter, but by then the earlier stages have
    already downloaded and dereplicated; a missing tree builder should fail in
    the first second, not after the genome download. ``auto`` choices are
    resolved from the genome count inside their stage and are skipped here.
    """
    from ..dereplicators.base import registry as derep_registry
    from ..treebuilders.base import InputKind
    from ..treebuilders.base import registry as tb_registry

    if derep_tool != "auto":
        derep_registry.create(derep_tool).preflight()
    if treebuilder == "auto":
        return
    builder = tb_registry.create(treebuilder)
    builder.preflight()
    if builder.input_kind != InputKind.MSA_FASTA:
        return
    if msa_source == "aligner":
        from ..aligners.base import registry as aligner_registry

        aligner_registry.create(aligner).preflight()
    elif msa_source == "snptype":
        from ..snptypers.base import registry as snp_registry

        snp_registry.create(snptyper).preflight()


def _msa_source_summary(treebuilder: str, msa_source: str, aligner: str, snptyper: str) -> str:
    """The MSA source a dry run will use, or nothing for an alignment-free builder."""
    from ..treebuilders.base import InputKind
    from ..treebuilders.base import registry as tb_registry

    if treebuilder != "auto" and tb_registry.create(treebuilder).input_kind == InputKind.GENOMES:
        return " (alignment-free)"
    if msa_source == "snptype":
        return f", snptyper={snptyper}"
    return f", aligner={aligner}"


@app.command()
def run(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory (created)."),
    viral: bool = typer.Option(
        False, "--viral", help="Run the viral chain (vmetadata -> vgenome) instead of bacterial."
    ),
    # --- selection: bacterial (GTDB) ---
    dataset: str = typer.Option("rep", "-d", "--dataset", help="all or rep (bacterial)."),
    level: str | None = typer.Option(None, "-l", "--level", help="family/genus/species."),
    target_family: str | None = typer.Option(
        None, "-tf", "--target-family", help=HELP_TARGET_FAMILY
    ),
    target_genus: str | None = typer.Option(None, "-tg", "--target-genus", help=HELP_TARGET_GENUS),
    target_species: str | None = typer.Option(
        None, "-ts", "--target-species", help=HELP_TARGET_SPECIES
    ),
    release: str | None = typer.Option(None, "-r", "--release", help="GTDB release (tsv source)."),
    gtdb_version: str | None = typer.Option(None, "--gtdb-version", help="bac120/ar53."),
    metadata_source: str = typer.Option("tsv", "--metadata-source", help="tsv or api."),
    outgroup_accession: str | None = typer.Option(
        None, "--outgroup-accession", help=HELP_OUTGROUP_ACCESSION
    ),
    limit: int | None = typer.Option(
        None,
        "--limit",
        min=1,
        help="Keep at most N genomes, round-robin over species by CheckM quality (bacterial).",
    ),
    # --- selection: viral (NCBI Virus) ---
    target: str | None = typer.Option(None, "--target", help="Virus taxon (viral)."),
    viral_source: str = typer.Option("ncbi_virus", "--viral-source", help="ncbi_virus or bvbrc."),
    complete_only: bool = typer.Option(
        False, "--complete-only", help="ncbi_virus: only COMPLETE sequences (viral)."
    ),
    host: str | None = typer.Option(
        None, "--host", help="ncbi_virus: restrict to a host species (viral)."
    ),
    released_after: str | None = typer.Option(
        None,
        "--released-after",
        callback=_validate_released_after,
        help="ncbi_virus: MM/DD/YYYY (viral).",
    ),
    group_segments: bool = typer.Option(False, "--group-segments", help="Group viral segments."),
    # --- dereplication ---
    derep_tool: str = typer.Option("skder", "--tool", help=_derep_help()),
    primary_ani: float = typer.Option(0.90, "--primary-ani", help=HELP_PRIMARY_ANI),
    secondary_ani: float = typer.Option(0.99, "--secondary-ani", help=HELP_SECONDARY_ANI),
    aligned_fraction: float = typer.Option(0.50, "--aligned-fraction", help=HELP_ALIGNED_FRACTION),
    keeper: str = typer.Option(
        "quality",
        "--keeper",
        help="Representative choice per cluster: quality (CheckM score from GTDB) "
        "or tool (adapter's own).",
    ),
    process_size: int | None = typer.Option(
        None, "-s", "--process-size", help="Chunk size for two-stage dereplication."
    ),
    num_processes: int = typer.Option(
        0, "-p", "--num-processes", help="Parallel chunk workers (0 = auto)."
    ),
    reduce: str = typer.Option(
        "none", "--reduce", help="Taxonomy-aware reduction after ANI: none, species or genus."
    ),
    target_reps: int = typer.Option(
        0, "--target-reps", help="Target representative count (0 = off)."
    ),
    tool_arg: list[str] = typer.Option(
        [], "--tool-arg", help="Dereplicator tuning as key=value (repeatable)."
    ),
    # --- phylogeny ---
    treebuilder: str = typer.Option("iqtree", "--treebuilder", help=_tree_help()),
    msa_source: str = typer.Option("aligner", "--msa-source", help="aligner or snptype."),
    aligner: str = typer.Option("progressivemauve", "--aligner", help=_aligner_help()),
    snptyper: str = typer.Option("simple", "--snptyper", help=_snp_help()),
    no_outgroup: bool = typer.Option(False, "--no-outgroup", help=HELP_NO_OUTGROUP),
    all_genomes: bool = typer.Option(
        False, "--all-genomes", help="Build the tree from all genomes, not the representatives."
    ),
    bootstrap: int = typer.Option(
        0, "-B", "--bootstrap", min=0, help="Bootstrap replicates (>=1000 for IQ-TREE)."
    ),
    reference: str | None = typer.Option(None, "--reference", help="Reference genome filename."),
    aligner_arg: list[str] = typer.Option(
        [], "--aligner-arg", help="Aligner tuning as key=value (repeatable)."
    ),
    mask: str = typer.Option(
        "none", "--mask", help="Recombination masking for --msa-source snptype."
    ),
    # --- taxonomy output ---
    include_dereplicated: bool = typer.Option(
        True,
        "--include-dereplicated/--no-include-dereplicated",
        help="List redundant genomes under their representative in tree2tax.",
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
    # --- common ---
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS),
    allow_incomplete: bool = typer.Option(
        False,
        "--allow-incomplete",
        help="Proceed with a warning when genomes/ is missing selected genomes.",
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Print the stages and key parameters, then exit."
    ),
) -> None:
    """Run the whole pipeline end to end (bacterial by default, --viral for viruses)."""
    from .param_builders import (
        METADATA_DATASETS,
        METADATA_LEVELS,
        METADATA_SOURCES,
        VIRAL_SOURCES,
        dereplicate_params,
        genome_params,
        metadata_params,
        phylo_params,
        tree2tax_params,
        vgenome_params,
        vmetadata_params,
    )

    logger = configure_logging(workdir if workdir.exists() else None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        # Fail fast before any stage runs; the same validation re-runs inside
        # each stage's builder. The selection sources are checked here under
        # the names `run` gives them (the builders know them as --source).
        _require_choice(metadata_source, METADATA_SOURCES, "--metadata-source")
        _require_choice(viral_source, VIRAL_SOURCES, "--viral-source")
        _require_choice(dataset, METADATA_DATASETS, "--dataset")
        if level is not None:
            _require_choice(level, METADATA_LEVELS, "--level")
        derep_extra = {
            **_parse_key_values(tool_arg, "--tool-arg"),
            **_virus_extra(derep_tool, viral),
        }
        phylo_extra = {
            **_parse_key_values(aligner_arg, "--aligner-arg"),
            **({"mask": mask} if mask != "none" else {}),
        }
        dereplicate_params(
            tool=derep_tool,
            primary_ani=primary_ani,
            secondary_ani=secondary_ani,
            aligned_fraction=aligned_fraction,
            keeper=keeper,
            reduce=reduce,
            target_reps=target_reps,
            extra=derep_extra,
        )
        phylo_params(
            treebuilder=treebuilder,
            msa_source=msa_source,
            aligner=aligner,
            snptyper=snptyper,
            extra=phylo_extra,
        )
        if not viral and not level:
            raise UserInputError("The bacterial chain needs -l/--level (family/genus/species).")

    if dry_run:
        chain = PIPELINE_VIRAL if viral else PIPELINE_BACTERIAL
        typer.echo(f"[dry-run] {'viral' if viral else 'bacterial'} pipeline in {workdir}:")
        for stage in chain:
            typer.echo(f"  - {stage}")
        selection = (
            f"target={target}, genus={target_genus}, species={target_species}"
            if viral
            else f"dataset={dataset}, level={level}, "
            f"family={target_family}, genus={target_genus}, species={target_species}"
        )
        typer.echo(f"selection: {selection}")
        typer.echo(
            f"dereplicate: tool={derep_tool}, primary_ani={primary_ani}, "
            f"secondary_ani={secondary_ani}; phylo: treebuilder={treebuilder}"
            f"{_msa_source_summary(treebuilder, msa_source, aligner, snptyper)}; "
            f"threads={threads}"
        )
        typer.echo("[dry-run] no work done.")
        return

    with stage_errors(logger):
        _preflight_tools(derep_tool, treebuilder, msa_source, aligner, snptyper)

    if viral:
        _run(
            "vmetadata",
            workdir,
            lambda: vmetadata_params(
                target=target,
                source=viral_source,
                complete_only=complete_only,
                host=host,
                released_after=released_after,
            ),
            create=True,
        )
        _run(
            "vgenome",
            workdir,
            lambda: vgenome_params(
                target_genus=target_genus,
                target_species=target_species,
                no_outgroup=no_outgroup,
                group_segments=group_segments,
            ),
        )
    else:
        _run(
            "metadata",
            workdir,
            lambda: metadata_params(
                dataset=dataset,
                level=level or "",
                source=metadata_source,
                release=release,
                version=gtdb_version,
                target_family=target_family,
                target_genus=target_genus,
                target_species=target_species,
                outgroup_accession=outgroup_accession,
                limit=limit,
            ),
            create=True,
        )
        _run("genome", workdir, lambda: genome_params())

    _run(
        "dereplicate",
        workdir,
        lambda: dereplicate_params(
            tool=derep_tool,
            primary_ani=primary_ani,
            secondary_ani=secondary_ani,
            aligned_fraction=aligned_fraction,
            threads=threads,
            process_size=process_size,
            num_processes=num_processes,
            reduce=reduce,
            target_reps=target_reps,
            extra=derep_extra,
            keeper=keeper,
            allow_incomplete=allow_incomplete,
        ),
    )
    _run(
        "phylo",
        workdir,
        lambda: phylo_params(
            treebuilder=treebuilder,
            msa_source=msa_source,
            aligner=aligner,
            snptyper=snptyper,
            no_outgroup=no_outgroup,
            all_genomes=all_genomes,
            bootstrap=bootstrap,
            reference=reference,
            threads=threads,
            extra=phylo_extra,
            allow_incomplete=allow_incomplete,
        ),
    )
    _run(
        "tree2tax",
        workdir,
        lambda: tree2tax_params(
            include_dereplicated=include_dereplicated,
            collapse_support=collapse_support,
            collapse_length=collapse_length,
        ),
    )

    typer.echo(
        f"\nPipeline complete. Deliverables in {workdir}: "
        "tree2tax.tsv, genomes_map.tsv, tree/tree.nwk"
    )
