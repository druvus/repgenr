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

from ..core.logging import configure_logging
from .base import (
    _RUN_STATE,
    DEFAULT_THREADS,
    PIPELINE_BACTERIAL,
    PIPELINE_VIRAL,
    _require_choice,
    _require_unit_interval,
    _run,
    app,
    stage_errors,
)


@app.command()
def run(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory (created)."),
    viral: bool = typer.Option(
        False, "--viral", help="Run the viral chain (vmetadata -> vgenome) instead of bacterial."
    ),
    # --- selection: bacterial (GTDB) ---
    dataset: str = typer.Option("rep", "-d", "--dataset", help="all or rep (bacterial)."),
    level: str | None = typer.Option(None, "-l", "--level", help="family/genus/species."),
    target_family: str | None = typer.Option(None, "-tf", "--target-family"),
    target_genus: str | None = typer.Option(None, "-tg", "--target-genus"),
    target_species: str | None = typer.Option(None, "-ts", "--target-species"),
    release: str | None = typer.Option(None, "-r", "--release", help="GTDB release (tsv source)."),
    gtdb_version: str | None = typer.Option(None, "--gtdb-version", help="bac120/ar53."),
    metadata_source: str = typer.Option("tsv", "--metadata-source", help="tsv or api."),
    outgroup_accession: str | None = typer.Option(None, "--outgroup-accession"),
    # --- selection: viral (NCBI Virus) ---
    target: str | None = typer.Option(None, "-t", "--target", help="Virus taxon (viral)."),
    viral_source: str = typer.Option("ncbi_virus", "--viral-source", help="ncbi_virus or bvbrc."),
    group_segments: bool = typer.Option(False, "--group-segments", help="Group viral segments."),
    # --- dereplication ---
    derep_tool: str = typer.Option("skder", "--tool", help="auto/skder/drep/galah/sourmash."),
    primary_ani: float = typer.Option(0.90, "--primary-ani"),
    secondary_ani: float = typer.Option(0.99, "--secondary-ani"),
    aligned_fraction: float = typer.Option(0.50, "--aligned-fraction"),
    # --- phylogeny ---
    treebuilder: str = typer.Option("iqtree", "--treebuilder"),
    aligner: str = typer.Option("progressivemauve", "--aligner"),
    no_outgroup: bool = typer.Option(False, "--no-outgroup"),
    # --- common ---
    threads: int = typer.Option(DEFAULT_THREADS, "--threads"),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Print the stages and key parameters, then exit."
    ),
) -> None:
    """Run the whole pipeline end to end (bacterial by default, --viral for viruses)."""
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.dereplicate import DereplicateParams
    from ..stages.phylo import PhyloParams
    from ..stages.tree2tax import Tree2taxParams

    logger = configure_logging(
        workdir if workdir.exists() else None, level=_RUN_STATE["log_level"]
    )
    with stage_errors(logger):
        _require_choice(derep_tool, {"auto", *_derep_registry.names()}, "--tool")
        _require_unit_interval(primary_ani, "--primary-ani")
        _require_unit_interval(secondary_ani, "--secondary-ani")
        _require_unit_interval(aligned_fraction, "--aligned-fraction")

    if dry_run:
        chain = PIPELINE_VIRAL if viral else PIPELINE_BACTERIAL
        typer.echo(
            f"[dry-run] {'viral' if viral else 'bacterial'} pipeline in {workdir}:"
        )
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
            f"secondary_ani={secondary_ani}; phylo: treebuilder={treebuilder}, "
            f"aligner={aligner}; threads={threads}"
        )
        typer.echo("[dry-run] no work done.")
        return

    if viral:
        from ..stages.vgenome import VgenomeParams
        from ..stages.vmetadata import VmetadataParams

        _run("vmetadata", workdir, lambda: VmetadataParams(
            target=target, source=viral_source,
        ), create=True)
        _run("vgenome", workdir, lambda: VgenomeParams(
            target_genus=target_genus, target_species=target_species,
            no_outgroup=no_outgroup, group_segments=group_segments,
        ))
    else:
        from ..stages.genome import GenomeParams
        from ..stages.metadata import MetadataParams

        _run("metadata", workdir, lambda: MetadataParams(
            dataset=dataset, level=level or "", source=metadata_source,
            release=release, version=gtdb_version,
            target_family=target_family, target_genus=target_genus,
            target_species=target_species, outgroup_accession=outgroup_accession,
        ), create=True)
        _run("genome", workdir, lambda: GenomeParams())

    _run("dereplicate", workdir, lambda: DereplicateParams(
        tool=derep_tool, primary_ani=primary_ani, secondary_ani=secondary_ani,
        aligned_fraction=aligned_fraction, threads=threads,
        extra={"virus": True} if viral else {},
    ))
    _run("phylo", workdir, lambda: PhyloParams(
        treebuilder=treebuilder, aligner=aligner, no_outgroup=no_outgroup, threads=threads,
    ))
    _run("tree2tax", workdir, lambda: Tree2taxParams(include_dereplicated=True))

    typer.echo(
        f"\nPipeline complete. Deliverables in {workdir}: "
        "tree2tax.tsv, genomes_map.tsv, tree/tree.nwk"
    )
