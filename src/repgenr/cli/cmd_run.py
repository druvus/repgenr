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
    PIPELINE_BACTERIAL,
    PIPELINE_VIRAL,
    _aligner_help,
    _derep_help,
    _run,
    _tree_help,
    app,
    gated_extra,
    stage_errors,
)


def _virus_extra(derep_tool: str, viral: bool) -> dict:
    """Virus-tuned extras, injected only when the tool actually reads them.

    Passing extra["virus"] to a tool that ignores it would silently change the
    resume fingerprint for nothing (and imply tuning that never happened).
    """
    if not viral:
        return {}
    from ..dereplicators.base import registry as _derep_registry

    return gated_extra(_derep_registry, derep_tool, "virus", True)


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
    derep_tool: str = typer.Option("skder", "--tool", help=_derep_help()),
    primary_ani: float = typer.Option(0.90, "--primary-ani"),
    secondary_ani: float = typer.Option(0.99, "--secondary-ani"),
    aligned_fraction: float = typer.Option(0.50, "--aligned-fraction"),
    # --- phylogeny ---
    treebuilder: str = typer.Option("iqtree", "--treebuilder", help=_tree_help()),
    msa_source: str = typer.Option("aligner", "--msa-source", help="aligner or snptype."),
    aligner: str = typer.Option("progressivemauve", "--aligner", help=_aligner_help()),
    snptyper: str = typer.Option("simple", "--snptyper", help="SNP typer for snptype source."),
    no_outgroup: bool = typer.Option(False, "--no-outgroup"),
    # --- taxonomy output ---
    include_dereplicated: bool = typer.Option(
        True, "--include-dereplicated/--no-include-dereplicated",
        help="List redundant genomes under their representative in tree2tax.",
    ),
    # --- common ---
    threads: int = typer.Option(DEFAULT_THREADS, "--threads", min=1),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Print the stages and key parameters, then exit."
    ),
) -> None:
    """Run the whole pipeline end to end (bacterial by default, --viral for viruses)."""
    from .param_builders import (
        dereplicate_params,
        genome_params,
        metadata_params,
        phylo_params,
        tree2tax_params,
        vgenome_params,
        vmetadata_params,
    )

    logger = configure_logging(
        workdir if workdir.exists() else None, level=_RUN_STATE["log_level"]
    )
    with stage_errors(logger):
        # Fail fast before any stage runs; the same validation re-runs inside
        # each stage's builder.
        dereplicate_params(
            tool=derep_tool, primary_ani=primary_ani, secondary_ani=secondary_ani,
            aligned_fraction=aligned_fraction,
        )
        phylo_params(
            treebuilder=treebuilder, msa_source=msa_source,
            aligner=aligner, snptyper=snptyper,
        )
        if not viral and not level:
            raise UserInputError("The bacterial chain needs -l/--level (family/genus/species).")

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
        _run("vmetadata", workdir, lambda: vmetadata_params(
            target=target, source=viral_source,
        ), create=True)
        _run("vgenome", workdir, lambda: vgenome_params(
            target_genus=target_genus, target_species=target_species,
            no_outgroup=no_outgroup, group_segments=group_segments,
        ))
    else:
        _run("metadata", workdir, lambda: metadata_params(
            dataset=dataset, level=level or "", source=metadata_source,
            release=release, version=gtdb_version,
            target_family=target_family, target_genus=target_genus,
            target_species=target_species, outgroup_accession=outgroup_accession,
        ), create=True)
        _run("genome", workdir, lambda: genome_params())

    _run("dereplicate", workdir, lambda: dereplicate_params(
        tool=derep_tool, primary_ani=primary_ani, secondary_ani=secondary_ani,
        aligned_fraction=aligned_fraction, threads=threads,
        extra=_virus_extra(derep_tool, viral),
    ))
    _run("phylo", workdir, lambda: phylo_params(
        treebuilder=treebuilder, msa_source=msa_source, aligner=aligner,
        snptyper=snptyper, no_outgroup=no_outgroup, threads=threads,
    ))
    _run("tree2tax", workdir, lambda: tree2tax_params(
        include_dereplicated=include_dereplicated,
    ))

    typer.echo(
        f"\nPipeline complete. Deliverables in {workdir}: "
        "tree2tax.tsv, genomes_map.tsv, tree/tree.nwk"
    )
