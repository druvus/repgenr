"""Viral lineage commands: vmetadata, vgenome."""

from __future__ import annotations

from pathlib import Path

import typer

from .base import _require_choice, _run, app


@app.command()
def vmetadata(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory (created)."),
    target: str | None = typer.Option(None, "-t", "--target", help="Virus taxon/group/family."),
    source: str = typer.Option(
        "ncbi_virus", "--source", help="ncbi_virus (NCBI Virus via datasets) or bvbrc."
    ),
    filter: str = typer.Option("complete genome", "-f", "--filter", help="BV-BRC header tag."),
    host: str | None = typer.Option(None, "--host", help="ncbi_virus: restrict to a host species."),
    complete_only: bool = typer.Option(
        False, "--complete-only", help="ncbi_virus: only COMPLETE sequences."
    ),
    released_after: str | None = typer.Option(
        None, "--released-after", help="ncbi_virus: MM/DD/YYYY."
    ),
    list_targets: bool = typer.Option(False, "-l", "--list", help="List BV-BRC targets and exit."),
) -> None:
    """Retrieve viral metadata from NCBI Virus (default) or BV-BRC."""
    from ..stages.vmetadata import VmetadataParams

    def build() -> VmetadataParams:
        _require_choice(source, {"ncbi_virus", "bvbrc"}, "--source")
        return VmetadataParams(
            target=target, filter=filter, list_targets=list_targets,
            source=source, host=host, complete_only=complete_only,
            released_after=released_after,
        )

    _run("vmetadata", workdir, build, create=True)


@app.command()
def vgenome(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    target_genus: str | None = typer.Option(None, "-tg", "--target-genus"),
    target_species: str | None = typer.Option(None, "-ts", "--target-species"),
    target_serotype: str | None = typer.Option(None, "-tse", "--target-serotype"),
    target_custom: str | None = typer.Option(None, "-tc", "--target-custom", help="key:value."),
    length_all: bool = typer.Option(False, "--length-all"),
    length_deviation: int = typer.Option(10, "--length-deviation"),
    length_method: str = typer.Option("median_of_medians", "--length-method"),
    length_range: str | None = typer.Option(None, "--length-range", help="e.g. 25000-35000."),
    discard: str | None = typer.Option(None, "--discard", help="Comma-separated header tags."),
    no_outgroup: bool = typer.Option(False, "--no-outgroup"),
    group_segments: bool = typer.Option(
        False, "--group-segments",
        help="ncbi_virus: combine an isolate's segments into one genome (segmented viruses).",
    ),
    min_outgroup_genomes: int = typer.Option(5, "--outgroup-candidates-taxid-min-genomes"),
    glance: bool = typer.Option(False, "--glance", help="Print selection and stop."),
    print_fasta_headers: bool = typer.Option(False, "--print-fasta-headers"),
    ignore_duplicates: bool = typer.Option(False, "--ignore-duplicates"),
    keep_files: bool = typer.Option(False, "--keep-files"),
) -> None:
    """Select and organize viral genomes (virus equivalent of genome)."""
    from ..stages.vgenome import VgenomeParams

    def build() -> VgenomeParams:
        return VgenomeParams(
            target_genus=target_genus, target_species=target_species,
            target_serotype=target_serotype, target_custom=target_custom,
            length_all=length_all, length_deviation=length_deviation,
            length_method=length_method, length_range=length_range, discard=discard,
            no_outgroup=no_outgroup, group_segments=group_segments,
            outgroup_candidates_taxid_min_genomes=min_outgroup_genomes,
            glance=glance, print_fasta_headers=print_fasta_headers,
            ignore_duplicates=ignore_duplicates, keep_files=keep_files,
        )

    _run("vgenome", workdir, build)
