"""RepGenR command-line interface.

Replaces the old ``repgenr.py`` string-rewriting dispatcher with a real Typer
app. Each subcommand parses arguments, builds a :class:`WorkdirContext`, and
calls the matching ``stages.<name>.run(ctx, params)``. Errors of type
:class:`RepGenRError` are logged cleanly and turned into a non-zero exit code.
"""

from __future__ import annotations

import sys
from pathlib import Path

import typer

from .. import __version__
from ..core.context import WorkdirContext
from ..core.errors import RepGenRError
from ..core.logging import configure_logging

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="RepGenR: modular genome dereplication, alignment, SNP typing and phylogenetics.",
)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"repgenr {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        False, "--version", callback=_version_callback, is_eager=True,
        help="Show version and exit.",
    ),
) -> None:
    """RepGenR top-level entry point."""


def _run(stage_name: str, workdir: Path, build_params, *, create: bool = False) -> None:
    """Common harness: context, dispatch, clean error handling."""
    logger = configure_logging(workdir if (create or workdir.exists()) else None)
    try:
        ctx = WorkdirContext(workdir, logger=logger, create=create)
        module = __import__(f"repgenr.stages.{stage_name}", fromlist=["run"])
        module.run(ctx, build_params())
    except RepGenRError as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=1) from exc


@app.command()
def metadata(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory (created)."),
    release: str = typer.Option(..., "-r", "--release", help="GTDB release, e.g. 207.0."),
    version: str = typer.Option(..., "-v", "--version", help="bac120 or ar53."),
    dataset: str = typer.Option(..., "-d", "--dataset", help="all or rep."),
    level: str = typer.Option(..., "-l", "--level", help="family, genus or species."),
    target_family: str | None = typer.Option(None, "-tf", "--target-family"),
    target_genus: str | None = typer.Option(None, "-tg", "--target-genus"),
    target_species: str | None = typer.Option(None, "-ts", "--target-species"),
    outgroup_accession: str | None = typer.Option(None, "--outgroup-accession"),
    metadata_path: str | None = typer.Option(None, "--metadata-path"),
    nodownload: bool = typer.Option(False, "--nodownload"),
    limit: int | None = typer.Option(None, "--limit"),
) -> None:
    """Select a taxon's genomes from GTDB metadata."""
    from ..stages.metadata import MetadataParams

    def build() -> MetadataParams:
        return MetadataParams(
            release=release, version=version, dataset=dataset, level=level,
            target_family=target_family, target_genus=target_genus,
            target_species=target_species, outgroup_accession=outgroup_accession,
            metadata_path=metadata_path, nodownload=nodownload, limit=limit,
        )

    _run("metadata", workdir, build, create=True)


@app.command()
def vmetadata(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory (created)."),
    target: str | None = typer.Option(None, "-t", "--target", help="Virus group/family."),
    filter: str = typer.Option("complete genome", "-f", "--filter", help="FASTA header tag."),
    list_targets: bool = typer.Option(False, "-l", "--list", help="List BV-BRC targets and exit."),
) -> None:
    """Retrieve viral metadata from BV-BRC and NCBI (virus equivalent of metadata)."""
    from ..stages.vmetadata import VmetadataParams

    def build() -> VmetadataParams:
        return VmetadataParams(target=target, filter=filter, list_targets=list_targets)

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
            no_outgroup=no_outgroup,
            outgroup_candidates_taxid_min_genomes=min_outgroup_genomes,
            glance=glance, print_fasta_headers=print_fasta_headers,
            ignore_duplicates=ignore_duplicates, keep_files=keep_files,
        )

    _run("vgenome", workdir, build)


@app.command()
def genome(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    accession_list_only: bool = typer.Option(False, "--accession-list-only"),
    keep_files: bool = typer.Option(False, "--keep-files"),
) -> None:
    """Download and organize genomes selected by the metadata stage."""
    from ..stages.genome import GenomeParams

    def build() -> GenomeParams:
        return GenomeParams(accession_list_only=accession_list_only, keep_files=keep_files)

    _run("genome", workdir, build)


@app.command()
def dereplicate(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    tool: str = typer.Option("skder", "--tool", help="Dereplicator: skder, drep, galah, sourmash."),
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani"),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani"),
    aligned_fraction: float = typer.Option(0.50, "-af", "--aligned-fraction"),
    threads: int = typer.Option(16, "-t", "--threads"),
    process_size: int | None = typer.Option(
        None, "-s", "--process-size", help="Chunk size for tools that don't scale natively."
    ),
    virus: bool = typer.Option(False, "--virus", help="Pass virus-tuned parameters to the tool."),
) -> None:
    """Cluster genomes by ANI and select representatives."""
    from ..stages.dereplicate import DereplicateParams

    def build() -> DereplicateParams:
        return DereplicateParams(
            tool=tool,
            primary_ani=primary_ani,
            secondary_ani=secondary_ani,
            aligned_fraction=aligned_fraction,
            threads=threads,
            process_size=process_size,
            extra={"virus": virus} if virus else {},
        )

    _run("dereplicate", workdir, build)


@app.command()
def snptype(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    tool: str = typer.Option("simple", "--tool", help="SNP typer: simple/snippy/parsnp/ksnp."),
    reference: str | None = typer.Option(None, "--reference", help="Reference genome filename."),
    all_genomes: bool = typer.Option(False, "--all-genomes", help="Use all genomes, not reps."),
    mask: str = typer.Option("none", "--mask", help="Recombination masking: none or gubbins."),
    threads: int = typer.Option(16, "-t", "--threads"),
) -> None:
    """Call SNPs and build a core-SNP alignment."""
    from ..stages.snptype import SnptypeParams

    def build() -> SnptypeParams:
        return SnptypeParams(
            tool=tool,
            threads=threads,
            reference=reference,
            all_genomes=all_genomes,
            mask=mask,
        )

    _run("snptype", workdir, build)


@app.command()
def phylo(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    treebuilder: str = typer.Option(
        "iqtree", "--treebuilder", help="iqtree, fasttree, raxmlng, mashtree, sourmash."
    ),
    msa_source: str = typer.Option("aligner", "--msa-source", help="aligner or snptype."),
    aligner: str = typer.Option(
        "progressivemauve", "--aligner", help="progressivemauve, cactus, sibeliaz."
    ),
    snptyper: str = typer.Option("simple", "--snptyper", help="SNP typer for snptype source."),
    all_genomes: bool = typer.Option(False, "--all-genomes", help="Use all genomes, not reps."),
    no_outgroup: bool = typer.Option(False, "--no-outgroup", help="Do not root with an outgroup."),
    bootstrap: int = typer.Option(0, "-B", "--bootstrap", help="Bootstrap replicates (>=1000)."),
    reference: str | None = typer.Option(None, "--reference", help="Reference genome filename."),
    threads: int = typer.Option(16, "-t", "--threads"),
) -> None:
    """Build a phylogenetic tree from an alignment, SNP alignment, or directly."""
    from ..stages.phylo import PhyloParams

    def build() -> PhyloParams:
        return PhyloParams(
            treebuilder=treebuilder,
            msa_source=msa_source,
            aligner=aligner,
            snptyper=snptyper,
            all_genomes=all_genomes,
            no_outgroup=no_outgroup,
            bootstrap=bootstrap,
            reference=reference,
            threads=threads,
        )

    _run("phylo", workdir, build)


@app.command()
def tree2tax(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    node_basename: str | None = typer.Option(None, "--node-basename", help="Prefix for nodes."),
    root_name: str = typer.Option("root", "-r", "--root-name", help="Name for the root node."),
    remove_outgroup: bool = typer.Option(False, "--remove-outgroup", help="Drop outgroup."),
    include_dereplicated: bool = typer.Option(
        False, "--include-dereplicated", help="List redundant genomes under their representative."
    ),
) -> None:
    """Emit FlexTaxD-compatible taxonomy relations from the tree."""
    from ..stages.tree2tax import Tree2taxParams

    def build() -> Tree2taxParams:
        return Tree2taxParams(
            node_basename=node_basename,
            root_name=root_name,
            remove_outgroup=remove_outgroup,
            include_dereplicated=include_dereplicated,
        )

    _run("tree2tax", workdir, build)


@app.command()
def glance(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    threads: int = typer.Option(24, "-t", "--threads"),
    plot_max: float = typer.Option(1.0, "--plot-max"),
    plot_min: float = typer.Option(0.0, "--plot-min"),
    keep_files: bool = typer.Option(False, "--keep-files"),
) -> None:
    """Quick all-vs-all ANI overview (dRep compare dendrogram + plots)."""
    from ..stages.glance import GlanceParams

    def build() -> GlanceParams:
        return GlanceParams(
            threads=threads, plot_max=plot_max, plot_min=plot_min, keep_files=keep_files
        )

    _run("glance", workdir, build)


@app.command(name="derep-unpack")
def derep_unpack(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    no_representant: bool = typer.Option(False, "--no-representant"),
) -> None:
    """Explode clusters into one directory per representative."""
    from ..stages.derep_unpack import DerepUnpackParams

    def build() -> DerepUnpackParams:
        return DerepUnpackParams(no_representant=no_representant)

    _run("derep_unpack", workdir, build)


@app.command(name="derep-stock")
def derep_stock(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    action: str = typer.Option(..., "--action", help="list, pack, unpack or delete."),
    name: str | None = typer.Option(None, "--name", help="Run name for pack/unpack/delete."),
) -> None:
    """Store, load, list or delete named dereplication runs."""
    from ..stages.derep_stock import DerepStockParams

    def build() -> DerepStockParams:
        return DerepStockParams(action=action, name=name)

    _run("derep_stock", workdir, build)


@app.command(name="list-tools")
def list_tools() -> None:
    """List the available pluggable tools in each family."""
    from ..aligners.base import registry as aligners
    from ..dereplicators.base import registry as dereplicators
    from ..snptypers.base import registry as snptypers
    from ..treebuilders.base import registry as treebuilders

    for label, reg in (
        ("dereplicators", dereplicators),
        ("aligners", aligners),
        ("snptypers", snptypers),
        ("treebuilders", treebuilders),
    ):
        typer.echo(f"{label}: {', '.join(reg.names()) or '(none)'}")


if __name__ == "__main__":  # pragma: no cover
    sys.exit(app())
