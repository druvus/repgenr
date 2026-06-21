"""Phylogenetics commands: snptype, phylo, tree2tax."""

from __future__ import annotations

from pathlib import Path

import typer

from .base import _parse_key_values, _require_choice, _run, app


@app.command()
def snptype(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    tool: str = typer.Option("simple", "--tool", help="SNP typer: simple/snippy/parsnp."),
    reference: str | None = typer.Option(None, "--reference", help="Reference genome filename."),
    all_genomes: bool = typer.Option(False, "--all-genomes", help="Use all genomes, not reps."),
    mask: str = typer.Option("none", "--mask", help="Recombination masking: none or gubbins."),
    threads: int = typer.Option(16, "-t", "--threads"),
) -> None:
    """Call SNPs and build a core-SNP alignment."""
    from ..snptypers.base import registry as _snp_registry
    from ..stages.snptype import SnptypeParams

    def build() -> SnptypeParams:
        _require_choice(tool, set(_snp_registry.names()), "--tool")
        _require_choice(mask, {"none", "gubbins"}, "--mask")
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
        "iqtree", "--treebuilder", help="auto/iqtree/fasttree/raxmlng/mashtree/sourmash."
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
    aligner_arg: list[str] = typer.Option(
        [], "--aligner-arg",
        help="Aligner tuning as key=value (repeatable), e.g. kmer=15 (sibeliaz) "
        "or seed_weight=11 (progressivemauve).",
    ),
    threads: int = typer.Option(16, "-t", "--threads"),
) -> None:
    """Build a phylogenetic tree from an alignment, SNP alignment, or directly."""
    from ..aligners.base import registry as _aln_registry
    from ..snptypers.base import registry as _snp_registry
    from ..stages.phylo import PhyloParams
    from ..treebuilders.base import registry as _tb_registry

    def build() -> PhyloParams:
        _require_choice(treebuilder, {"auto", *_tb_registry.names()}, "--treebuilder")
        _require_choice(msa_source, {"aligner", "snptype"}, "--msa-source")
        if msa_source == "aligner":
            _require_choice(aligner, set(_aln_registry.names()), "--aligner")
        else:
            _require_choice(snptyper, set(_snp_registry.names()), "--snptyper")
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
            extra=_parse_key_values(aligner_arg, "--aligner-arg"),
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
