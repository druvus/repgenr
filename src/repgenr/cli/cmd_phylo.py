"""Phylogenetics commands: snptype, phylo, tree2tax."""

from __future__ import annotations

from pathlib import Path

import typer

from .base import (
    DEFAULT_THREADS,
    _aligner_help,
    _mask_help,
    _parse_key_values,
    _require_choice,
    _run,
    _snp_help,
    _tree_help,
    app,
)


@app.command()
def snptype(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    tool: str = typer.Option("simple", "--tool", help=_snp_help()),
    reference: str | None = typer.Option(None, "--reference", help="Reference genome filename."),
    all_genomes: bool = typer.Option(False, "--all-genomes", help="Use all genomes, not reps."),
    mask: str = typer.Option("none", "--mask", help=_mask_help()),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", min=1),
    tool_arg: list[str] = typer.Option(
        [], "--tool-arg", help="Tool tuning as key=value (repeatable)."
    ),
    allow_incomplete: bool = typer.Option(
        False, "--allow-incomplete",
        help="Proceed with a warning when the input genome set is incomplete.",
    ),
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
            allow_incomplete=allow_incomplete,
            extra=_parse_key_values(tool_arg, "--tool-arg"),
        )

    _run("snptype", workdir, build)


@app.command()
def phylo(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    treebuilder: str = typer.Option("iqtree", "--treebuilder", help=_tree_help()),
    msa_source: str = typer.Option("aligner", "--msa-source", help="aligner or snptype."),
    aligner: str = typer.Option("progressivemauve", "--aligner", help=_aligner_help()),
    snptyper: str = typer.Option("simple", "--snptyper", help=_snp_help()),
    all_genomes: bool = typer.Option(False, "--all-genomes", help="Use all genomes, not reps."),
    no_outgroup: bool = typer.Option(False, "--no-outgroup", help="Do not root with an outgroup."),
    bootstrap: int = typer.Option(
        0, "-B", "--bootstrap", min=0, help="Bootstrap replicates (>=1000)."
    ),
    reference: str | None = typer.Option(None, "--reference", help="Reference genome filename."),
    aligner_arg: list[str] = typer.Option(
        [], "--aligner-arg",
        help="Aligner tuning as key=value (repeatable), e.g. kmer=15 (sibeliaz) "
        "or seed_weight=11 (progressivemauve).",
    ),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", min=1),
    mask: str = typer.Option(
        "none", "--mask", help="Recombination masking for --msa-source snptype.",
    ),
    allow_incomplete: bool = typer.Option(
        False, "--allow-incomplete",
        help="Proceed with a warning when the input genome set is incomplete.",
    ),
) -> None:
    """Build a phylogenetic tree from an alignment, SNP alignment, or directly."""
    from .param_builders import phylo_params

    def build():
        return phylo_params(
            treebuilder=treebuilder,
            msa_source=msa_source,
            aligner=aligner,
            snptyper=snptyper,
            all_genomes=all_genomes,
            no_outgroup=no_outgroup,
            bootstrap=bootstrap,
            reference=reference,
            threads=threads,
            extra={
                **_parse_key_values(aligner_arg, "--aligner-arg"),
                **({"mask": mask} if mask != "none" else {}),
            },
            allow_incomplete=allow_incomplete,
        )

    _run("phylo", workdir, build)


@app.command()
def tree2tax(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    node_basename: str | None = typer.Option(None, "--node-basename", help="Prefix for nodes."),
    root_name: str = typer.Option("root", "-r", "--root-name", help="Name for the root node."),
    remove_outgroup: bool = typer.Option(False, "--remove-outgroup", help="Drop outgroup."),
    include_dereplicated: bool = typer.Option(
        True, "--include-dereplicated/--no-include-dereplicated",
        help="List redundant genomes under their representative.",
    ),
) -> None:
    """Emit FlexTaxD-compatible taxonomy relations from the tree."""
    from .param_builders import tree2tax_params

    def build():
        return tree2tax_params(
            node_basename=node_basename,
            root_name=root_name,
            remove_outgroup=remove_outgroup,
            include_dereplicated=include_dereplicated,
        )

    _run("tree2tax", workdir, build)
