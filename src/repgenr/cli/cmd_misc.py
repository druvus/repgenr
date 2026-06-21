"""Auxiliary commands: glance, derep-unpack, derep-stock, list-tools."""

from __future__ import annotations

from pathlib import Path

import typer

from .base import _run, app


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
