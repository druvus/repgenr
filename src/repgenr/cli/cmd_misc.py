"""Auxiliary commands: status, glance, derep-unpack, derep-stock, list-tools."""

from __future__ import annotations

from pathlib import Path

import typer

from .base import DEFAULT_THREADS, PIPELINE_BACTERIAL, PIPELINE_VIRAL, _run, app


@app.command()
def versions(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    versions_out: Path | None = typer.Option(
        None, "--versions-out", help="Write a versions.yml fragment here instead of stdout."
    ),
) -> None:
    """Print the external-tool versions recorded in a workdir's repgenr.yaml.

    Lets the Nextflow bridge modules (which run a full stage in a scratch workdir)
    surface the resolved tool versions into versions.yml.
    """
    from ..core.config import Config
    from ..core.versions import write_versions_fragment

    cfg = Config.load(workdir)
    merged: dict[str, str] = {}
    for record in cfg.stages.values():
        merged.update(record.tool_versions)
    if versions_out is not None:
        write_versions_fragment(versions_out, merged)
    else:
        for tool, ver in sorted(merged.items()):
            typer.echo(f"{tool}: {ver}")


@app.command()
def status(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
) -> None:
    """Show which pipeline stages have completed in a working directory."""
    from ..core.config import CONFIG_FILENAME, Config

    if not (workdir / CONFIG_FILENAME).exists():
        typer.echo(f"No RepGenR run found at {workdir} (no {CONFIG_FILENAME}).")
        typer.echo("Start with 'repgenr metadata -wd <wd> ...' (or 'vmetadata' for viruses).")
        raise typer.Exit()

    cfg = Config.load(workdir)
    recorded = cfg.stages
    viral = any(name in recorded for name in ("vmetadata", "vgenome"))
    chain = PIPELINE_VIRAL if viral else PIPELINE_BACTERIAL

    typer.echo(f"RepGenR workdir: {workdir}")
    typer.echo(f"Pipeline: {'viral' if viral else 'bacterial'}\n")

    next_stage: str | None = None
    for stage in chain:
        rec = recorded.get(stage)
        if rec is not None and rec.completed:
            tool = f" [{rec.tool}]" if rec.tool else ""
            typer.echo(f"  [done]    {stage}{tool}  {rec.completed}")
        else:
            marker = "next" if next_stage is None else "    "
            typer.echo(f"  [{marker}] {stage}")
            if next_stage is None:
                next_stage = stage

    extras = [s for s in recorded if s not in chain]
    if extras:
        typer.echo("\n  optional stages run:")
        for stage in extras:
            rec = recorded[stage]
            tool = f" [{rec.tool}]" if rec.tool else ""
            typer.echo(f"    {stage}{tool}  {rec.completed or '(incomplete)'}")

    if next_stage is None:
        typer.echo("\nAll stages complete. Deliverables: tree2tax.tsv, genomes_map.tsv.")
    else:
        typer.echo(f"\nNext: repgenr {next_stage} -wd {workdir} ...")


@app.command()
def glance(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads"),
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
