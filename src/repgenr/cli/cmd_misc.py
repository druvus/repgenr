"""Auxiliary commands: status, doctor, glance, derep-unpack, derep-stock, list-tools."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import typer

from ..core.errors import UserInputError
from .base import (
    DEFAULT_THREADS,
    HELP_KEEP_FILES,
    HELP_THREADS,
    HELP_WORKDIR,
    PANEL_ENV,
    PANEL_INSPECT,
    PANEL_PIPELINE,
    PIPELINE_BACTERIAL,
    PIPELINE_LOCAL,
    PIPELINE_READS,
    PIPELINE_VIRAL,
    _require_choice,
    _run,
    app,
    require_existing_workdir,
)


@app.command(rich_help_panel=PANEL_ENV)
def versions(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    versions_out: Path | None = typer.Option(
        None, "--versions-out", help="Write a versions.yml fragment here instead of stdout."
    ),
) -> None:
    """Print the external-tool versions recorded in a workdir's repgenr.yaml.

    Lets the Nextflow bridge modules (which run a full stage in a scratch workdir)
    surface the resolved tool versions into versions.yml.
    """
    from ..core.config import CONFIG_FILENAME, Config
    from ..core.errors import WorkdirError
    from ..core.versions import write_versions_fragment

    if not (workdir / CONFIG_FILENAME).exists():
        # A wrong -wd would otherwise print nothing and exit 0.
        err = WorkdirError(f"No RepGenR run found at {workdir} (no {CONFIG_FILENAME}).")
        typer.echo(f"ERROR {err}", err=True)
        raise typer.Exit(code=err.exit_code)
    cfg = Config.load(workdir)
    merged: dict[str, str] = {}
    for record in cfg.stages.values():
        merged.update(record.tool_versions)
    if versions_out is not None:
        write_versions_fragment(versions_out, merged)
    else:
        for tool, ver in sorted(merged.items()):
            typer.echo(f"{tool}: {ver}")


@app.command(rich_help_panel=PANEL_PIPELINE)
def status(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
) -> None:
    """Show which pipeline stages have completed in a working directory.

    A -wd that does not exist exits 3; an existing directory without
    repgenr.yaml prints which entry stage to run first and exits 0.
    """
    from ..core.config import CONFIG_FILENAME, Config

    require_existing_workdir(workdir)
    if not (workdir / CONFIG_FILENAME).exists():
        typer.echo(f"No RepGenR run found at {workdir} (no {CONFIG_FILENAME}).")
        typer.echo(
            "Start with 'repgenr metadata' (bacteria), 'vmetadata' (viruses), "
            "'ingest' (local genomes) or 'reads' (sequencing reads), with -wd <wd>."
        )
        raise typer.Exit()

    cfg = Config.load(workdir)
    recorded = cfg.stages
    chain: tuple[str, ...]
    if "reads" in recorded and "metadata" not in recorded:
        lineage, chain = "reads", PIPELINE_READS
    elif "ingest" in recorded:
        lineage, chain = "local", PIPELINE_LOCAL
    elif any(name in recorded for name in ("vmetadata", "vgenome")):
        lineage, chain = "viral", PIPELINE_VIRAL
    else:
        lineage, chain = "bacterial", PIPELINE_BACTERIAL

    typer.echo(f"RepGenR workdir: {workdir}")
    typer.echo(f"Pipeline: {lineage}\n")

    next_stage: str | None = None
    for stage in chain:
        rec = recorded.get(stage)
        if rec is not None and rec.completed:
            tool = f" [{rec.tool}]" if rec.tool else ""
            typer.echo(f"  [done]    {stage}{tool}  {rec.completed}")
        elif rec is not None and (rec.params or rec.tool):
            # A record without a completed stamp but with provenance: the stage
            # started a (re-)run and failed or was killed; outputs may be partial.
            typer.echo(
                f"  [interrupted] {stage}  "
                "(did not finish; outputs may be partial; see repgenr.log)"
            )
            if next_stage is None:
                next_stage = stage
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


@app.command(rich_help_panel=PANEL_ENV)
def doctor(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
) -> None:
    """Verify a workdir's outputs against its records (read-only health check).

    `status` reports what repgenr.yaml claims; `doctor` checks the claims
    against the filesystem and the manifest: interrupted stages, missing or
    corrupt genomes, manifest drift, representative/cluster mismatches,
    truncated or missing deliverables, and stages whose inputs changed since
    completion.
    Exits 1 when any failure is found and 3 when the workdir does not exist.
    """
    from ..core.doctor import diagnose

    require_existing_workdir(workdir)

    findings = diagnose(workdir)
    label = {"ok": "OK  ", "warn": "WARN", "fail": "FAIL"}
    order = {"fail": 0, "warn": 1, "ok": 2}
    for f in sorted(findings, key=lambda f: (order[f.level], f.area)):
        typer.echo(f"[{label[f.level]}] {f.area}: {f.message}")
    failures = sum(1 for f in findings if f.level == "fail")
    warnings = sum(1 for f in findings if f.level == "warn")
    typer.echo(f"\n{failures} failure(s), {warnings} warning(s).")
    if failures:
        raise typer.Exit(code=1)


def _glance_tool_help() -> str:
    """Dereplicators that implement the comparison capability, from the registry."""
    from ..dereplicators.base import Dereplicator, registry

    names = [
        n
        for n in registry.names()
        if not registry.is_broken(n) and registry.get(n).compare is not Dereplicator.compare
    ]
    return f"Dereplicator with comparison support: {', '.join(names) or '(none registered)'}."


@app.command(rich_help_panel=PANEL_INSPECT)
def glance(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    tool: str = typer.Option("drep", "--tool", help=_glance_tool_help()),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS),
    plot_max: float = typer.Option(
        1.0, "--plot-max", help="Upper similarity bound of the values plotted."
    ),
    plot_min: float = typer.Option(
        0.0, "--plot-min", help="Lower similarity bound of the values plotted."
    ),
    keep_files: bool = typer.Option(False, "--keep-files", help=HELP_KEEP_FILES),
) -> None:
    """Quick all-vs-all ANI overview (dRep compare dendrogram + plots)."""
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.glance import GlanceParams

    def build() -> GlanceParams:
        _require_choice(tool, set(_derep_registry.names()), "--tool")
        if not 0.0 <= plot_min <= plot_max <= 1.0:
            raise UserInputError(
                "--plot-min and --plot-max are Mash ANI fractions with "
                f"0 <= --plot-min <= --plot-max <= 1; got {plot_min} and {plot_max}."
            )
        return GlanceParams(
            tool=tool,
            threads=threads,
            plot_max=plot_max,
            plot_min=plot_min,
            keep_files=keep_files,
        )

    _run("glance", workdir, build)


@app.command(name="derep-unpack", rich_help_panel=PANEL_INSPECT)
def derep_unpack(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    no_representant: bool = typer.Option(
        False, "--no-representant", help="Leave the representative out of its cluster directory."
    ),
) -> None:
    """Explode clusters into one directory per representative."""
    from ..stages.derep_unpack import DerepUnpackParams

    def build() -> DerepUnpackParams:
        return DerepUnpackParams(no_representant=no_representant)

    _run("derep_unpack", workdir, build)


@app.command(name="cluster-summary", rich_help_panel=PANEL_INSPECT)
def cluster_summary(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
) -> None:
    """Regenerate derep/cluster_summary.tsv (size, species, keeper quality per cluster)."""
    from ..stages.cluster_summary import ClusterSummaryParams

    _run("cluster_summary", workdir, ClusterSummaryParams)


@app.command(name="derep-stock", rich_help_panel=PANEL_INSPECT)
def derep_stock(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    action: str = typer.Option(..., "--action", help="list, pack, unpack or delete."),
    name: str | None = typer.Option(None, "--name", help="Run name for pack/unpack/delete."),
) -> None:
    """Store, load, list or delete named dereplication runs."""
    from ..stages.derep_stock import DerepStockParams
    from .param_builders import DEREP_STOCK_ACTIONS

    def build() -> DerepStockParams:
        _require_choice(action, DEREP_STOCK_ACTIONS, "--action")
        return DerepStockParams(action=action, name=name)

    _run("derep_stock", workdir, build)


def _tool_label(reg: Any, name: str) -> str:
    """The name, plus the adapter's declared genome limit when it has one."""
    from ..core.plugins import _capabilities_of

    if reg.is_broken(name):
        return f"{name} (broken)"
    cap = _capabilities_of(reg, name)
    limit = None if cap is None else cap.recommended_max_genomes
    return name if limit is None else f"{name} (up to {limit} genomes)"


@app.command(name="list-tools", rich_help_panel=PANEL_ENV)
def list_tools(
    check: bool = typer.Option(
        False,
        "--check",
        help="Run each adapter's preflight and report whether its binaries are present.",
    ),
) -> None:
    """List the available pluggable tools in each family.

    A tool that declares a recommended scale is shown as 'name (up to N
    genomes)'; auto-selection and the scale warnings use the same limit.
    With --check, every adapter's required binaries are looked up (version
    floors included) and reported per tool, so an environment can be
    verified before a run without a working directory.
    """
    from ..aligners.base import registry as aligners
    from ..assemblers.base import registry as assemblers
    from ..classifiers.base import registry as classifiers
    from ..dereplicators.base import registry as dereplicators
    from ..maskers.base import registry as maskers
    from ..polishers.base import registry as polishers
    from ..snptypers.base import registry as snptypers
    from ..treebuilders.base import registry as treebuilders

    for label, reg in (
        ("dereplicators", dereplicators),
        ("aligners", aligners),
        ("snptypers", snptypers),
        ("maskers", maskers),
        ("treebuilders", treebuilders),
        ("assemblers", assemblers),
        ("classifiers", classifiers),
        ("polishers", polishers),
    ):
        entries = [_tool_label(reg, name) for name in reg.names()]
        typer.echo(f"{label}: {', '.join(entries) or '(none)'}")
        if not check:
            continue
        for name in reg.names():
            typer.echo(f"  {name}: {_preflight_summary(reg, name)}")


def _one_line(exc: Exception) -> str:
    """Flatten a multi-line preflight message to one line, dropping its preamble."""
    text = str(exc).replace("Required external tools are missing or outdated:", "")
    return "; ".join(part.strip() for part in text.splitlines() if part.strip())


def _preflight_summary(reg, name: str) -> str:
    """One line per adapter for `list-tools --check`: ok with versions, or why not."""
    from ..core.errors import MissingBinaryError, RepGenRError

    if reg.is_broken(name):
        return "broken (see list-tools)"
    try:
        versions = reg.create(name).preflight()
    except MissingBinaryError as exc:
        return f"missing ({_one_line(exc)})"
    except RepGenRError as exc:
        return f"error ({_one_line(exc)})"
    shown = ", ".join(f"{k} {v}" for k, v in sorted(versions.items())) or "no binaries declared"
    return f"ok ({shown})"
