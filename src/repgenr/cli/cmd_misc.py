"""Auxiliary commands: status, doctor, glance, derep-unpack, derep-stock, list-tools."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import typer

from ..core.errors import UserInputError
from .base import (
    DEFAULT_THREADS,
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

    One 'tool: version' line per tool. A tool that stages recorded with
    different versions (an image in one, the host binary in another) is listed
    once per stage as 'tool (stage): version'. A containerized tool's version
    is its image reference. A stage that did not finish is named on stderr.

    Lets the Nextflow bridge modules (which run a full stage in a scratch
    workdir) surface the resolved tool versions into versions.yml.
    """
    from ..core.config import CONFIG_FILENAME, Config
    from ..core.errors import WorkdirError
    from ..core.versions import merge_stage_versions, write_versions_fragment

    if not (workdir / CONFIG_FILENAME).exists():
        # A wrong -wd would otherwise print nothing and exit 0.
        err = WorkdirError(f"No RepGenR run found at {workdir} (no {CONFIG_FILENAME}).")
        typer.echo(f"ERROR {err}", err=True)
        raise typer.Exit(code=err.exit_code)
    from ..stages.metadata import gtdb_provenance

    try:
        cfg = Config.load(workdir)
    except WorkdirError as exc:
        typer.echo(f"ERROR {exc}", err=True)
        raise typer.Exit(code=exc.exit_code) from exc
    for name, record in cfg.stages.items():
        if record.tool_versions and not record.completed:
            typer.echo(
                f"WARNING {name} did not finish; its versions are from its last recorded run.",
                err=True,
            )
    merged = merge_stage_versions({n: r.tool_versions for n, r in cfg.stages.items()})
    metadata_record = cfg.stages.get("metadata")
    if metadata_record is not None:
        merged.update(gtdb_provenance(metadata_record.params))
    if versions_out is not None:
        write_versions_fragment(versions_out, merged)
    else:
        for tool, ver in sorted(merged.items()):
            typer.echo(f"{tool}: {ver}")


def _gtdb_note(params: dict) -> str:
    """The GTDB release, or the date of the API query, after the metadata line."""
    from ..stages.metadata import gtdb_provenance

    gtdb = gtdb_provenance(params)
    if "gtdb_release" in gtdb:
        return f"  (GTDB release {gtdb['gtdb_release']})"
    if "gtdb_api_query_date" in gtdb:
        return f"  (GTDB API queried {gtdb['gtdb_api_query_date']})"
    return ""


@app.command(rich_help_panel=PANEL_PIPELINE)
def status(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
) -> None:
    """Show which pipeline stages have completed in a working directory.

    A completed stage is listed as stale when one of its inputs changed or
    one of its outputs is missing since it finished (it re-runs on its next
    invocation), and as interrupted when it did not finish. A -wd that does
    not exist exits 3; an existing directory without repgenr.yaml prints
    which entry stage to run first and exits 0.
    """
    from ..core.config import CONFIG_FILENAME, Config
    from ..core.errors import WorkdirError

    require_existing_workdir(workdir)
    if not (workdir / CONFIG_FILENAME).exists():
        typer.echo(f"No RepGenR run found at {workdir} (no {CONFIG_FILENAME}).")
        typer.echo(
            "Start with 'repgenr metadata' (bacteria), 'vmetadata' (viruses), "
            "'ingest' (local genomes) or 'reads' (sequencing reads), with -wd <wd>."
        )
        raise typer.Exit()

    try:
        cfg = Config.load(workdir)
    except WorkdirError as exc:
        typer.echo(f"ERROR {exc}", err=True)
        raise typer.Exit(code=exc.exit_code) from exc
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

    from ..core.doctor import stale_stages

    try:
        stale = stale_stages(workdir, cfg)
        unchecked = None
    except Exception as exc:  # a damaged artifact must not stop the report
        stale, unchecked = {}, str(exc)

    typer.echo(f"RepGenR workdir: {workdir}")
    typer.echo(f"Pipeline: {lineage}\n")

    next_stage: str | None = None
    for stage in chain:
        rec = recorded.get(stage)
        tool = f" [{rec.tool}]" if rec is not None and rec.tool else ""
        finished = rec is not None and not rec.interrupted
        note = (
            _gtdb_note(rec.params) if rec is not None and finished and stage == "metadata" else ""
        )
        if rec is not None and not rec.interrupted and stage not in stale:
            typer.echo(f"  [done]    {stage}{tool}  {rec.completed}{note}")
            continue
        if rec is not None and not rec.interrupted:
            # Completed, but an input changed or an output is missing since:
            # the stage re-runs on its next invocation.
            typer.echo(f"  [stale]   {stage}{tool}  {rec.completed}{note}  ({stale[stage]})")
        elif rec is not None:
            # The stage started a (re-)run and failed or was killed; outputs
            # may be partial.
            typer.echo(
                f"  [interrupted] {stage}  "
                "(did not finish; outputs may be partial; see repgenr.log)"
            )
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
            if rec.interrupted:
                when = "[interrupted] (did not finish; outputs may be partial; see repgenr.log)"
            elif stage in stale:
                when = f"{rec.completed}  [stale] ({stale[stage]})"
            else:
                when = rec.completed or ""
            typer.echo(f"    {stage}{tool}  {when}")

    if unchecked is not None:
        typer.echo(f"\nStale stages were not checked ({unchecked}); run repgenr doctor.")
    if next_stage is None:
        typer.echo("\nAll stages complete. Deliverables: tree2tax.tsv, genomes_map.tsv.")
    else:
        typer.echo(f"\nNext: repgenr {next_stage} -wd {workdir} ...")
        hint = _next_stage_note(workdir, next_stage, recorded)
        if hint:
            typer.echo(hint)


def _next_stage_note(workdir: Path, stage: str, recorded: dict[str, Any]) -> str | None:
    """A known refusal of the suggested stage, said before the user runs it."""
    if stage != "phylo" or "dereplicate" not in recorded:
        return None
    from ..core.contracts import list_fasta
    from ..stages.phylo import MIN_TREE_GENOMES

    phylo = recorded.get("phylo")
    if phylo is not None and phylo.params.get("all_genomes"):
        return None
    count = len(list_fasta(workdir / "derep" / "representatives"))
    if count >= MIN_TREE_GENOMES:
        return None
    return (
        f"Note: derep/representatives holds {count} genome(s) and a tree needs at least "
        f"{MIN_TREE_GENOMES}; run phylo with --all-genomes, or dereplicate again with a "
        "higher --secondary-ani."
    )


@app.command(rich_help_panel=PANEL_ENV)
def doctor(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
) -> None:
    """Verify a workdir's outputs against its records (read-only health check).

    `status` lists each stage as done, stale or interrupted; `doctor` also
    checks the outputs themselves: missing, corrupt or untracked genomes,
    dangling links, manifest drift, representative/cluster mismatches,
    truncated tree and tree2tax tables, missing deliverables, stages whose
    inputs changed since completion, and leftover temp files.
    Exits 0 when only warnings are found (a stale stage re-runs on its next
    invocation), 1 when any failure is found, and 3 when the workdir does
    not exist.
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
    from ..dereplicators.base import compare_supporters

    names = compare_supporters()
    return (
        f"Dereplicator with comparison support: auto, {', '.join(names) or '(none registered)'}. "
        "auto uses dRep when it can run (on the PATH or via the container backend), "
        "otherwise sourmash."
    )


@app.command(rich_help_panel=PANEL_INSPECT)
def glance(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    tool: str = typer.Option("auto", "--tool", help=_glance_tool_help()),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS),
    plot_max: float = typer.Option(
        1.0,
        "--plot-max",
        help="Upper bound of the ANI values plotted, as a fraction from 0 to 1.",
    ),
    plot_min: float = typer.Option(
        0.0,
        "--plot-min",
        help="Lower bound of the ANI values plotted, as a fraction from 0 to 1.",
    ),
    keep_files: bool = typer.Option(
        False,
        "--keep-files",
        help="Keep the comparison tool's working directory glance_wd/.",
    ),
) -> None:
    """Quick all-vs-all ANI overview (dendrogram + similarity plots).

    Needs dRep or sourmash (on the PATH, or via the container backend) and
    no dereplication; compares every genome in genomes/ and writes
    glance_clustering_dendrogram.pdf and two ANI plots,
    glance_MASH_ANI_similarity_boxplot.png and
    glance_MASH_ANI_similarity_histogram.png. dRep reports Mash ANI;
    sourmash reports the k-mer ANI estimate that dereplicate --tool
    sourmash thresholds.
    """
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.glance import GlanceParams, log_auto_choice, resolve_auto_tool

    def build() -> GlanceParams:
        _require_choice(tool, {"auto", *_derep_registry.names()}, "--tool")
        if not 0.0 <= plot_min <= plot_max <= 1.0:
            raise UserInputError(
                "--plot-min and --plot-max are ANI fractions with "
                f"0 <= --plot-min <= --plot-max <= 1; got {plot_min} and {plot_max}."
            )
        chosen = tool
        if tool == "auto":
            # Resolved here, before the resume fingerprint is taken, so the
            # fingerprint and the stage record name the concrete tool. When
            # nothing can run, 'auto' passes through and the stage reports it
            # after its workdir and genome checks.
            resolved = resolve_auto_tool()
            if resolved is not None:
                log_auto_choice(logging.getLogger("repgenr"), resolved)
                chosen = resolved
        return GlanceParams(
            tool=chosen,
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
    """Regenerate derep/cluster_summary.tsv (size, species, keeper quality per cluster).

    n_members counts the genomes under the representative and excludes it;
    n_genomes includes it. Species come from the manifest taxonomy, or from
    canonical filenames when the manifest has none, and at most five are
    listed.
    """
    from ..stages.cluster_summary import ClusterSummaryParams

    _run("cluster_summary", workdir, ClusterSummaryParams)


@app.command(name="derep-stock", rich_help_panel=PANEL_INSPECT)
def derep_stock(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    action: str = typer.Option(..., "--action", help="list, pack, unpack or delete."),
    name: str | None = typer.Option(
        None,
        "--name",
        help="Run name for pack/unpack/delete: up to 100 letters, digits, '.', '_' or '-', "
        "starting with a letter or digit.",
    ),
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
    The last line names the dereplicators that glance can run.
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
    from ..dereplicators.base import compare_supporters

    # glance is not a family of its own: it runs any dereplicator with compare().
    typer.echo(
        f"glance (dereplicators with compare): {', '.join(compare_supporters()) or '(none)'}"
    )


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
