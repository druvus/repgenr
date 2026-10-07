"""derep_stock stage: store/load named dereplication runs.

Ports ``derep_stocker.py`` to the new ``derep/`` contract. A packed run keeps
``clusters.tsv`` + ``genome_status.tsv`` and symlinks the representative genome
files; unpacking restores them into the working directory.
"""

from __future__ import annotations

import re
import shutil
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import (
    CLUSTER_SUMMARY_TSV,
    CLUSTERS_TSV,
    GENOME_STATUS_TSV,
    list_fasta,
    read_clusters,
    read_genome_status,
    write_cluster_summary,
)
from ..core.errors import UserInputError, WorkdirError
from ..core.process import remove_tree
from ..dereplicators.base import DerepResult
from .dereplicate import _update_manifest

_FLAT_FILES = (CLUSTERS_TSV, GENOME_STATUS_TSV, CLUSTER_SUMMARY_TSV)
_NAME_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")


@dataclass
class DerepStockParams:
    action: str  # list | pack | unpack | delete
    name: str | None = None


def run(ctx: WorkdirContext, params: DerepStockParams) -> None:
    if not ctx.workdir.is_dir():
        raise WorkdirError(f"Working directory not found: {ctx.workdir}")
    store = ctx.derep_dir / "stock"
    if params.action == "list":
        _list(store, ctx.logger)
        return
    if not params.name:
        raise UserInputError("pack/unpack/delete require --name")
    # The name becomes a directory under the store and is passed to rmtree on
    # pack/delete, so it must be a plain single-component name.
    if not _NAME_RE.fullmatch(params.name):
        raise UserInputError(
            f"Invalid run name '{params.name}': use letters, digits, '.', '_' or '-' "
            "(no leading '.', no path separators)"
        )
    run_path = store / params.name
    match params.action:
        case "pack":
            _pack(ctx, run_path)
        case "unpack":
            _unpack(ctx, run_path)
        case "delete":
            # Not recorded: a delete leaves nothing to resume, and the CLI runs
            # it as a query so a repeat delete is checked instead of skipped.
            _delete(run_path)
            return
        case _:
            raise UserInputError(f"Unknown action '{params.action}'")
    ctx.config.record_stage(
        "derep_stock", params=asdict(params), completed=datetime.now(UTC).isoformat()
    )
    ctx.save_config()


def _list(store: Path, logger) -> None:
    if not store.exists() or not any(store.iterdir()):
        logger.info("No stored runs")
        return
    for run_dir in _stored_runs(store):
        logger.info(run_dir)


def _pack(ctx: WorkdirContext, run_path: Path) -> None:
    # Check before touching the store: a workdir without a dereplication
    # would otherwise be stored as an empty run.
    clusters = ctx.derep_dir / CLUSTERS_TSV
    if not clusters.is_file():
        raise WorkdirError(f"Missing {clusters}. Run the dereplicate stage first.")
    reps = list_fasta(ctx.representatives_dir)
    if not reps:
        raise WorkdirError(
            f"No representative genomes under {ctx.representatives_dir}. "
            "Run the dereplicate stage first."
        )
    if run_path.exists():
        remove_tree(run_path)
    run_path.mkdir(parents=True)
    for name in _FLAT_FILES:
        src = ctx.derep_dir / name
        if src.exists():
            shutil.copy2(src, run_path / name)
    reps_dir = run_path / "representatives"
    reps_dir.mkdir()
    for rep in reps:
        (reps_dir / rep.name).symlink_to((ctx.genomes_dir / rep.name).resolve())
    ctx.logger.info("Packed run to %s", run_path)


def _unpack(ctx: WorkdirContext, run_path: Path) -> None:
    if not run_path.exists():
        raise UserInputError(f"No stored run named '{run_path.name}'")
    # Validate the stored run in full before the current dereplication is
    # replaced, so an incomplete run leaves the workdir unchanged.
    stored_clusters = run_path / CLUSTERS_TSV
    if not stored_clusters.is_file():
        raise WorkdirError(f"Stored run '{run_path.name}' is incomplete: missing {stored_clusters}")
    stored_reps = run_path / "representatives"
    if not stored_reps.is_dir():
        raise WorkdirError(f"Stored run '{run_path.name}' is incomplete: missing {stored_reps}")
    rep_names = [rep.name for rep in list_fasta(stored_reps)]
    absent = [name for name in rep_names if not (ctx.genomes_dir / name).is_file()]
    if absent:
        raise WorkdirError(
            f"Stored run '{run_path.name}' names {len(absent)} representative(s) "
            f"not found under {ctx.genomes_dir}, e.g. {absent[0]}"
        )
    for name in _FLAT_FILES:
        src = run_path / name
        if src.exists():
            shutil.copy2(src, ctx.derep_dir / name)
        elif (ctx.derep_dir / name).exists():
            # Do not leave a file of the replaced dereplication beside the
            # restored ones; the summary is rebuilt below.
            (ctx.derep_dir / name).unlink()
    if ctx.representatives_dir.exists():
        remove_tree(ctx.representatives_dir)
    ctx.representatives_dir.mkdir(parents=True)
    for name in rep_names:
        shutil.copy2(ctx.genomes_dir / name, ctx.representatives_dir / name)
    # The derep contract now describes the stored run: bring the manifest's
    # per-genome status in line with it and re-stamp the dereplicate record
    # without a fingerprint, so `status` reports the run on disk and the next
    # `dereplicate` recomputes instead of skipping on a stale fingerprint.
    clusters = read_clusters(ctx.derep_dir / CLUSTERS_TSV)
    status_path = ctx.derep_dir / GENOME_STATUS_TSV
    genome_status = read_genome_status(status_path) if status_path.exists() else {}
    _update_manifest(
        ctx, DerepResult(representatives=[], clusters=clusters, genome_status=genome_status)
    )
    summary = ctx.derep_dir / CLUSTER_SUMMARY_TSV
    if not summary.exists():
        # A run packed before the summary existed: rebuild it from the
        # restored clusters, as the dereplicate stage would have written it.
        from .cluster_summary import summarise_clusters
        from .dereplicate import _quality_lookup

        write_cluster_summary(summary, summarise_clusters(clusters, _quality_lookup(ctx)))
        ctx.logger.info(
            "Stored run '%s' has no %s; rebuilt it from the restored clusters",
            run_path.name,
            CLUSTER_SUMMARY_TSV,
        )
    prior = ctx.config.stages.get("dereplicate")
    ctx.config.record_stage(
        "dereplicate",
        tool=prior.tool if prior else None,
        params={**(prior.params if prior else {}), "stock": run_path.name},
        tool_versions=prior.tool_versions if prior else None,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.logger.info("Unpacked run from %s", run_path)


def _stored_runs(store: Path) -> list[str]:
    if not store.is_dir():
        return []
    return sorted(p.name for p in store.iterdir() if p.is_dir())


def _delete(run_path: Path) -> None:
    if not run_path.exists():
        stored = _stored_runs(run_path.parent)
        raise WorkdirError(
            f"No stored run named '{run_path.name}'; stored runs: "
            f"{', '.join(stored) if stored else 'none'}."
        )
    remove_tree(run_path)
