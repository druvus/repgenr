"""derep_stock stage: store/load named dereplication runs.

Ports ``derep_stocker.py`` to the new ``derep/`` contract. A packed run keeps
``clusters.tsv`` + ``genome_status.tsv`` and symlinks the representative genome
files; unpacking restores them into the working directory.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import CLUSTERS_TSV, GENOME_STATUS_TSV
from ..core.errors import UserInputError

_FLAT_FILES = (CLUSTERS_TSV, GENOME_STATUS_TSV)


@dataclass
class DerepStockParams:
    action: str  # list | pack | unpack | delete
    name: str | None = None


def run(ctx: WorkdirContext, params: DerepStockParams) -> None:
    store = ctx.derep_dir / "stock"
    if params.action == "list":
        _list(store, ctx.logger)
        return
    if not params.name:
        raise UserInputError("pack/unpack/delete require --name")
    run_path = store / params.name
    match params.action:
        case "pack":
            _pack(ctx, run_path)
        case "unpack":
            _unpack(ctx, run_path)
        case "delete":
            _delete(run_path)
        case _:
            raise UserInputError(f"Unknown action '{params.action}'")


def _list(store: Path, logger) -> None:
    if not store.exists() or not any(store.iterdir()):
        logger.info("No stored runs")
        return
    for run_dir in sorted(p.name for p in store.iterdir() if p.is_dir()):
        logger.info(run_dir)


def _pack(ctx: WorkdirContext, run_path: Path) -> None:
    if run_path.exists():
        shutil.rmtree(run_path)
    run_path.mkdir(parents=True)
    for name in _FLAT_FILES:
        src = ctx.derep_dir / name
        if src.exists():
            shutil.copy2(src, run_path / name)
    reps = run_path / "representatives"
    reps.mkdir()
    for rep in ctx.representatives_dir.iterdir():
        (reps / rep.name).symlink_to((ctx.genomes_dir / rep.name).resolve())
    ctx.logger.info("Packed run to %s", run_path)


def _unpack(ctx: WorkdirContext, run_path: Path) -> None:
    if not run_path.exists():
        raise UserInputError(f"No stored run named '{run_path.name}'")
    for name in _FLAT_FILES:
        src = run_path / name
        if src.exists():
            shutil.copy2(src, ctx.derep_dir / name)
    if ctx.representatives_dir.exists():
        shutil.rmtree(ctx.representatives_dir)
    ctx.representatives_dir.mkdir(parents=True)
    for rep in (run_path / "representatives").iterdir():
        shutil.copy2(ctx.genomes_dir / rep.name, ctx.representatives_dir / rep.name)
    ctx.logger.info("Unpacked run from %s", run_path)


def _delete(run_path: Path) -> None:
    if not run_path.exists():
        raise UserInputError(f"No stored run named '{run_path.name}'")
    shutil.rmtree(run_path)
