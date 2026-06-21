"""Shared CLI app, callback and helpers.

The Typer ``app``, the top-level callback (container/logging setup) and the
common stage harness (:func:`_run`) live here so the per-domain command modules
(``cmd_*.py``) can register against a single app without circular imports.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import logging
import os
from pathlib import Path
from typing import Any

import typer

from .. import __version__
from ..core.context import WorkdirContext
from ..core.errors import RepGenRError, UserInputError
from ..core.logging import configure_logging

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="RepGenR: modular genome dereplication, alignment, SNP typing and phylogenetics.",
)

# Top-level run options shared by every subcommand (set in the callback).
_RUN_STATE: dict[str, Any] = {"force": False, "log_level": logging.INFO}


def _require_choice(value: str, choices: set[str], label: str) -> None:
    if value not in choices:
        raise UserInputError(
            f"Invalid {label} {value!r}. Choose from: {', '.join(sorted(choices))}."
        )


def _require_unit_interval(value: float | None, label: str) -> None:
    if value is not None and not (0.0 < value <= 1.0):
        raise UserInputError(f"{label} must be in (0, 1], got {value}.")


def _stage_fingerprint(stage_name: str, params: object) -> str:
    """Stable hash of a stage invocation, used to skip already-completed work.

    Built from the stage name plus the parameter object (a dataclass), so the
    same invocation produces the same fingerprint across runs. Paths and other
    non-JSON values are stringified.
    """
    if dataclasses.is_dataclass(params) and not isinstance(params, type):
        payload: object = dataclasses.asdict(params)
    else:
        payload = vars(params)
    blob = json.dumps({"stage": stage_name, "params": payload}, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


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
    container: str = typer.Option(
        "none", "--container", envvar="REPGENR_CONTAINER",
        help="Run external tools in containers: none, docker, or singularity.",
    ),
    container_engine: str | None = typer.Option(
        None, "--container-engine", envvar="REPGENR_CONTAINER_ENGINE",
        help="Engine binary override (e.g. apptainer, podman).",
    ),
    container_cache: str | None = typer.Option(
        None, "--container-cache", envvar="REPGENR_CONTAINER_CACHE",
        help="Directory for Singularity .sif images / Wave cache (large; can be external).",
    ),
    platform: str | None = typer.Option(
        None, "--platform", envvar="REPGENR_CONTAINER_PLATFORM",
        help="Container platform, e.g. linux/amd64 for emulated BioContainers on arm64.",
    ),
    wave: bool = typer.Option(
        False, "--wave/--no-wave", envvar="REPGENR_WAVE",
        help="Resolve images for multi-tool adapters via the Seqera Wave CLI.",
    ),
    force: bool = typer.Option(
        False, "--force/--no-force", "-f", envvar="REPGENR_FORCE",
        help="Re-run a stage even if it already completed with the same parameters.",
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose (DEBUG) logging."),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Only warnings and errors."),
) -> None:
    """RepGenR top-level entry point."""
    from ..core.containers import configure_container

    _RUN_STATE["force"] = force
    if quiet:
        level = logging.WARNING
    elif verbose:
        level = logging.DEBUG
    else:
        env = os.environ.get("REPGENR_LOG_LEVEL")
        level = getattr(logging, env.upper(), logging.INFO) if env else logging.INFO
    _RUN_STATE["log_level"] = level
    configure_container(
        backend=container, engine=container_engine, platform=platform,
        cache_dir=container_cache, wave_enabled=wave,
    )


def _run(stage_name: str, workdir: Path, build_params, *, create: bool = False) -> None:
    """Common harness: context, dispatch, clean error handling.

    Resume: a stage that already completed with the same parameters is skipped
    (fingerprint match), unless ``--force`` is set. A stage that crashed before
    recording completion has no ``completed`` stamp and so always re-runs.
    """
    logger = configure_logging(
        workdir if (create or workdir.exists()) else None, level=_RUN_STATE["log_level"]
    )
    try:
        ctx = WorkdirContext(workdir, logger=logger, create=create)
        params = build_params()
        fingerprint = _stage_fingerprint(stage_name, params)
        prior = ctx.config.stages.get(stage_name)
        if (
            not _RUN_STATE["force"]
            and prior is not None
            and prior.completed
            and prior.fingerprint == fingerprint
        ):
            logger.info(
                "Stage '%s' already completed with the same parameters; skipping "
                "(use --force to re-run).", stage_name,
            )
            return
        module = __import__(f"repgenr.stages.{stage_name}", fromlist=["run"])
        module.run(ctx, params)
        # Stamp the fingerprint on the record the stage just wrote, so the next
        # invocation with identical params can skip.
        record = ctx.config.stages.get(stage_name)
        if record is not None:
            record.fingerprint = fingerprint
            ctx.save_config()
    except RepGenRError as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=1) from exc


def _parse_key_values(items: list[str], label: str) -> dict[str, str]:
    """Parse repeated ``key=value`` options into a dict (used for tool extras)."""
    out: dict[str, str] = {}
    for item in items:
        if "=" not in item:
            raise UserInputError(f"{label} must be key=value, got '{item}'.")
        key, value = item.split("=", 1)
        key = key.strip()
        if not key:
            raise UserInputError(f"{label} has an empty key in '{item}'.")
        out[key] = value.strip()
    return out


def _read_path_fofn(path: Path) -> list[Path]:
    """Read a file-of-filenames (one path per line; blank lines ignored)."""
    if not path.exists():
        raise UserInputError(f"File not found: {path}")
    return [Path(line.strip()) for line in path.read_text().splitlines() if line.strip()]
