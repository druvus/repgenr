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
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import typer

from .. import __version__
from ..core.context import WorkdirContext
from ..core.contracts import CLUSTERS_TSV, SELECTION_TSV, TREE_NWK
from ..core.errors import RepGenRError, ToolExecutionError, UserInputError
from ..core.inputs import inputs_digest, manifest_digest_for_stage
from ..core.logging import configure_logging

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="RepGenR: modular genome dereplication, alignment, SNP typing and phylogenetics.",
)

# Top-level run options shared by every subcommand (set in the callback).
_RUN_STATE: dict[str, Any] = {"force": False, "log_level": logging.INFO}

# One default thread count for every stage's -t/--threads, so the CLI is
# consistent (stages previously mixed 16 and 24).
DEFAULT_THREADS = 16

# Canonical stage order per lineage. Used to show progress (`status`) and by
# `run --dry-run` to print the chain.
PIPELINE_BACTERIAL = ("metadata", "genome", "dereplicate", "phylo", "tree2tax")
PIPELINE_VIRAL = ("vmetadata", "vgenome", "dereplicate", "phylo", "tree2tax")


def _phylo_inputs(ctx: WorkdirContext, params: Any) -> list[Path]:
    # snp/core_snp.fasta is deliberately NOT declared for msa_source=snptype:
    # phylo regenerates it from the same genome set, so declaring it would make
    # the stage's fingerprint depend on its own output and force a spurious
    # rerun on every second invocation.
    return [
        ctx.genomes_dir if getattr(params, "all_genomes", False)
        else ctx.representatives_dir,
        ctx.outgroup_dir,
    ]


def _tree2tax_inputs(ctx: WorkdirContext, params: Any) -> list[Path]:
    paths = [ctx.tree_dir / TREE_NWK]
    if getattr(params, "include_dereplicated", False):
        paths.append(ctx.derep_dir / CLUSTERS_TSV)
    return paths


# What each stage reads, for the input digests in the resume fingerprint:
# stage -> callable(ctx, params) -> paths (directories are digested from file
# metadata, files by content; see core.inputs). Conditional edges (--all-genomes,
# --msa-source snptype, --include-dereplicated) live in the helpers above.
# Stages not listed digest no inputs and fingerprint on params alone.
STAGE_INPUTS: dict[str, Any] = {
    "metadata": lambda ctx, p: [],  # network-only
    "vmetadata": lambda ctx, p: [],
    "genome": lambda ctx, p: [ctx.workdir / SELECTION_TSV],
    # vgenome WRITES selection.tsv, so its inputs are the vmetadata download
    # artifacts (records path and legacy BV-BRC tables; absent ones digest to
    # the stable sentinel).
    "vgenome": lambda ctx, p: [
        ctx.workdir / "virus_download_wd" / "download.fa",
        ctx.workdir / "virus_download_wd" / "virus_records.json",
        ctx.workdir / "virus_download_wd" / "metadata_base.tsv",
        ctx.workdir / "virus_download_wd" / "metadata_ncbi.tsv",
    ],
    "dereplicate": lambda ctx, p: [ctx.genomes_dir],
    "snptype": lambda ctx, p: [
        ctx.genomes_dir if getattr(p, "all_genomes", False) else ctx.representatives_dir
    ],
    "phylo": _phylo_inputs,
    "tree2tax": _tree2tax_inputs,
}

# Stages whose result also depends on the manifest's genome rows (taxonomy,
# derep status, CheckM quality), digested from ordered query results.
# "dereplicate" reads the manifest for the quality-aware keeper and --reduce
# taxonomy grouping, so a manifest-only edit (no genome file touched) must
# still invalidate a prior resume.
_MANIFEST_INPUT_STAGES = frozenset({"tree2tax", "dereplicate"})

# Param flags that turn a stage invocation into a pure query (list/preview
# modes that write no pipeline outputs). Such invocations bypass the resume
# machinery entirely: no skip check, no dirty marker, no fingerprint stamp --
# otherwise a query would either be wrongly skipped or would restamp/dirty the
# record of the last real run.
QUERY_ONLY_FLAGS: dict[str, tuple[str, ...]] = {
    "genome": ("accession_list_only",),
    "vmetadata": ("list_targets",),
    "vgenome": ("glance",),
}


def _is_query_only(stage_name: str, params: Any) -> bool:
    return any(
        getattr(params, flag, False) for flag in QUERY_ONLY_FLAGS.get(stage_name, ())
    )


def _stage_input_digests(ctx: WorkdirContext, stage_name: str, params: Any) -> dict[str, str]:
    """Digest a stage's declared inputs; empty for stages with no declaration."""
    spec = STAGE_INPUTS.get(stage_name)
    if spec is None:
        return {}
    digests = inputs_digest(ctx.workdir, spec(ctx, params))
    if stage_name in _MANIFEST_INPUT_STAGES:
        digests["manifest"] = manifest_digest_for_stage(stage_name, ctx.manifest)
    return digests


def _env_fragment() -> dict[str, Any]:
    """Result-affecting execution environment for the fingerprint.

    The container backend/platform/wave selection changes which tool builds run
    (and so can change results); the engine binary, cache directory, and extra
    mounts are plumbing and deliberately excluded, like _NON_RESULT_PARAMS.
    """
    from ..core.containers import get_config

    config = get_config()
    return {"container": [config.backend, config.platform, config.wave_enabled]}


def _derep_help(*, auto: bool = True) -> str:
    from ..core.plugins import tool_choices_help
    from ..dereplicators.base import registry

    return tool_choices_help(registry, auto=auto)


def _tree_help(*, auto: bool = True) -> str:
    from ..core.plugins import tool_choices_help
    from ..treebuilders.base import registry

    return tool_choices_help(registry, auto=auto)


def _aligner_help() -> str:
    from ..aligners.base import registry
    from ..core.plugins import tool_choices_help

    return tool_choices_help(registry, auto=False)


def _snp_help() -> str:
    from ..core.plugins import tool_choices_help
    from ..snptypers.base import registry

    return tool_choices_help(registry, auto=False, prefix="SNP typer: ")


def _mask_help() -> str:
    from ..core.plugins import tool_choices_help
    from ..maskers.base import registry

    return tool_choices_help(registry, auto=False, prefix="Recombination masking: none, ")


def _require_choice(value: str, choices: set[str], label: str) -> None:
    if value not in choices:
        raise UserInputError(
            f"Invalid {label} {value!r}. Choose from: {', '.join(sorted(choices))}."
        )


def _require_unit_interval(value: float | None, label: str) -> None:
    if value is not None and not (0.0 < value <= 1.0):
        raise UserInputError(f"{label} must be in (0, 1], got {value}.")


# Parameters that change how work is scheduled but not the result, so they are
# excluded from the resume fingerprint -- changing --threads / --num-processes
# must not force an otherwise-identical stage to recompute from scratch.
# allow_incomplete only gates the input-completeness refusal; on complete
# inputs it changes nothing, so it must not invalidate the resume cache.
_NON_RESULT_PARAMS = frozenset({"threads", "num_processes", "allow_incomplete"})


# Fingerprint format version. Bumping it guarantees fingerprints from older
# releases never false-match, so old workdirs rerun once instead of skipping
# against semantics they were not computed under.
_FINGERPRINT_VERSION = 2


def _stage_fingerprint(
    stage_name: str, params: object, inputs: dict[str, str], env: dict[str, Any]
) -> str:
    """Stable hash of a stage invocation, used to skip already-completed work.

    Built from the stage name, the parameter object (a dataclass), the digests
    of the stage's declared inputs, and the result-affecting environment
    (container identity), so a skip means "same request, same inputs, same
    execution environment". Paths and other non-JSON values are stringified.
    Non-result params (thread/worker counts) are excluded so they do not
    invalidate the resume cache.
    """
    if dataclasses.is_dataclass(params) and not isinstance(params, type):
        payload: dict = dataclasses.asdict(params)
    else:
        payload = dict(vars(params))
    payload = {k: v for k, v in payload.items() if k not in _NON_RESULT_PARAMS}
    blob = json.dumps(
        {
            "fpv": _FINGERPRINT_VERSION,
            "stage": stage_name,
            "params": payload,
            "inputs": inputs,
            "env": env,
        },
        sort_keys=True, default=str,
    )
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


def _tool_exit_code(returncode: int) -> int:
    """Exit code for a failed external tool.

    Interactive use keeps the historical exit 1. Under
    ``REPGENR_PROPAGATE_TOOL_EXIT=1`` (set by the Nextflow modules) the tool's
    code is forwarded, with a signal kill mapped to 128+signum (SIGKILL -> 137),
    so Nextflow's retry-on-exitStatus rule can react to e.g. an OOM kill.
    """
    if os.environ.get("REPGENR_PROPAGATE_TOOL_EXIT", "") in ("", "0"):
        return 1
    if returncode < 0:
        code = 128 - returncode
    else:
        code = returncode
    return code if 0 < code <= 255 else 1


@contextmanager
def stage_errors(logger: logging.Logger) -> Iterator[None]:
    """Turn errors into clean CLI exits instead of raw tracebacks.

    A :class:`RepGenRError` (expected, user-facing) is logged concisely. Any
    other exception is unexpected: a concise message goes to the console and the
    full traceback is captured -- in the run log when one exists (DEBUG), or on
    the console otherwise (e.g. a data-channel step with no workdir). Both exit
    non-zero. ``--verbose`` shows the traceback on the console too.
    """
    try:
        yield
    except typer.Exit:
        raise
    except ToolExecutionError as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=_tool_exit_code(exc.returncode)) from exc
    except RepGenRError as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=1) from exc
    except Exception as exc:
        logger.error("Unexpected error: %s", exc)
        if any(isinstance(h, logging.FileHandler) for h in logger.handlers):
            logger.debug("Full traceback:", exc_info=True)
            logger.error("See the run log for the full traceback, or re-run with --verbose.")
        else:
            # No persistent log (data-channel step): surface the traceback now.
            logger.error("Full traceback:", exc_info=True)
        raise typer.Exit(code=1) from exc


def _run(stage_name: str, workdir: Path, build_params, *, create: bool = False) -> None:
    """Common harness: context, dispatch, clean error handling.

    Resume: a stage that already completed with the same parameters, the same
    inputs (per STAGE_INPUTS digests), and the same container identity is
    skipped, unless ``--force`` is set. Re-running an upstream stage changes a
    downstream stage's input digests, so the downstream stage reruns
    automatically. A stage that crashed before recording completion has no
    ``completed`` stamp and so always re-runs.
    """
    logger = configure_logging(
        workdir if (create or workdir.exists()) else None, level=_RUN_STATE["log_level"]
    )
    with stage_errors(logger):
        ctx = WorkdirContext(workdir, logger=logger, create=create)
        try:
            _run_stage(stage_name, ctx, build_params, logger)
        finally:
            ctx.close()


def _run_stage(stage_name: str, ctx: WorkdirContext, build_params, logger) -> None:
    params = build_params()
    if _is_query_only(stage_name, params):
        # Pure query (list/preview): run the body, leave the resume record
        # of the last real run untouched.
        module = __import__(f"repgenr.stages.{stage_name}", fromlist=["run"])
        module.run(ctx, params)
        return
    # Digested once: upstream inputs are stable while this stage executes,
    # so the same digests are stamped onto the record after the run.
    digests = _stage_input_digests(ctx, stage_name, params)
    fingerprint = _stage_fingerprint(stage_name, params, digests, _env_fragment())
    prior = ctx.config.stages.get(stage_name)
    if not _RUN_STATE["force"] and prior is not None and prior.completed:
        if prior.fingerprint == fingerprint:
            logger.info(
                "Stage '%s' already completed with the same parameters and "
                "inputs; skipping (use --force to re-run).", stage_name,
            )
            return
        changed = sorted(
            key for key in {*prior.inputs, *digests}
            if prior.inputs.get(key) != digests.get(key)
        )
        if changed and prior.inputs:
            logger.info(
                "Stage '%s': input %s changed since last completion; re-running.",
                stage_name, ", ".join(f"'{c}'" for c in changed),
            )
    if prior is not None and prior.completed:
        # Dirty the record before the stage body runs: a crash mid-stage
        # must not leave a completed-looking record over partial outputs.
        prior.completed = None
        prior.fingerprint = None
        ctx.save_config()
    module = __import__(f"repgenr.stages.{stage_name}", fromlist=["run"])
    module.run(ctx, params)
    # Stamp fingerprint + input digests on the record the stage just wrote,
    # so the next invocation can skip.
    record = ctx.config.stages.get(stage_name)
    if record is not None:
        record.fingerprint = fingerprint
        record.inputs = digests
        ctx.save_config()


def gated_extra(registry, tool: str, key: str, value: object) -> dict:
    """Return ``{key: value}`` only when ``tool`` reads that extra.

    Injecting a key a tool ignores would change the resume fingerprint without
    changing the result. ``auto`` passes the key through; the stage warns after
    it has picked a concrete tool.
    """
    if tool != "auto":
        caps = registry.get(tool).capabilities
        if key not in caps.accepted_extras:
            return {}
    return {key: value}


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
    lines = path.read_text(encoding="utf-8").splitlines()
    return [Path(line.strip()) for line in lines if line.strip()]
