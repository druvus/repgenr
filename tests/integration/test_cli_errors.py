"""CLI error handling: clean exits instead of raw tracebacks."""

from __future__ import annotations

import logging
import sys
import types
from pathlib import Path

import pytest
import typer

from repgenr.cli import base as cli
from repgenr.core.errors import (
    MissingBinaryError,
    PluginError,
    ToolExecutionError,
    UserInputError,
    WorkdirError,
)
from repgenr.core.logging import configure_logging


@pytest.mark.parametrize(
    ("exc", "code"),
    [
        (UserInputError("bad flag"), 2),
        (WorkdirError("missing genomes"), 3),
        (MissingBinaryError("skder: not found"), 4),
        (PluginError("no such tool"), 5),
        (ToolExecutionError(["skder"], 137), 6),
    ],
    ids=["input", "workdir", "binary", "plugin", "tool"],
)
def test_stage_errors_exit_code_names_the_error_class(exc, code, monkeypatch) -> None:
    """A caller can tell a usage error from a missing tool or a failed run."""
    monkeypatch.delenv("REPGENR_PROPAGATE_TOOL_EXIT", raising=False)
    logger = logging.getLogger("repgenr")
    with pytest.raises(typer.Exit) as ei:  # noqa: PT012
        with cli.stage_errors(logger):
            raise exc
    assert ei.value.exit_code == code


def test_tool_exit_is_forwarded_when_propagation_is_requested(monkeypatch) -> None:
    monkeypatch.setenv("REPGENR_PROPAGATE_TOOL_EXIT", "1")
    logger = logging.getLogger("repgenr")
    with pytest.raises(typer.Exit) as ei:  # noqa: PT012
        with cli.stage_errors(logger):
            raise ToolExecutionError(["skder"], -9)  # killed by SIGKILL
    assert ei.value.exit_code == 137


def test_stage_errors_unexpected_exits_and_logs_traceback(workdir: Path) -> None:
    # With a workdir log present, the concise message is on console and the full
    # traceback is captured in repgenr.log.
    workdir.mkdir(parents=True)
    logger = configure_logging(workdir, level=logging.INFO)
    with pytest.raises(typer.Exit) as ei:  # noqa: PT012
        with cli.stage_errors(logger):
            raise ValueError("boom")
    assert ei.value.exit_code == 1
    for h in logger.handlers:  # flush file handler
        h.flush()
    log_text = (workdir / "repgenr.log").read_text()
    assert "Unexpected error: boom" in log_text
    assert "Traceback" in log_text  # full traceback captured to the file


def _install_fake_stage(monkeypatch, exc: Exception) -> None:
    fake = types.ModuleType("repgenr.stages.faketest")

    def run(ctx, params):  # noqa: ANN001
        raise exc

    fake.run = run  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "repgenr.stages.faketest", fake)


def test_run_unexpected_exception_is_clean_exit(tmp_path: Path, monkeypatch) -> None:
    _install_fake_stage(monkeypatch, KeyError("changed_api_field"))
    monkeypatch.setitem(cli._RUN_STATE, "force", False)

    class _P:
        a = 1

    with pytest.raises(typer.Exit) as ei:
        cli._run("faketest", tmp_path, lambda: _P(), create=True)
    assert ei.value.exit_code == 1
    # the workdir log holds the traceback for diagnosis
    assert "Traceback" in (tmp_path / "repgenr.log").read_text()
