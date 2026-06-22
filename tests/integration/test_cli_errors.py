"""CLI error handling: clean exits instead of raw tracebacks."""

from __future__ import annotations

import logging
import sys
import types
from pathlib import Path

import pytest
import typer

from repgenr.cli import base as cli
from repgenr.core.errors import WorkdirError
from repgenr.core.logging import configure_logging


def test_stage_errors_repgenr_error_exits(caplog) -> None:
    logger = logging.getLogger("repgenr")
    with pytest.raises(typer.Exit) as ei:  # noqa: PT012
        with cli.stage_errors(logger):
            raise WorkdirError("missing genomes")
    assert ei.value.exit_code == 1


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
