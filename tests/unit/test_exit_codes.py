"""Exit-code mapping at the CLI boundary.

By default every failure exits 1. Under REPGENR_PROPAGATE_TOOL_EXIT=1 (set by
the Nextflow modules), a failed external tool's exit code is forwarded -- with
signal kills mapped to 128+signum -- so Nextflow's retry-on-exitStatus rule can
distinguish an OOM kill (137) from an ordinary error.
"""

from __future__ import annotations

import logging

import pytest
import typer

from repgenr.cli.base import stage_errors
from repgenr.core.errors import ToolExecutionError, UserInputError

_LOG = logging.getLogger("test")


def _exit_code(exc: BaseException) -> int:
    with pytest.raises(typer.Exit) as ei:
        with stage_errors(_LOG):
            raise exc
    return ei.value.exit_code


def test_default_collapses_to_one(monkeypatch) -> None:
    monkeypatch.delenv("REPGENR_PROPAGATE_TOOL_EXIT", raising=False)
    assert _exit_code(ToolExecutionError(["tool"], 137, output="oom")) == 1


def test_propagates_tool_exit_code(monkeypatch) -> None:
    monkeypatch.setenv("REPGENR_PROPAGATE_TOOL_EXIT", "1")
    assert _exit_code(ToolExecutionError(["tool"], 137, output="oom")) == 137
    assert _exit_code(ToolExecutionError(["tool"], 2, output="bad args")) == 2


def test_signal_kill_maps_to_128_plus_signum(monkeypatch) -> None:
    monkeypatch.setenv("REPGENR_PROPAGATE_TOOL_EXIT", "1")
    assert _exit_code(ToolExecutionError(["tool"], -9, output="killed")) == 137


def test_out_of_range_code_falls_back_to_one(monkeypatch) -> None:
    monkeypatch.setenv("REPGENR_PROPAGATE_TOOL_EXIT", "1")
    assert _exit_code(ToolExecutionError(["tool"], 300, output="odd")) == 1
    assert _exit_code(ToolExecutionError(["tool"], 0, output="odd")) == 1


def test_other_errors_still_exit_one(monkeypatch) -> None:
    monkeypatch.setenv("REPGENR_PROPAGATE_TOOL_EXIT", "1")
    assert _exit_code(UserInputError("bad flag")) == 1
    assert _exit_code(RuntimeError("boom")) == 1
