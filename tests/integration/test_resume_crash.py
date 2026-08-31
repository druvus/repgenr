"""Crash-and-restart correctness of the resume harness.

A stage that previously completed and then crashes while re-running must not be
skipped on the next invocation: `_run` dirties the prior record (completed and
fingerprint cleared) before the stage body executes. Query-only invocations
(--accession-list-only / --list / --glance) must not touch the record at all.
"""

from __future__ import annotations

import sys
import types
from dataclasses import dataclass
from pathlib import Path

import pytest
import typer
import yaml

from repgenr.cli import base as cli
from repgenr.core.config import CONFIG_FILENAME


@dataclass
class _P:
    a: int = 1
    query: bool = False


def _install_fake_stage(monkeypatch, calls: list[int], *, crash_on: set[int] | None = None):
    fake = types.ModuleType("repgenr.stages.crashtest")
    crash_on = crash_on or set()

    def run(ctx, params):  # noqa: ANN001
        call_index = len(calls)
        calls.append(params.a)
        if call_index in crash_on:
            raise RuntimeError("simulated crash mid-stage")
        ctx.config.record_stage(
            "crashtest", tool="x", params={"a": params.a}, completed="2026-01-01T00:00:00"
        )
        ctx.save_config()

    fake.run = run  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "repgenr.stages.crashtest", fake)


def _record_from_disk(workdir: Path) -> dict:
    data = yaml.safe_load((workdir / CONFIG_FILENAME).read_text(encoding="utf-8"))
    return data["stages"]["crashtest"]


def test_crash_on_forced_rerun_is_not_skipped_afterwards(tmp_path: Path, monkeypatch) -> None:
    calls: list[int] = []
    _install_fake_stage(monkeypatch, calls, crash_on={1})
    monkeypatch.setitem(cli._RUN_STATE, "force", False)

    cli._run("crashtest", tmp_path, lambda: _P(), create=True)  # completes
    assert calls == [1]

    monkeypatch.setitem(cli._RUN_STATE, "force", True)
    with pytest.raises(typer.Exit):
        cli._run("crashtest", tmp_path, lambda: _P(), create=True)  # crashes mid-stage
    assert calls == [1, 1]

    # the crashed stage's record is dirty on disk: it must not look completed
    record = _record_from_disk(tmp_path)
    assert not record.get("completed")

    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    cli._run("crashtest", tmp_path, lambda: _P(), create=True)  # must re-run, not skip
    assert calls == [1, 1, 1]


def test_crash_on_input_change_rerun_is_not_skipped_afterwards(tmp_path: Path, monkeypatch) -> None:
    calls: list[int] = []
    _install_fake_stage(monkeypatch, calls, crash_on={1})
    input_dir = tmp_path / "genomes"
    input_dir.mkdir(parents=True)
    (input_dir / "g1.fasta").write_text(">g1\nACGT\n", encoding="utf-8")
    monkeypatch.setitem(cli.STAGE_INPUTS, "crashtest", lambda ctx, p: [ctx.genomes_dir])
    monkeypatch.setitem(cli._RUN_STATE, "force", False)

    cli._run("crashtest", tmp_path, lambda: _P(), create=True)  # completes
    (input_dir / "g2.fasta").write_text(">g2\nACGT\n", encoding="utf-8")
    with pytest.raises(typer.Exit):
        cli._run("crashtest", tmp_path, lambda: _P(), create=True)  # rerun on input change, crashes
    assert calls == [1, 1]

    cli._run("crashtest", tmp_path, lambda: _P(), create=True)  # must re-run
    assert calls == [1, 1, 1]


def test_skipped_invocation_does_not_dirty_the_record(tmp_path: Path, monkeypatch) -> None:
    calls: list[int] = []
    _install_fake_stage(monkeypatch, calls)
    monkeypatch.setitem(cli._RUN_STATE, "force", False)

    cli._run("crashtest", tmp_path, lambda: _P(), create=True)  # completes
    cli._run("crashtest", tmp_path, lambda: _P(), create=True)  # skips
    assert calls == [1]
    record = _record_from_disk(tmp_path)
    assert record.get("completed")
    assert record.get("fingerprint")


def test_query_only_invocation_leaves_record_untouched(tmp_path: Path, monkeypatch) -> None:
    calls: list[int] = []
    _install_fake_stage(monkeypatch, calls)
    monkeypatch.setitem(cli.QUERY_ONLY_FLAGS, "crashtest", ("query",))
    monkeypatch.setitem(cli._RUN_STATE, "force", False)

    cli._run("crashtest", tmp_path, lambda: _P(), create=True)  # real run completes
    before = _record_from_disk(tmp_path)

    # A query-only invocation runs the stage body but records nothing --
    # the fake's record_stage call is what a real query path would skip, so
    # install a variant that early-returns like the real query paths do.
    fake = types.ModuleType("repgenr.stages.crashtest")

    def run(ctx, params):  # noqa: ANN001
        calls.append(params.a)  # no record_stage, like --accession-list-only

    fake.run = run  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "repgenr.stages.crashtest", fake)

    cli._run("crashtest", tmp_path, lambda: _P(a=7, query=True), create=True)
    assert calls == [1, 7]  # ran despite the completed record (no skip check)
    after = _record_from_disk(tmp_path)
    assert after == before  # record untouched: not dirtied, not restamped

    # a second query invocation is NOT skipped either
    cli._run("crashtest", tmp_path, lambda: _P(a=8, query=True), create=True)
    assert calls == [1, 7, 8]


def test_query_only_flag_false_behaves_normally(tmp_path: Path, monkeypatch) -> None:
    calls: list[int] = []
    _install_fake_stage(monkeypatch, calls)
    monkeypatch.setitem(cli.QUERY_ONLY_FLAGS, "crashtest", ("query",))
    monkeypatch.setitem(cli._RUN_STATE, "force", False)

    cli._run("crashtest", tmp_path, lambda: _P(query=False), create=True)
    cli._run("crashtest", tmp_path, lambda: _P(query=False), create=True)  # skips normally
    assert calls == [1]
