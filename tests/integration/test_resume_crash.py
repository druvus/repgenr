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
from repgenr.core.errors import (
    MissingBinaryError,
    ToolExecutionError,
    UserInputError,
    WorkdirError,
)


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


# -- first run: a provisional record marks a stage that fails ----------------


def _install_named_fake(monkeypatch, name: str, *, fail: bool, record: bool = True) -> None:
    fake = types.ModuleType(f"repgenr.stages.{name}")

    def run(ctx, params):  # noqa: ANN001
        (ctx.workdir / "partial.txt").write_text("partial\n", encoding="utf-8")
        if fail:
            raise ToolExecutionError(["faketool"], 1, "simulated tool failure")
        if record:
            ctx.config.record_stage(
                name, tool="x", params={"a": params.a}, completed="2026-01-01T00:00:00"
            )
            ctx.save_config()

    fake.run = run  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, f"repgenr.stages.{name}", fake)


def test_first_run_failure_leaves_an_interrupted_record(tmp_path: Path, monkeypatch) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app
    from repgenr.core.doctor import diagnose

    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    _install_named_fake(monkeypatch, "metadata", fail=True)
    with pytest.raises(typer.Exit):
        cli._run("metadata", tmp_path, lambda: _P(), create=True)

    data = yaml.safe_load((tmp_path / CONFIG_FILENAME).read_text(encoding="utf-8"))
    record = data["stages"]["metadata"]
    assert not record.get("completed")
    assert not record.get("fingerprint")
    assert record["params"] == {"a": 1, "query": False}

    result = CliRunner().invoke(app, ["status", "-wd", str(tmp_path)])
    assert result.exit_code == 0
    assert "[interrupted] metadata" in result.stdout
    findings = diagnose(tmp_path)
    assert any(f.level == "fail" and f.area == "metadata" for f in findings)

    # a successful rerun replaces the provisional record with a completed one
    _install_named_fake(monkeypatch, "metadata", fail=False)
    cli._run("metadata", tmp_path, lambda: _P(), create=True)
    data = yaml.safe_load((tmp_path / CONFIG_FILENAME).read_text(encoding="utf-8"))
    assert data["stages"]["metadata"]["completed"]
    result = CliRunner().invoke(app, ["status", "-wd", str(tmp_path)])
    assert "[interrupted]" not in result.stdout


def test_successful_stage_without_its_own_record_leaves_none(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    _install_named_fake(monkeypatch, "crashtest", fail=False, record=False)
    cli._run("crashtest", tmp_path, lambda: _P(), create=True)
    config = tmp_path / CONFIG_FILENAME
    if config.exists():
        data = yaml.safe_load(config.read_text(encoding="utf-8")) or {}
        assert "crashtest" not in (data.get("stages") or {})


def test_parameter_validation_failure_writes_no_record(tmp_path: Path, monkeypatch) -> None:
    calls: list[int] = []
    _install_fake_stage(monkeypatch, calls)
    monkeypatch.setitem(cli._RUN_STATE, "force", False)

    def bad_params() -> _P:
        raise UserInputError("bad flag")

    with pytest.raises(typer.Exit):
        cli._run("crashtest", tmp_path, bad_params, create=True)
    assert calls == []
    config = tmp_path / CONFIG_FILENAME
    if config.exists():
        data = yaml.safe_load(config.read_text(encoding="utf-8")) or {}
        assert "crashtest" not in (data.get("stages") or {})


def test_query_only_first_invocation_writes_no_record(tmp_path: Path, monkeypatch) -> None:
    _install_named_fake(monkeypatch, "crashtest", fail=False, record=False)
    monkeypatch.setitem(cli.QUERY_ONLY_FLAGS, "crashtest", ("query",))
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    cli._run("crashtest", tmp_path, lambda: _P(query=True), create=True)
    assert not (tmp_path / CONFIG_FILENAME).exists()


def test_stuck_interrupted_record_is_cleared_by_a_later_success(
    tmp_path: Path, monkeypatch
) -> None:
    """A stage that failed once and then succeeds without writing a record of
    its own must not leave the earlier [interrupted] record behind."""
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    _install_named_fake(monkeypatch, "crashtest", fail=True)
    monkeypatch.setattr(
        sys.modules["repgenr.stages.crashtest"],
        "run",
        _raising(RuntimeError("tool crashed")),
    )
    with pytest.raises(typer.Exit):
        cli._run("crashtest", tmp_path, lambda: _P(), create=True)
    assert "crashtest" in yaml.safe_load((tmp_path / CONFIG_FILENAME).read_text())["stages"]

    _install_named_fake(monkeypatch, "crashtest", fail=False, record=False)
    cli._run("crashtest", tmp_path, lambda: _P(), create=True)
    data = yaml.safe_load((tmp_path / CONFIG_FILENAME).read_text(encoding="utf-8")) or {}
    assert "crashtest" not in (data.get("stages") or {})


def test_clean_refusal_without_deliverable_change_leaves_no_record(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    monkeypatch.setitem(
        cli.STAGE_DELIVERABLES, "crashtest", lambda ctx, p: [ctx.workdir / "out.txt"]
    )
    fake = types.ModuleType("repgenr.stages.crashtest")
    fake.run = _raising(UserInputError("refused"))  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "repgenr.stages.crashtest", fake)
    with pytest.raises(typer.Exit):
        cli._run("crashtest", tmp_path, lambda: _P(), create=True)
    config = tmp_path / CONFIG_FILENAME
    data = yaml.safe_load(config.read_text(encoding="utf-8")) if config.exists() else {}
    assert "crashtest" not in ((data or {}).get("stages") or {})


def _raising(exc: Exception):
    def run(ctx, params):  # noqa: ANN001
        raise exc

    return run


@pytest.mark.parametrize(
    "exc",
    [UserInputError("refused"), WorkdirError("too few"), MissingBinaryError("tool missing")],
)
def test_refused_rerun_keeps_the_finished_record(tmp_path: Path, monkeypatch, exc) -> None:
    # A re-run refused before it touched any deliverable (a tool missing at
    # preflight, exit 4; too few genomes, exit 3) used to leave the finished
    # record as interrupted although its outputs were intact.
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    monkeypatch.setitem(
        cli.STAGE_DELIVERABLES, "crashtest", lambda ctx, p: [ctx.workdir / "out.txt"]
    )
    calls: list[int] = []
    _install_fake_stage(monkeypatch, calls)
    cli._run("crashtest", tmp_path, lambda: _P(), create=True)
    (tmp_path / "out.txt").write_text("result\n", encoding="utf-8")
    finished = _record_from_disk(tmp_path)
    assert finished["completed"]

    fake = types.ModuleType("repgenr.stages.crashtest")
    fake.run = _raising(exc)  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "repgenr.stages.crashtest", fake)
    with pytest.raises(typer.Exit) as raised:
        cli._run("crashtest", tmp_path, lambda: _P(a=2), create=True)
    assert raised.value.exit_code == exc.exit_code
    assert _record_from_disk(tmp_path) == finished


def test_refusal_that_changed_a_deliverable_stays_interrupted(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    monkeypatch.setitem(
        cli.STAGE_DELIVERABLES, "crashtest", lambda ctx, p: [ctx.workdir / "out.txt"]
    )
    calls: list[int] = []
    _install_fake_stage(monkeypatch, calls)
    cli._run("crashtest", tmp_path, lambda: _P(), create=True)
    (tmp_path / "out.txt").write_text("result\n", encoding="utf-8")

    def run(ctx, params):  # noqa: ANN001
        (ctx.workdir / "out.txt").write_text("partial, longer\n", encoding="utf-8")
        raise WorkdirError("refused after writing")

    fake = types.ModuleType("repgenr.stages.crashtest")
    fake.run = run  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "repgenr.stages.crashtest", fake)
    with pytest.raises(typer.Exit):
        cli._run("crashtest", tmp_path, lambda: _P(a=2), create=True)
    assert _record_from_disk(tmp_path)["completed"] is None


def test_tool_missing_on_a_first_run_leaves_no_record(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    monkeypatch.setitem(
        cli.STAGE_DELIVERABLES, "crashtest", lambda ctx, p: [ctx.workdir / "out.txt"]
    )
    fake = types.ModuleType("repgenr.stages.crashtest")
    fake.run = _raising(MissingBinaryError("gubbins: not found on PATH"))  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "repgenr.stages.crashtest", fake)
    with pytest.raises(typer.Exit) as raised:
        cli._run("crashtest", tmp_path, lambda: _P(), create=True)
    assert raised.value.exit_code == 4
    data = yaml.safe_load((tmp_path / CONFIG_FILENAME).read_text(encoding="utf-8")) or {}
    assert "crashtest" not in (data.get("stages") or {})
