"""Resume/idempotency: a completed stage with unchanged params is skipped."""

from __future__ import annotations

import sys
import types
from dataclasses import dataclass
from pathlib import Path

from repgenr.cli import base as cli


@dataclass
class _P:
    a: int = 1


def _install_fake_stage(monkeypatch, calls: list[int]) -> None:
    fake = types.ModuleType("repgenr.stages.faketest")

    def run(ctx, params):  # noqa: ANN001
        calls.append(params.a)
        ctx.config.record_stage("faketest", tool="x", params={"a": params.a}, completed="t")
        ctx.save_config()

    fake.run = run  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "repgenr.stages.faketest", fake)


def test_resume_skips_then_force_and_param_change_rerun(tmp_path: Path, monkeypatch) -> None:
    calls: list[int] = []
    _install_fake_stage(monkeypatch, calls)

    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    cli._run("faketest", tmp_path, lambda: _P(), create=True)  # runs
    cli._run("faketest", tmp_path, lambda: _P(), create=True)  # skipped (same params)
    assert calls == [1]

    monkeypatch.setitem(cli._RUN_STATE, "force", True)  # --force re-runs
    cli._run("faketest", tmp_path, lambda: _P(), create=True)
    assert calls == [1, 1]

    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    cli._run("faketest", tmp_path, lambda: _P(a=2), create=True)  # changed params re-runs
    assert calls == [1, 1, 2]


def test_incomplete_stage_reruns(tmp_path: Path, monkeypatch) -> None:
    # A stage that never recorded `completed` (e.g. crashed) must re-run.
    calls: list[int] = []
    fake = types.ModuleType("repgenr.stages.faketest2")

    def run(ctx, params):  # noqa: ANN001
        calls.append(params.a)
        ctx.config.record_stage("faketest2", tool="x", params={})  # no completed stamp
        ctx.save_config()

    fake.run = run  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "repgenr.stages.faketest2", fake)
    monkeypatch.setitem(cli._RUN_STATE, "force", False)

    cli._run("faketest2", tmp_path, lambda: _P(), create=True)
    cli._run("faketest2", tmp_path, lambda: _P(), create=True)
    assert calls == [1, 1]  # not skipped, because completed was never set


def test_registered_input_change_reruns(tmp_path: Path, monkeypatch) -> None:
    """A stage with a STAGE_INPUTS entry reruns when its input dir changes."""
    calls: list[int] = []
    _install_fake_stage(monkeypatch, calls)
    input_dir = tmp_path / "genomes"
    input_dir.mkdir(parents=True)
    (input_dir / "g1.fasta").write_text(">g1\nACGT\n", encoding="utf-8")
    monkeypatch.setitem(cli.STAGE_INPUTS, "faketest", lambda ctx, p: [ctx.genomes_dir])
    monkeypatch.setitem(cli._RUN_STATE, "force", False)

    cli._run("faketest", tmp_path, lambda: _P(), create=True)  # runs
    cli._run("faketest", tmp_path, lambda: _P(), create=True)  # unchanged input -> skip
    assert calls == [1]

    (input_dir / "g2.fasta").write_text(">g2\nACGT\n", encoding="utf-8")
    cli._run("faketest", tmp_path, lambda: _P(), create=True)  # new input file -> rerun
    assert calls == [1, 1]

    cli._run("faketest", tmp_path, lambda: _P(), create=True)  # stable again -> skip
    assert calls == [1, 1]


def test_unregistered_stage_still_skips_on_params(tmp_path: Path, monkeypatch) -> None:
    """A stage without a STAGE_INPUTS entry fingerprints on params alone."""
    calls: list[int] = []
    _install_fake_stage(monkeypatch, calls)
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    assert "faketest" not in cli.STAGE_INPUTS

    cli._run("faketest", tmp_path, lambda: _P(), create=True)
    cli._run("faketest", tmp_path, lambda: _P(), create=True)
    assert calls == [1]
