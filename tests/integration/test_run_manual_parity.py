"""`repgenr run` and the manual per-stage commands fingerprint identically.

Both entry points are driven through CliRunner with the `_run` harness replaced
by a recorder; for equivalent CLI intent, every stage's params object -- and
therefore its resume fingerprint -- must be equal, so a `run` followed by a
manual command (or vice versa) skips instead of silently re-running.
"""

from __future__ import annotations

import pytest
from typer.testing import CliRunner

from repgenr.cli import cmd_bacterial, cmd_phylo, cmd_run, cmd_viral
from repgenr.cli.base import _stage_fingerprint
from repgenr.cli.main import app

_runner = CliRunner()


@pytest.fixture()
def dispatched(monkeypatch) -> list[tuple]:
    calls: list[tuple] = []

    def fake_run(stage, workdir, build, *, create=False):
        calls.append((stage, build()))

    for mod in (cmd_run, cmd_bacterial, cmd_phylo, cmd_viral):
        monkeypatch.setattr(mod, "_run", fake_run)
    monkeypatch.setattr(cmd_run, "_preflight_tools", lambda *a, **k: None)
    return calls


def _fingerprints(calls: list[tuple]) -> dict[str, str]:
    return {stage: _stage_fingerprint(stage, params, {}, {}) for stage, params in calls}


def test_bacterial_run_matches_manual_commands(dispatched, tmp_path) -> None:
    wd = str(tmp_path)
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            wd,
            "-d",
            "rep",
            "-l",
            "genus",
            "-tg",
            "francisella",
            "-r",
            "232.0",
            "--gtdb-version",
            "bac120",
            "--tool",
            "skder",
            "--treebuilder",
            "mashtree",
            "--keeper",
            "tool",
        ],
    )
    assert result.exit_code == 0, result.output
    run_fps = _fingerprints(dispatched)
    dispatched.clear()

    for args in (
        [
            "metadata",
            "-wd",
            wd,
            "-d",
            "rep",
            "-l",
            "genus",
            "-tg",
            "francisella",
            "-r",
            "232.0",
            "--gtdb-version",
            "bac120",
        ],
        ["genome", "-wd", wd],
        ["dereplicate", "-wd", wd, "--tool", "skder", "--keeper", "tool"],
        ["phylo", "-wd", wd, "--treebuilder", "mashtree"],
        ["tree2tax", "-wd", wd],
    ):
        result = _runner.invoke(app, args)
        assert result.exit_code == 0, result.output
    manual_fps = _fingerprints(dispatched)

    assert set(run_fps) == {"metadata", "genome", "dereplicate", "phylo", "tree2tax"}
    for stage, fp in run_fps.items():
        assert manual_fps[stage] == fp, f"run vs manual fingerprint diverges for {stage}"


def test_viral_run_matches_manual_commands(dispatched, tmp_path) -> None:
    wd = str(tmp_path)
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            wd,
            "--viral",
            "-t",
            "adenoviridae",
            "-tg",
            "mastadenovirus",
            "--tool",
            "skder",
            "--treebuilder",
            "mashtree",
        ],
    )
    assert result.exit_code == 0, result.output
    run_fps = _fingerprints(dispatched)
    dispatched.clear()

    for args in (
        ["vmetadata", "-wd", wd, "-t", "adenoviridae"],
        ["vgenome", "-wd", wd, "-tg", "mastadenovirus"],
        # skder ignores extra["virus"], so `run --viral` no longer injects it;
        # the equivalent manual invocation omits --virus.
        ["dereplicate", "-wd", wd, "--tool", "skder"],
        ["phylo", "-wd", wd, "--treebuilder", "mashtree"],
        ["tree2tax", "-wd", wd],
    ):
        result = _runner.invoke(app, args)
        assert result.exit_code == 0, result.output
    manual_fps = _fingerprints(dispatched)

    assert set(run_fps) == {"vmetadata", "vgenome", "dereplicate", "phylo", "tree2tax"}
    for stage, fp in run_fps.items():
        assert manual_fps[stage] == fp, f"run vs manual fingerprint diverges for {stage}"


def test_run_include_dereplicated_flag(dispatched, tmp_path) -> None:
    wd = str(tmp_path)
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            wd,
            "-l",
            "genus",
            "-tg",
            "x",
            "-r",
            "232.0",
            "--gtdb-version",
            "bac120",
            "--no-include-dereplicated",
        ],
    )
    assert result.exit_code == 0, result.output
    params = dict(dispatched)["tree2tax"]
    assert params.include_dereplicated is False


def test_tree2tax_defaults_to_include_dereplicated(dispatched, tmp_path) -> None:
    result = _runner.invoke(app, ["tree2tax", "-wd", str(tmp_path)])
    assert result.exit_code == 0, result.output
    params = dict(dispatched)["tree2tax"]
    assert params.include_dereplicated is True


def test_bacterial_run_requires_level(dispatched, tmp_path) -> None:
    result = _runner.invoke(app, ["run", "-wd", str(tmp_path), "-tg", "francisella"])
    assert result.exit_code != 0
    assert dispatched == []


def test_run_snptype_msa_source_matches_manual_phylo(dispatched, tmp_path) -> None:
    wd = str(tmp_path)
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            wd,
            "-l",
            "genus",
            "-tg",
            "x",
            "-r",
            "232.0",
            "--gtdb-version",
            "bac120",
            "--msa-source",
            "snptype",
            "--snptyper",
            "simple",
            "--treebuilder",
            "iqtree",
        ],
    )
    assert result.exit_code == 0, result.output
    run_phylo = dict(dispatched)["phylo"]
    assert (run_phylo.msa_source, run_phylo.snptyper) == ("snptype", "simple")
    run_fp = _stage_fingerprint("phylo", run_phylo, {}, {})
    dispatched.clear()

    result = _runner.invoke(
        app,
        [
            "phylo",
            "-wd",
            wd,
            "--msa-source",
            "snptype",
            "--snptyper",
            "simple",
            "--treebuilder",
            "iqtree",
        ],
    )
    assert result.exit_code == 0, result.output
    manual_fp = _stage_fingerprint("phylo", dict(dispatched)["phylo"], {}, {})
    assert manual_fp == run_fp


def test_run_keeper_flag_reaches_dereplicate_params(dispatched, tmp_path) -> None:
    wd = str(tmp_path)
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            wd,
            "-l",
            "genus",
            "-tg",
            "x",
            "-r",
            "232.0",
            "--gtdb-version",
            "bac120",
            "--keeper",
            "tool",
        ],
    )
    assert result.exit_code == 0, result.output
    params = dict(dispatched)["dereplicate"]
    assert params.keeper == "tool"


def test_run_keeper_defaults_to_quality(dispatched, tmp_path) -> None:
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            str(tmp_path),
            "-l",
            "genus",
            "-tg",
            "x",
            "-r",
            "232.0",
            "--gtdb-version",
            "bac120",
        ],
    )
    assert result.exit_code == 0, result.output
    params = dict(dispatched)["dereplicate"]
    assert params.keeper == "quality"


def test_run_bad_msa_source_fails(dispatched, tmp_path) -> None:
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            str(tmp_path),
            "-l",
            "genus",
            "-tg",
            "x",
            "-r",
            "232.0",
            "--gtdb-version",
            "bac120",
            "--msa-source",
            "nosuch",
        ],
    )
    assert result.exit_code != 0
    assert dispatched == []
