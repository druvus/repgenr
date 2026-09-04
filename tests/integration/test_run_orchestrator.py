"""The `repgenr run` orchestrator chains the canonical stages in order."""

from __future__ import annotations

from typer.testing import CliRunner

from repgenr.cli import cmd_run
from repgenr.cli.main import app

_runner = CliRunner()


def _record(monkeypatch) -> list[str]:
    calls: list[str] = []

    def fake_run(stage, workdir, build, *, create=False):
        build()  # exercise the param builder (catches bad kwargs)
        calls.append(stage)

    monkeypatch.setattr(cmd_run, "_run", fake_run)
    # The wiring tests run without external tools; the preflight has its own test.
    monkeypatch.setattr(cmd_run, "_preflight_tools", lambda *a, **k: None)
    return calls


def test_run_bacterial_chain(monkeypatch, tmp_path) -> None:
    calls = _record(monkeypatch)
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            str(tmp_path),
            "-d",
            "rep",
            "-l",
            "genus",
            "-tg",
            "francisella",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert calls == ["metadata", "genome", "dereplicate", "phylo", "tree2tax"]
    assert "Pipeline complete" in result.stdout


def test_run_viral_chain(monkeypatch, tmp_path) -> None:
    calls = _record(monkeypatch)
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            str(tmp_path),
            "--viral",
            "-t",
            "mastadenovirus",
            "-tg",
            "Mastadenovirus",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert calls == ["vmetadata", "vgenome", "dereplicate", "phylo", "tree2tax"]


def test_run_validates_tool(monkeypatch, tmp_path) -> None:
    _record(monkeypatch)
    result = _runner.invoke(app, ["run", "-wd", str(tmp_path), "--tool", "bogus"])
    assert result.exit_code != 0


def test_run_dry_run_previews_without_executing(monkeypatch, tmp_path) -> None:
    calls = _record(monkeypatch)
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            str(tmp_path),
            "--dry-run",
            "-l",
            "genus",
            "-tg",
            "francisella",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert calls == []  # no stage executed
    assert "[dry-run]" in result.stdout
    assert "dereplicate" in result.stdout and "tree2tax" in result.stdout


def test_run_preflights_every_tool_before_the_first_stage(
    monkeypatch, tmp_path, register_tool
) -> None:
    # A missing tree builder must surface before metadata/genome download, not
    # after dereplication has finished.
    from repgenr.core.errors import MissingBinaryError
    from repgenr.core.plugins import ToolCapabilities
    from repgenr.dereplicators.base import Dereplicator
    from repgenr.dereplicators.base import registry as derep_registry
    from repgenr.treebuilders.base import InputKind, TreeBuilder
    from repgenr.treebuilders.base import registry as tb_registry

    class OkDerep(Dereplicator):
        capabilities = ToolCapabilities(name="okderep")

        def preflight(self):
            return {"okderep": "1.0"}

        def dereplicate(self, genomes, out_dir, params, logger):  # noqa: ANN001
            raise AssertionError("must not run")

    class AbsentBuilder(TreeBuilder):
        capabilities = ToolCapabilities(name="absenttree")
        input_kind = InputKind.GENOMES

        def preflight(self):
            raise MissingBinaryError("absenttree: not found on PATH")

        def build(self, msa_or_genomes, out_dir, params, logger):  # noqa: ANN001
            raise AssertionError("must not run")

    register_tool(derep_registry, "okderep", OkDerep)
    register_tool(tb_registry, "absenttree", AbsentBuilder)
    calls: list[str] = []
    monkeypatch.setattr(cmd_run, "_run", lambda stage, *a, **k: calls.append(stage))
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            str(tmp_path),
            "-l",
            "genus",
            "-tg",
            "francisella",
            "--tool",
            "okderep",
            "--treebuilder",
            "absenttree",
        ],
    )
    assert result.exit_code != 0
    assert calls == []  # nothing downloaded, nothing dereplicated
    assert "absenttree" in result.output
