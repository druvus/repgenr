"""--bin-dir TOOL=DIR puts one directory first on PATH for one tool."""

from __future__ import annotations

import logging
import os
import stat
from pathlib import Path

import pytest
from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core import bindirs
from repgenr.core.binaries import BinarySpec
from repgenr.core.config import StageRecord
from repgenr.core.containers import configure_container, run_tool
from repgenr.core.errors import MissingBinaryError, UserInputError
from repgenr.core.plugins import ToolCapabilities, preflight, tool_available
from repgenr.maskers.gubbins import resolve_tree_builder

_runner = CliRunner()
_LOG = logging.getLogger("test-bin-dirs")


@pytest.fixture(autouse=True)
def _reset():
    yield
    bindirs.configure_bin_dirs({})
    configure_container("none")


def _script(directory: Path, name: str, body: str) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(f"#!/bin/sh\n{body}\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


def _caps(name: str, binary: str, floor: str | None = None) -> ToolCapabilities:
    return ToolCapabilities(name=name, required_binaries=(BinarySpec(binary, min_version=floor),))


def test_a_tool_is_found_only_through_its_directory(tmp_path) -> None:
    sat = tmp_path / "sat" / "bin"
    _script(sat, "fakegubbins", 'echo "fakegubbins 3.4.1"')
    caps = _caps("fakegub", "fakegubbins", "3.0")
    with pytest.raises(MissingBinaryError):
        preflight(caps)
    assert not tool_available(caps)

    bindirs.configure_bin_dirs({"fakegub": sat})
    assert preflight(caps) == {"fakegubbins": "3.4.1"}
    assert tool_available(caps)
    # Another tool's lookup does not see the directory.
    assert not tool_available(_caps("other", "fakegubbins"))


def test_the_version_query_runs_with_the_directory_first(tmp_path) -> None:
    # A wrapper that calls a helper by name (as run_gubbins.py calls python)
    # must get the helper from its own directory.
    sat = tmp_path / "sat"
    _script(sat, "helper", 'echo "helper 9.1"')
    _script(sat, "wrapped", "helper")
    bindirs.configure_bin_dirs({"w": sat})
    assert preflight(_caps("w", "wrapped", "9.0")) == {"wrapped": "9.1.0"}


def test_the_tool_subprocess_sees_its_directory_first_and_others_do_not(tmp_path) -> None:
    sat = tmp_path / "sat"
    _script(sat, "printpath", 'echo "$PATH" > "$1"')
    bindirs.configure_bin_dirs({"mine": sat})

    out = tmp_path / "mine.txt"
    run_tool(_caps("mine", "printpath"), ["printpath", str(out)], logger=_LOG)
    assert out.read_text().split(os.pathsep)[0] == str(sat)

    other = tmp_path / "other.txt"
    run_tool(_caps("other", "sh"), ["sh", "-c", f'echo "$PATH" > {other}'], logger=_LOG)
    assert str(sat) not in other.read_text().split(os.pathsep)


def test_an_env_path_given_by_the_adapter_keeps_the_directory_first(tmp_path) -> None:
    sat = tmp_path / "sat"
    _script(sat, "printpath", 'echo "$PATH" > "$1"')
    bindirs.configure_bin_dirs({"mine": sat})
    out = tmp_path / "out.txt"
    env = {**os.environ, "OMP_NUM_THREADS": "2"}
    run_tool(_caps("mine", "printpath"), ["printpath", str(out)], logger=_LOG, env=env)
    assert out.read_text().split(os.pathsep)[0] == str(sat)


def test_unset_leaves_lookups_unchanged() -> None:
    assert bindirs.tool_path("gubbins") is None
    assert bindirs.host_which("sh", "gubbins") is not None


def test_used_directories_are_collected_per_tool(tmp_path) -> None:
    sat = tmp_path / "sat"
    _script(sat, "fakegubbins", 'echo "fakegubbins 3.4.1"')
    bindirs.configure_bin_dirs({"fakegub": sat, "unused": tmp_path})
    bindirs.reset_used()
    preflight(_caps("fakegub", "fakegubbins"))
    assert bindirs.used() == {"fakegub": str(sat)}
    bindirs.reset_used()
    assert bindirs.used() == {}


def test_the_stage_record_keeps_the_directories() -> None:
    record = StageRecord(tool="gubbins", bin_dirs={"gubbins": "/env/bin"})
    assert StageRecord.from_dict(record.to_dict()).bin_dirs == {"gubbins": "/env/bin"}
    # Without directories the record is written as before.
    assert "bin_dirs" not in StageRecord(tool="x").to_dict()


def test_gubbins_tree_builder_choice_uses_the_tool_path(tmp_path, monkeypatch) -> None:
    sat = tmp_path / "gub"
    _script(sat, "iqtree2", 'echo "IQ-TREE 2.4.0"')
    # No RAxML or IQ-TREE on the global PATH.
    monkeypatch.setenv("PATH", "/usr/bin:/bin")
    bindirs.configure_bin_dirs({"gubbins": sat})
    builder, threads = resolve_tree_builder(
        None, 4, _LOG, on_host=True, path=bindirs.tool_path("gubbins")
    )
    assert (builder, threads) == ("iqtree", 4)
    assert resolve_tree_builder(None, 4, _LOG, on_host=True) == (None, 1)


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("gubbins=/p/bin,mashtree=/q/bin", {"gubbins": "/p/bin", "mashtree": "/q/bin"}),
        (" gubbins = /p/bin , ", {"gubbins": "/p/bin"}),
        ("", {}),
    ],
)
def test_the_environment_variable_is_parsed(text, expected) -> None:
    assert bindirs.parse_entries(text.split(",")) == expected


@pytest.mark.parametrize("entry", ["gubbins", "=/p", "gubbins="])
def test_a_malformed_entry_is_rejected(entry) -> None:
    with pytest.raises(UserInputError, match="TOOL=DIR"):
        bindirs.parse_entries([entry])


def _list_tools(args, env=None):
    return _runner.invoke(app, [*args, "list-tools"], env=env)


def test_the_cli_option_and_variable_configure_the_directories(tmp_path) -> None:
    a, b = tmp_path / "a", tmp_path / "b"
    a.mkdir()
    b.mkdir()
    result = _list_tools(
        ["--bin-dir", f"gubbins={b}"],
        env={"REPGENR_BIN_DIRS": f"gubbins={a},mashtree={a}"},
    )
    assert result.exit_code == 0, result.output
    # The option overrides the variable for the same tool.
    assert bindirs.get_bin_dirs() == {"gubbins": b.resolve(), "mashtree": a.resolve()}


@pytest.mark.parametrize("via_env", [False, True])
def test_an_unknown_tool_is_a_usage_error(tmp_path, via_env) -> None:
    if via_env:
        result = _list_tools([], env={"REPGENR_BIN_DIRS": f"nosuchtool={tmp_path}"})
    else:
        result = _list_tools(["--bin-dir", f"nosuchtool={tmp_path}"])
    assert result.exit_code == 2
    assert "nosuchtool" in result.output


def test_a_missing_directory_is_a_usage_error(tmp_path) -> None:
    result = _list_tools(["--bin-dir", f"gubbins={tmp_path / 'absent'}"])
    assert result.exit_code == 2
    assert "absent" in result.output


def test_non_registry_tools_are_accepted(tmp_path) -> None:
    result = _list_tools(["--bin-dir", f"checkm2={tmp_path}", "--bin-dir", f"datasets={tmp_path}"])
    assert result.exit_code == 0, result.output


def test_a_directory_for_a_tool_that_runs_in_an_image_is_reported(tmp_path) -> None:
    result = _list_tools(["--container", "docker", "--bin-dir", f"gubbins={tmp_path}"])
    assert result.exit_code == 0, result.output
    assert "--bin-dir gubbins has no effect" in result.stderr
    # skder has no pinned image, so without --wave it runs on the host.
    result = _list_tools(["--container", "docker", "--bin-dir", f"skder={tmp_path}"])
    assert "no effect" not in result.stderr


def test_bin_dirs_raise_on_bad_input_outside_the_cli(tmp_path) -> None:
    with pytest.raises(UserInputError, match="nosuch"):
        bindirs.validate({"nosuch": str(tmp_path)})


def test_the_stage_record_names_the_directory_a_tool_used(
    tmp_path, workdir, genome_files, register_tool
) -> None:
    from repgenr.core.config import Config
    from repgenr.snptypers.base import SnpResult, SnpTyper
    from repgenr.snptypers.base import registry as snptypers

    sat = tmp_path / "sat"
    _script(sat, "dirtyper", 'echo "dirtyper 1.2"')

    class _DirTyper(SnpTyper):
        capabilities = _caps("dirtyper", "dirtyper", "1.0")
        requires_reference = False

        def call(self, genomes, reference, out_dir, params, logger) -> SnpResult:  # noqa: ANN001
            core = out_dir / "core.fasta"
            core.write_text("".join(f">{g.stem}\nACGT\n" for g in genomes))
            return SnpResult(core_snp_fasta=core)

    register_tool(snptypers, "dirtyper", _DirTyper)
    result = _runner.invoke(
        app,
        [
            "--bin-dir",
            f"dirtyper={sat}",
            "snptype",
            "-wd",
            str(workdir),
            "--tool",
            "dirtyper",
            "--all-genomes",
        ],
    )
    assert result.exit_code == 0, result.output
    record = Config.load(workdir).stages["snptype"]
    assert record.bin_dirs == {"dirtyper": str(sat.resolve())}
    assert record.tool_versions == {"dirtyper": "1.2.0"}
