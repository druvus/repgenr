"""Long command lines and name lists are shortened on the console only.

The console handler prints a short form of a message that carries one, and
the run log keeps the full text.
"""

from __future__ import annotations

import logging
import os
import stat
from pathlib import Path

import pytest

from repgenr.core import process
from repgenr.core.logging import capped_names, configure_logging, console_extra

_CONSOLE_LIMIT = 200


@pytest.fixture(autouse=True)
def _reset_repgenr_logger():
    yield
    logger = logging.getLogger("repgenr")
    for handler in logger.handlers:
        handler.close()
    logger.handlers.clear()


def _run_logger(tmp_path: Path, level: int = logging.INFO) -> tuple[logging.Logger, Path]:
    # Configured inside the test body: the console handler binds the stderr
    # that capsys replaces for the call phase.
    workdir = tmp_path / "wd"
    return configure_logging(workdir, level=level), workdir / "repgenr.log"


def _fake_tool(tmp_path: Path) -> Path:
    tool = tmp_path / "bin" / "faketool"
    tool.parent.mkdir()
    tool.write_text("#!/bin/sh\necho 'fake tool output line'\n")
    tool.chmod(tool.stat().st_mode | stat.S_IXUSR)
    return tool


def _inputs(tmp_path: Path, n: int) -> list[str]:
    deep = tmp_path / "a_long_directory_name" / "genomes"
    deep.mkdir(parents=True)
    return [str(deep / f"Genus_species_g{i:06d}_GCF_{i:09d}.1.fasta") for i in range(n)]


def test_command_line_is_short_on_console_and_full_in_log(tmp_path, capsys) -> None:
    logger, log_file = _run_logger(tmp_path)
    tool = _fake_tool(tmp_path)
    cmd = [str(tool), "-g", *_inputs(tmp_path, 500), "-o", str(tmp_path / "out"), "-c", "4"]
    process.run(cmd, logger=logger, log_prefix="faketool")
    for handler in logger.handlers:
        handler.flush()

    err = capsys.readouterr().err
    command_lines = [ln for ln in err.splitlines() if "[faketool] $ " in ln]
    assert len(command_lines) == 1
    line = command_lines[0]
    assert len(line) < _CONSOLE_LIMIT, line
    assert "500" in line  # the inputs are counted
    assert "repgenr.log" in line  # and the reader is told where the full line is
    assert "-c 4" in line  # options after the inputs are kept
    assert "fake tool output line" not in err  # tool output stays off the console

    log_text = log_file.read_text()
    assert " ".join(cmd) in log_text
    assert "fake tool output line" in log_text


def test_verbose_console_also_gets_the_short_command(tmp_path, capsys) -> None:
    logger, log_file = _run_logger(tmp_path, logging.DEBUG)
    tool = _fake_tool(tmp_path)
    cmd = [str(tool), *_inputs(tmp_path, 500)]
    process.run(cmd, logger=logger, log_prefix="faketool")
    for handler in logger.handlers:
        handler.flush()
    err = capsys.readouterr().err
    assert all(len(ln) < _CONSOLE_LIMIT for ln in err.splitlines()), err
    assert " ".join(cmd) in log_file.read_text()


def test_short_command_is_unchanged_on_console(tmp_path, capsys) -> None:
    logger, _ = _run_logger(tmp_path)
    process.run(["sh", "-c", "exit 0"], logger=logger, log_prefix="faketool")
    err = capsys.readouterr().err
    assert "[faketool] $ sh -c exit 0\n" in err
    assert "repgenr.log" not in err


def test_without_a_run_log_the_console_keeps_the_full_text(capsys) -> None:
    logger = configure_logging(None, level=logging.INFO)
    logger.info("full text", extra=console_extra("short"))
    err = capsys.readouterr().err
    assert "full text" in err and "short" not in err


def test_shorten_command_counts_a_run_of_paths() -> None:
    paths = [f"/data/set/genome_{i}.fasta" for i in range(50)]
    short = process.shorten_command(["skder", "-g", *paths, "-o", "/tmp/out", "-d", "greedy"])
    assert short.startswith("skder -g ")
    assert "50 paths" in short
    assert short.endswith("-o /tmp/out -d greedy")
    assert len(short) <= process.CONSOLE_COMMAND_CHARS


def test_shorten_command_caps_a_single_long_argument() -> None:
    short = process.shorten_command(["tool", "--list", ",".join(["x" * 20] * 100)])
    assert len(short) <= process.CONSOLE_COMMAND_CHARS
    assert short.endswith("...")


def test_capped_names() -> None:
    names = [f"g{i}" for i in range(12)]
    assert capped_names(names) == "g0, g1, g2, g3, g4 and 7 more"
    assert capped_names(names[:5]) == "g0, g1, g2, g3, g4"
    assert capped_names([]) == ""


def test_name_list_is_capped_on_console_and_complete_in_log(tmp_path, capsys) -> None:
    logger, log_file = _run_logger(tmp_path)
    names = [f"genome_{i:04d}.fasta" for i in range(300)]
    logger.warning(
        "Left out: %s",
        ", ".join(names),
        extra=console_extra(f"Left out: {capped_names(names)}"),
    )
    for handler in logger.handlers:
        handler.flush()
    err = capsys.readouterr().err
    assert "and 295 more" in err
    assert len(err.strip()) < _CONSOLE_LIMIT
    assert ", ".join(names) in log_file.read_text()


def test_console_extra_does_not_alter_the_record_for_other_handlers(tmp_path) -> None:
    # The console formatter must work on a copy: the file handler formats the
    # same record afterwards.
    logger, log_file = _run_logger(tmp_path)
    logger.info("long %s", "form", extra=console_extra("short form"))
    for handler in logger.handlers:
        handler.flush()
    assert "long form" in log_file.read_text()
    assert os.path.getsize(log_file) > 0


def test_unpack_missing_members_are_capped_on_console(tmp_path, capsys) -> None:
    from repgenr.stages.derep_unpack import _warn_missing

    logger, log_file = _run_logger(tmp_path)
    missing = [f"genome_{i:04d}.fasta" for i in range(40)]
    _warn_missing(missing, tmp_path / "genomes", logger)
    for handler in logger.handlers:
        handler.flush()
    err = capsys.readouterr().err
    assert "and 35 more" in err and "genome_0039.fasta" not in err
    assert ", ".join(missing) in log_file.read_text()
