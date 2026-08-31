"""process.run stdout_path is atomic: visible only when the command succeeds."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import pytest

from repgenr.core import process
from repgenr.core.errors import ToolExecutionError

_LOG = logging.getLogger("test")


def test_stdout_path_written_on_success(tmp_path: Path) -> None:
    out = tmp_path / "tree.nwk"
    process.run(
        [sys.executable, "-c", "print('(a,b);')"], logger=_LOG, stdout_path=out
    )
    assert out.read_text(encoding="utf-8").strip() == "(a,b);"
    assert list(tmp_path.iterdir()) == [out]  # no .part leftovers


def test_failed_command_preserves_previous_stdout_file(tmp_path: Path) -> None:
    out = tmp_path / "tree.nwk"
    out.write_text("(previous,good);\n", encoding="utf-8")
    with pytest.raises(ToolExecutionError):
        process.run(
            [sys.executable, "-c", "import sys; print('(part'); sys.exit(3)"],
            logger=_LOG, stdout_path=out,
        )
    assert out.read_text(encoding="utf-8") == "(previous,good);\n"
    assert list(tmp_path.iterdir()) == [out]


def test_timed_out_command_preserves_previous_stdout_file(tmp_path: Path) -> None:
    out = tmp_path / "tree.nwk"
    out.write_text("(previous,good);\n", encoding="utf-8")
    with pytest.raises(ToolExecutionError, match="timeout"):
        process.run(
            [sys.executable, "-c", "import time; print('(part'); time.sleep(30)"],
            logger=_LOG, stdout_path=out, timeout=0.5,
        )
    assert out.read_text(encoding="utf-8") == "(previous,good);\n"
    assert list(tmp_path.iterdir()) == [out]


def test_failed_command_with_no_previous_file_leaves_nothing(tmp_path: Path) -> None:
    out = tmp_path / "tree.nwk"
    with pytest.raises(ToolExecutionError):
        process.run(
            [sys.executable, "-c", "import sys; sys.exit(1)"],
            logger=_LOG, stdout_path=out,
        )
    assert not out.exists()
    assert list(tmp_path.iterdir()) == []


def test_check_false_nonzero_exit_still_publishes_stdout(tmp_path: Path) -> None:
    """check=False callers read the exit code themselves; stdout must be there."""
    out = tmp_path / "probe.txt"
    rc = process.run(
        [sys.executable, "-c", "import sys; print('partial-info'); sys.exit(5)"],
        logger=_LOG, stdout_path=out, check=False,
    )
    assert rc == 5
    assert out.read_text(encoding="utf-8").strip() == "partial-info"
