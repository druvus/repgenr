"""Streamed tool output goes to DEBUG; the command line stays at INFO."""

from __future__ import annotations

import logging
import sys

from repgenr.core import process

_LOG = logging.getLogger("test.process.logging")


def test_tool_output_lines_are_debug_and_command_is_info(caplog) -> None:
    with caplog.at_level(logging.DEBUG, logger=_LOG.name):
        process.run(
            [sys.executable, "-c", "import sys; print('step 1/3'); print('done', file=sys.stderr)"],
            logger=_LOG,
            log_prefix="fake",
        )
    by_level = {(r.levelname, r.message) for r in caplog.records}
    info = [msg for lvl, msg in by_level if lvl == "INFO"]
    assert len(info) == 1 and info[0].startswith("[fake] $ ")  # only the command line
    assert ("DEBUG", "[fake] step 1/3") in by_level
    assert ("DEBUG", "[fake] done") in by_level


def test_carriage_return_redraws_are_reduced_to_the_last_frame(caplog) -> None:
    # A progress bar redraws one line with \r; only its final state is worth a log line.
    script = "import sys; sys.stdout.write('a 10%\\ra 50%\\ra 100%\\n'); sys.stdout.flush()"
    with caplog.at_level(logging.DEBUG, logger=_LOG.name):
        process.run([sys.executable, "-c", script], logger=_LOG, log_prefix="bar")
    messages = [r.message for r in caplog.records if r.levelname == "DEBUG"]
    assert messages == ["[bar] a 100%"]


def test_tools_get_no_stdin(monkeypatch) -> None:
    # A tool that prompts must read end-of-file, not wait on the terminal.
    # pytest already replaces stdin, so assert on the Popen argument itself.
    import subprocess

    seen: dict = {}
    real_popen = subprocess.Popen

    def spy(*args, **kwargs):
        seen.update(kwargs)
        return real_popen(*args, **kwargs)

    monkeypatch.setattr(process.subprocess, "Popen", spy)
    process.run([sys.executable, "-c", "pass"], logger=_LOG, log_prefix="fake")
    assert seen.get("stdin") is subprocess.DEVNULL
