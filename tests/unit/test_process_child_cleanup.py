"""A failure inside the output read loop must not leave the tool running."""

from __future__ import annotations

import logging
import sys
import time
from pathlib import Path

import pytest

from repgenr.core import process


class _RaisingLogger(logging.Logger):
    """Raises on the first streamed line, standing in for any consumer failure."""

    def debug(self, msg, *args, **kwargs):  # noqa: ANN001
        if args and args[-1] == "started":
            raise RuntimeError("consumer failed")
        super().debug(msg, *args, **kwargs)


def test_read_loop_failure_kills_the_child(tmp_path: Path) -> None:
    marker = tmp_path / "child_finished"
    script = (
        "import sys, time, pathlib; print('started', flush=True); "
        f"time.sleep(2); pathlib.Path({str(marker)!r}).write_text('done')"
    )
    logger = _RaisingLogger("test.child.cleanup")
    logger.setLevel(logging.DEBUG)
    started = time.monotonic()
    with pytest.raises(RuntimeError, match="consumer failed"):
        process.run([sys.executable, "-c", script], logger=logger)
    # The exception propagated promptly and the orphan never got to finish.
    assert time.monotonic() - started < 1.5
    time.sleep(2.5)
    assert not marker.exists(), "child kept running after the read loop failed"


def test_sigterm_stops_the_running_tool(tmp_path: Path) -> None:
    """A repgenr process ended by SIGTERM stops its tool and removes the partial
    output (live: FastTree kept running after its phylo run was terminated)."""
    import os
    import signal
    import subprocess

    started = tmp_path / "started"
    finished = tmp_path / "finished"
    tool = (
        "import pathlib, time; "
        f"pathlib.Path({str(started)!r}).write_text('x'); "
        f"time.sleep(3); pathlib.Path({str(finished)!r}).write_text('x')"
    )
    out = tmp_path / "tree.nwk"
    driver = (
        "import logging, sys\n"
        "from repgenr.core import process\n"
        "process.install_termination_handler()\n"
        f"process.run([sys.executable, '-c', {tool!r}], logger=logging.getLogger('t'),"
        f" stdout_path={str(out)!r})\n"
    )
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(sys.path)}
    proc = subprocess.Popen(
        [sys.executable, "-c", driver], stderr=subprocess.PIPE, env=env, text=True
    )
    deadline = time.monotonic() + 10
    while not started.exists():
        assert time.monotonic() < deadline, "the tool never started"
        time.sleep(0.05)
    proc.send_signal(signal.SIGTERM)
    _, stderr = proc.communicate(timeout=10)
    assert proc.returncode == 128 + signal.SIGTERM
    assert "stopped 1 running tool" in stderr
    time.sleep(3.5)
    assert not finished.exists(), "the tool kept running after repgenr was terminated"
    assert not out.exists() and not (tmp_path / "tree.nwk.part").exists()
