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
