"""Subprocess timeout handling in core.process.run."""

from __future__ import annotations

import logging
import sys
import time

import pytest

from repgenr.core import process
from repgenr.core.errors import ToolExecutionError

_LOG = logging.getLogger("test")


def test_timeout_kills_and_raises() -> None:
    start = time.monotonic()
    with pytest.raises(ToolExecutionError, match="timeout"):
        process.run(
            [sys.executable, "-c", "import time; time.sleep(30)"],
            logger=_LOG, timeout=0.5,
        )
    # killed promptly, nowhere near the 30s sleep
    assert time.monotonic() - start < 10


def test_no_timeout_completes() -> None:
    rc = process.run([sys.executable, "-c", "print('ok')"], logger=_LOG, timeout=10)
    assert rc == 0


def test_env_default_timeout(monkeypatch) -> None:
    monkeypatch.setenv("REPGENR_SUBPROCESS_TIMEOUT", "0.5")
    with pytest.raises(ToolExecutionError, match="timeout"):
        process.run(
            [sys.executable, "-c", "import time; time.sleep(30)"], logger=_LOG
        )


def test_env_default_unset_means_no_timeout(monkeypatch) -> None:
    monkeypatch.delenv("REPGENR_SUBPROCESS_TIMEOUT", raising=False)
    assert process._default_timeout() is None
    # a quick command still runs fine with no timeout configured
    assert process.run([sys.executable, "-c", "pass"], logger=_LOG) == 0
