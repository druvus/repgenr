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


@pytest.mark.parametrize("raw", ["abc", "-5", "0", "1.5s"])
def test_env_malformed_value_warns_once_and_proceeds(monkeypatch, caplog, raw) -> None:
    monkeypatch.setenv("REPGENR_SUBPROCESS_TIMEOUT", raw)
    monkeypatch.setattr(process, "_warned_bad_timeout", False)
    with caplog.at_level(logging.WARNING, logger="repgenr.core.process"):
        assert process._default_timeout() is None
        assert process._default_timeout() is None
    warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert "REPGENR_SUBPROCESS_TIMEOUT" in warnings[0].getMessage()
    assert raw in warnings[0].getMessage()


def test_env_valid_value_does_not_warn(monkeypatch, caplog) -> None:
    monkeypatch.setenv("REPGENR_SUBPROCESS_TIMEOUT", "12.5")
    monkeypatch.setattr(process, "_warned_bad_timeout", False)
    with caplog.at_level(logging.WARNING, logger="repgenr.core.process"):
        assert process._default_timeout() == 12.5
    assert not [r for r in caplog.records if r.levelno == logging.WARNING]
