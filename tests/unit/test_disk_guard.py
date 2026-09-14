"""A shared free-disk guard: refuse under the floor, warn when the estimate exceeds free space."""

from __future__ import annotations

import logging
import shutil
from collections import namedtuple

import pytest

from repgenr.core.errors import WorkdirError
from repgenr.core.process import MIN_FREE_BYTES, check_free_disk

_Usage = namedtuple("_Usage", "total used free")


def test_refuses_below_the_floor(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(shutil, "disk_usage", lambda p: _Usage(10, 9, MIN_FREE_BYTES - 1))
    with pytest.raises(WorkdirError, match="free under"):
        check_free_disk(tmp_path, 10, logging.getLogger("t"), what="download 3 runs")


def test_warns_when_the_estimate_exceeds_free_space(tmp_path, monkeypatch, caplog) -> None:
    monkeypatch.setattr(shutil, "disk_usage", lambda p: _Usage(10, 0, MIN_FREE_BYTES + 10))
    with caplog.at_level(logging.WARNING):
        check_free_disk(tmp_path, MIN_FREE_BYTES * 5, logging.getLogger("t"), what="assembly")
    assert any("Low disk" in r.message for r in caplog.records)
