"""Tests for the in-stage parallel_map helper."""

from __future__ import annotations

import threading
import time

import pytest

from repgenr.core.executors import parallel_map


def test_order_preserved() -> None:
    assert parallel_map(lambda x: x * 2, [1, 2, 3, 4], workers=3) == [2, 4, 6, 8]


def test_sequential_when_single_worker() -> None:
    assert parallel_map(str, [1, 2, 3], workers=1) == ["1", "2", "3"]


def test_empty() -> None:
    assert parallel_map(lambda x: x, [], workers=4) == []


def test_exception_propagates() -> None:
    def boom(x):
        if x == 2:
            raise ValueError("boom")
        return x

    with pytest.raises(ValueError, match="boom"):
        parallel_map(boom, [1, 2, 3], workers=2)


def test_actually_concurrent() -> None:
    # 4 tasks that each sleep; with >=4 workers they overlap, so max concurrency > 1.
    active = 0
    peak = 0
    lock = threading.Lock()

    def task(_):
        nonlocal active, peak
        with lock:
            active += 1
            peak = max(peak, active)
        time.sleep(0.05)
        with lock:
            active -= 1
        return 1

    parallel_map(task, list(range(4)), workers=4)
    assert peak >= 2
