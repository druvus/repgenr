"""The branchwater availability probe is memoized (one container start, not N)."""

from __future__ import annotations

import logging

from repgenr.dereplicators import sourmash
from repgenr.dereplicators.sourmash import SourmashDereplicator, _branchwater_available

_LOG = logging.getLogger("test")


def test_branchwater_probe_memoized(monkeypatch) -> None:
    sourmash._BRANCHWATER_CACHE.clear()
    calls = {"n": 0}

    def fake_run_tool(caps, cmd, **kwargs):
        calls["n"] += 1
        return 0

    monkeypatch.setattr(sourmash, "run_tool", fake_run_tool)
    caps = SourmashDereplicator().capabilities

    assert _branchwater_available(caps, _LOG) is True
    assert _branchwater_available(caps, _LOG) is True
    assert _branchwater_available(caps, _LOG) is True
    assert calls["n"] == 1  # probed once, cached thereafter
