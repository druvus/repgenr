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


def test_probe_cache_keyed_by_config_values(monkeypatch) -> None:
    # Rebinding the global container config to an equal-valued object must hit
    # the cache: identity-based keys (id()) can produce stale hits after the
    # old object's address is reused, and spurious misses otherwise.
    from repgenr.core.containers import configure_container

    sourmash._BRANCHWATER_CACHE.clear()
    calls = {"n": 0}

    def fake_run_tool(caps, cmd, **kwargs):
        calls["n"] += 1
        return 0

    monkeypatch.setattr(sourmash, "run_tool", fake_run_tool)
    caps = SourmashDereplicator().capabilities
    try:
        configure_container("none")
        assert _branchwater_available(caps, _LOG) is True
        configure_container("none")  # new object, same values
        assert _branchwater_available(caps, _LOG) is True
        assert calls["n"] == 1
        configure_container("docker")  # different values -> re-probe
        assert _branchwater_available(caps, _LOG) is True
        assert calls["n"] == 2
    finally:
        configure_container("none")
        sourmash._BRANCHWATER_CACHE.clear()
