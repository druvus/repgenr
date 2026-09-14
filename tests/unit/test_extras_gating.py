"""Extras reach adapters only when they read them; unread keys are named."""

from __future__ import annotations

import logging

from repgenr.cli.base import gated_extra
from repgenr.core.plugins import (
    Registry,
    ToolCapabilities,
    warn_ignored_params,
    warn_unconsumed_extras,
)
from repgenr.dereplicators.base import Dereplicator, DerepParams


class _Reads(Dereplicator):
    capabilities = ToolCapabilities(name="reads", accepted_extras=frozenset({"virus"}))

    def dereplicate(self, genomes, out_dir, params, logger):  # noqa: ANN001
        raise NotImplementedError


class _Ignores(Dereplicator):
    capabilities = ToolCapabilities(name="ignores")

    def dereplicate(self, genomes, out_dir, params, logger):  # noqa: ANN001
        raise NotImplementedError


def _registry(register_tool):
    reg: Registry = Registry("repgenr.dereplicators")
    register_tool(reg, "reads", _Reads)
    register_tool(reg, "ignores", _Ignores)
    return reg


def test_gated_extra_injects_for_reader(register_tool) -> None:
    reg = _registry(register_tool)
    assert gated_extra(reg, "reads", "virus", True) == {"virus": True}


def test_gated_extra_skips_for_ignorer(register_tool) -> None:
    reg = _registry(register_tool)
    assert gated_extra(reg, "ignores", "virus", True) == {}


def test_gated_extra_passes_through_for_auto(register_tool) -> None:
    reg = _registry(register_tool)
    assert gated_extra(reg, "auto", "virus", True) == {"virus": True}


def test_warn_names_unread_keys(caplog) -> None:
    caps = ToolCapabilities(name="t", accepted_extras=frozenset({"a"}))
    with caplog.at_level(logging.WARNING):
        warn_unconsumed_extras(
            caps, {"a": 1, "b": 2, "c": 3}, logging.getLogger("x"), family="Aligner"
        )
    msgs = [r.message for r in caplog.records]
    assert any("Aligner 't' ignores" in m and "b, c" in m for m in msgs)


def test_warn_silent_when_all_read(caplog) -> None:
    caps = ToolCapabilities(name="t", accepted_extras=frozenset({"a"}))
    with caplog.at_level(logging.WARNING):
        warn_unconsumed_extras(caps, {"a": 1}, logging.getLogger("x"), family="Aligner")
    assert not caplog.records


def test_warn_ignored_params_names_non_default_values(caplog) -> None:
    caps = ToolCapabilities(name="t", ignored_params=frozenset({"primary_ani", "threads"}))
    params = DerepParams(primary_ani=0.95, threads=16)  # threads at its default
    with caplog.at_level(logging.WARNING):
        warn_ignored_params(caps, params, logging.getLogger("x"), family="Dereplicator")
    msgs = [r.message for r in caplog.records]
    assert len(msgs) == 1
    assert "Dereplicator 't' does not use primary_ani" in msgs[0]
    assert "0.95" in msgs[0]
    assert "threads" not in msgs[0]


def test_warn_ignored_params_silent_at_defaults(caplog) -> None:
    caps = ToolCapabilities(name="t", ignored_params=frozenset({"primary_ani"}))
    with caplog.at_level(logging.WARNING):
        warn_ignored_params(caps, DerepParams(), logging.getLogger("x"), family="Dereplicator")
    assert not caplog.records


def test_warn_ignored_params_silent_when_nothing_declared(caplog) -> None:
    caps = ToolCapabilities(name="t")
    with caplog.at_level(logging.WARNING):
        warn_ignored_params(
            caps, DerepParams(primary_ani=0.5), logging.getLogger("x"), family="Dereplicator"
        )
    assert not caplog.records
