"""Broken plugins are visible: load warnings, auto-select logging, listings (WP1-6)."""

from __future__ import annotations

import logging

import pytest

from repgenr.core import plugins
from repgenr.core.errors import PluginError
from repgenr.core.plugins import Registry, ToolCapabilities


class _GoodAdapter:
    capabilities = ToolCapabilities(name="goodtool")


class _GoodEP:
    name = "goodtool"

    def load(self):
        return _GoodAdapter


class _BadEP:
    name = "badtool"

    def load(self):
        raise ImportError("boom: missing dependency")


@pytest.fixture()
def registry(monkeypatch) -> Registry:
    monkeypatch.setattr(plugins, "entry_points", lambda group: [_GoodEP(), _BadEP()])
    return Registry("repgenr.test_tools")


def test_broken_plugin_import_logged_at_warning(registry, caplog) -> None:
    with caplog.at_level(logging.WARNING, logger="repgenr"):
        names = registry.names()
    assert names == ["badtool", "goodtool"]
    warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert any("badtool" in r.getMessage() for r in warnings)


def test_registry_reports_broken_names(registry) -> None:
    registry.names()  # trigger load
    assert registry.is_broken("badtool")
    assert not registry.is_broken("goodtool")
    assert not registry.is_broken("nonexistent")


def test_get_broken_plugin_still_raises(registry) -> None:
    with pytest.raises(PluginError, match="badtool"):
        registry.get("badtool")


def test_auto_select_logs_skipped_candidates(registry, caplog) -> None:
    registry.names()
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="repgenr"):
        chosen = plugins.auto_select(registry, 10)
    assert chosen == "goodtool"
    skipped = [r for r in caplog.records if "badtool" in r.getMessage()]
    assert skipped, "auto_select should log the skipped broken candidate"


def test_list_tools_marks_broken_entries(registry, monkeypatch) -> None:
    from typer.testing import CliRunner

    import repgenr.dereplicators.base as derep_base
    from repgenr.cli.main import app

    monkeypatch.setattr(derep_base, "registry", registry)
    result = CliRunner().invoke(app, ["list-tools"])
    assert result.exit_code == 0
    assert "badtool (broken)" in result.output
    assert "goodtool" in result.output
