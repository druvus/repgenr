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
    # Scope the fake entry points to this test group only: a test that loads
    # the REAL registries while this patch is active (e.g. list-tools) must not
    # cache fake adapters into them for the rest of the pytest process.
    real_entry_points = plugins.entry_points

    def fake_entry_points(group):
        if group == "repgenr.test_tools":
            return [_GoodEP(), _BadEP()]
        return real_entry_points(group=group)

    monkeypatch.setattr(plugins, "entry_points", fake_entry_points)
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


def test_register_and_unregister_roundtrip(register_tool) -> None:
    from repgenr.core.errors import PluginError
    from repgenr.core.plugins import Registry, ToolCapabilities

    reg: Registry = Registry("repgenr.test_register")
    reg._loaded = True  # bare registry: skip entry-point discovery

    class Adapter:
        capabilities = ToolCapabilities(name="mytool")

    register_tool(reg, "mytool", Adapter)
    assert "mytool" in reg.names()
    assert reg.get("mytool") is Adapter
    with pytest.raises(PluginError, match="already registered"):
        reg.register("mytool", Adapter)
    reg.register("mytool", Adapter, replace=True)  # explicit override allowed


def test_third_party_adapter_without_run_tool_survives_contract_fixture(
    register_tool, monkeypatch, tmp_path
) -> None:
    """An external adapter module lacking a module-level run_tool symbol must
    not break the contract suite's patching (raising=False)."""
    import types

    from repgenr.dereplicators.base import Dereplicator, registry

    mod = types.ModuleType("thirdparty_derep")

    class ThirdParty(Dereplicator):
        from repgenr.core.plugins import ToolCapabilities

        capabilities = ToolCapabilities(name="thirdparty")

        def dereplicate(self, genomes, out_dir, params, logger):
            raise NotImplementedError

    ThirdParty.__module__ = "thirdparty_derep"
    import sys

    monkeypatch.setitem(sys.modules, "thirdparty_derep", mod)
    register_tool(registry, "thirdparty", ThirdParty)

    # mimic the contract fixture's patch loop over every registered adapter
    for name in registry.names():
        module = sys.modules[registry.get(name).__module__]
        monkeypatch.setattr(module, "run_tool", lambda *a, **k: 0, raising=False)


def test_glance_rejects_tool_without_compare(tmp_path) -> None:
    from repgenr.core.context import WorkdirContext
    from repgenr.core.errors import UserInputError
    from repgenr.stages.glance import GlanceParams
    from repgenr.stages.glance import run as glance_run

    ctx = WorkdirContext(tmp_path, create=True)
    ctx.genomes_dir.mkdir()
    (ctx.genomes_dir / "g.fasta").write_text(">g\nAC\n", encoding="utf-8")
    with pytest.raises(UserInputError, match="drep"):
        glance_run(ctx, GlanceParams(tool="skder"))  # skder has no compare()


def test_outgroup_builder_must_support_distance_matrix() -> None:
    from repgenr.core.errors import UserInputError
    from repgenr.viral._outgroup import resolve_outgroup_builder

    builder = resolve_outgroup_builder("mashtree")
    assert builder.capabilities.name == "mashtree"
    with pytest.raises(UserInputError, match="mashtree"):
        resolve_outgroup_builder("iqtree")  # no distance_matrix support


def test_masker_family_is_registered() -> None:
    from repgenr.maskers.base import registry as masker_registry

    assert "gubbins" in masker_registry.names()


def test_snptype_dispatches_masker_via_registry(tmp_path, monkeypatch, register_tool) -> None:
    """The snptype stage's --mask goes through the masker registry, so a new
    masker is one adapter + one entry point away."""

    from repgenr.core.plugins import ToolCapabilities
    from repgenr.maskers.base import Masker, registry
    from repgenr.snptypers.base import SnpResult
    from repgenr.stages import snptype as snptype_mod
    from repgenr.stages.snptype import SnptypeParams, snptype_core

    class FakeMasker(Masker):
        capabilities = ToolCapabilities(name="fakemask")

        def mask(self, full_alignment, out_dir, params, logger):
            out_dir.mkdir(parents=True, exist_ok=True)
            out = out_dir / "filtered.fasta"
            out.write_text(">m\nAC\n", encoding="utf-8")
            return out

    register_tool(registry, "fakemask", FakeMasker)

    class FakeTyper:
        requires_reference = False
        capabilities = ToolCapabilities(name="faketyper")

        def preflight(self):
            return {}

        def call(self, genomes, ref, out_dir, params, logger):
            core = out_dir / "raw.fasta"
            full = out_dir / "full.fasta"
            core.write_text(">a\nACGT\n", encoding="utf-8")
            full.write_text(">a\nACGTACGT\n", encoding="utf-8")
            return SnpResult(core_snp_fasta=core, full_alignment=full)

    monkeypatch.setattr(snptype_mod.snp_registry, "create", lambda name: FakeTyper())
    genomes = [tmp_path / "g1.fasta"]
    genomes[0].write_text(">g1\nACGT\n", encoding="utf-8")
    result, versions = snptype_core(
        genomes,
        None,
        tmp_path / "snp",
        tmp_path / "scratch",
        SnptypeParams(tool="whatever", mask="fakemask"),
        logging.getLogger("test"),
    )
    assert result.masked is True
    assert (tmp_path / "snp" / "core_snp.fasta").read_text(encoding="utf-8") == ">m\nAC\n"


def test_unknown_masker_lists_available(tmp_path, monkeypatch) -> None:
    from repgenr.core.errors import PluginError
    from repgenr.core.plugins import ToolCapabilities
    from repgenr.snptypers.base import SnpResult
    from repgenr.stages import snptype as snptype_mod
    from repgenr.stages.snptype import SnptypeParams, snptype_core

    class FakeTyper:
        requires_reference = False
        capabilities = ToolCapabilities(name="faketyper")

        def preflight(self):
            return {}

        def call(self, genomes, ref, out_dir, params, logger):
            core = out_dir / "raw.fasta"
            core.write_text(">a\nACGT\n", encoding="utf-8")
            return SnpResult(core_snp_fasta=core)

    monkeypatch.setattr(snptype_mod.snp_registry, "create", lambda name: FakeTyper())
    genomes = [tmp_path / "g1.fasta"]
    genomes[0].write_text(">g1\nACGT\n", encoding="utf-8")
    with pytest.raises(PluginError, match="gubbins"):
        snptype_core(
            genomes,
            None,
            tmp_path / "snp",
            tmp_path / "scratch",
            SnptypeParams(tool="x", mask="nosuchmask"),
            logging.getLogger("test"),
        )
