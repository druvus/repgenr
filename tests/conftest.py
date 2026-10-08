"""Shared pytest fixtures."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    # Registered here (an initial conftest under `testpaths`) so the option is
    # known from any invocation; the live suite's own conftest consumes it.
    parser.addoption(
        "--live-config",
        default=None,
        help="TOML mapping tools to bin dirs and setting cache_dir for tests/live "
        "(see tests/live/live.example.toml).",
    )


def pytest_runtest_setup(item: pytest.Item) -> None:
    """Skip tests marked ``requires_binary("name")`` when the tool is absent."""
    for marker in item.iter_markers(name="requires_binary"):
        for binary in marker.args:
            if shutil.which(binary) is None:
                pytest.skip(f"requires external binary '{binary}'")


class _ReachableSession:
    """Stands in for the reachability probe's session: every host answers."""

    def get(self, url, **kwargs):
        return None

    def close(self) -> None:
        pass


@pytest.fixture(autouse=True)
def _hermetic_reachability_probe(request: pytest.FixtureRequest, monkeypatch) -> None:
    """Keep the core.http reachability probe off the network outside the live suite.

    Stages probe the NCBI datasets host before calling ``datasets``; unit and
    integration tests fake ``datasets`` itself, so the probe answers as if the
    host were reachable. Tests of the probe replace ``_probe_session`` again.
    """
    if request.node.get_closest_marker("live"):
        return
    from repgenr.core import http

    monkeypatch.setattr(http, "_probe_session", _ReachableSession)


@pytest.fixture
def workdir(tmp_path: Path) -> Path:
    return tmp_path / "wd"


@pytest.fixture
def genome_files(workdir: Path) -> list[Path]:
    """Create a small genomes/ dir with RepGenR-style filenames."""
    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    names = [
        "Francisellaceae_francisella_tularensis_GCA_000001.fasta",
        "Francisellaceae_francisella_tularensis_GCA_000002.fasta",
        "Francisellaceae_francisella_tularensis_GCA_000003.fasta",
    ]
    out = []
    for i, name in enumerate(names):
        p = gdir / name
        p.write_text(f">seq{i}\n{'ACGT' * 10}\n")
        out.append(p)
    return out


@pytest.fixture()
def register_tool():
    """Register an adapter class on a registry for one test, restoring after."""
    registered: list[tuple] = []

    def _register(registry, name, cls):
        registry.register(name, cls, replace=True)
        registered.append((registry, name))
        return cls

    yield _register
    for registry, name in registered:
        registry.unregister(name)


@pytest.fixture()
def write_deliverables():
    """Create placeholder outputs for each recorded stage of a workdir.

    Tests that build repgenr.yaml by hand and check `status` need the
    declared deliverables on disk, or status reports the stages as stale.
    """
    from types import SimpleNamespace

    from repgenr.cli.base import STAGE_DELIVERABLES
    from repgenr.core.config import Config
    from repgenr.core.doctor import _layout

    def _write(wd: Path) -> None:
        ctx = _layout(wd)
        for name, record in Config.load(wd).stages.items():
            spec = STAGE_DELIVERABLES.get(name)
            for path in spec(ctx, SimpleNamespace(**record.params)) if spec else []:
                if path.suffix:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    if not path.exists():
                        path.write_text("x\n", encoding="utf-8")
                else:
                    path.mkdir(parents=True, exist_ok=True)
                    (path / "placeholder.fasta").write_text(">x\nA\n", encoding="utf-8")

    return _write
