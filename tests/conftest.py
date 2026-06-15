"""Shared pytest fixtures."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest


def pytest_runtest_setup(item: pytest.Item) -> None:
    """Skip tests marked ``requires_binary("name")`` when the tool is absent."""
    for marker in item.iter_markers(name="requires_binary"):
        for binary in marker.args:
            if shutil.which(binary) is None:
                pytest.skip(f"requires external binary '{binary}'")


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
