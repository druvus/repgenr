"""The reads chain against real services and tools: ENA discovery, an ENA
download with checksum verification, and assembly in the pinned images."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from repgenr.core.contracts import read_assembly_stats, read_selection

pytestmark = [pytest.mark.live, pytest.mark.network, pytest.mark.container]

DOCKER = ["--container", "docker", "--platform", "linux/amd64"]


@pytest.fixture(scope="module", autouse=True)
def docker_preflight() -> None:
    if shutil.which("docker") is None:
        pytest.skip("docker binary not found")
    if subprocess.run(["docker", "info"], capture_output=True).returncode != 0:
        pytest.skip("docker daemon not running")


def test_run_reads_assembles_a_public_illumina_run(run_repgenr, tmp_path: Path) -> None:
    """SRR25474756: Mycoplasmopsis arginini, MiSeq paired, 134 MB, about 300x."""
    wd = tmp_path / "wd"
    listing = tmp_path / "runs.txt"
    listing.write_text("SRR25474756\n", encoding="utf-8")
    run_repgenr(
        *DOCKER,
        "run",
        "-wd",
        wd,
        "--reads",
        "--accession-file",
        listing,
        "--assembler",
        "skesa",
        "-t",
        "4",
        "--treebuilder",
        "mashtree",
        "--no-outgroup",
        "--tool",
        "sourmash",
        timeout=3600,
    )
    rows = read_selection(wd / "selection.tsv")
    assert [r.accession for r in rows] == ["SRR25474756"]
    assert rows[0].genus == "Mycoplasmopsis" and rows[0].species == "arginini"
    stats = read_assembly_stats(wd / "assembly_stats.tsv")[0]
    assert 500_000 < stats.total_length < 1_000_000  # M. arginini is about 0.7 Mb
    assert stats.n50 > 20_000 and stats.est_coverage > 100
    assert (wd / "tree" / "tree.nwk").exists() and (wd / "tree2tax.tsv").exists()
