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


def test_reads_and_assemble_a_public_illumina_run(run_repgenr, tmp_path: Path) -> None:
    """SRR25474756: Mycoplasmopsis arginini, MiSeq paired, 134 MB, about 300x.

    Ends after assembly: one genome is not a tree, and the tail of the chain
    is covered by the other live suites.
    """
    wd = tmp_path / "wd"
    listing = tmp_path / "runs.txt"
    listing.write_text("SRR25474756\n", encoding="utf-8")
    run_repgenr("reads", "-wd", wd, "--accession-file", listing)
    run_repgenr(*DOCKER, "assemble", "-wd", wd, "--assembler", "skesa", "-t", "4", timeout=3600)
    rows = read_selection(wd / "selection.tsv")
    assert [r.accession for r in rows] == ["SRR25474756"]
    assert rows[0].genus == "Mycoplasmopsis" and rows[0].species == "arginini"
    stats = read_assembly_stats(wd / "assembly_stats.tsv")[0]
    assert 500_000 < stats.total_length < 1_000_000  # M. arginini is about 0.7 Mb
    assert stats.n50 > 20_000 and stats.est_coverage > 100
    assert (wd / "genomes" / rows[0].filename).exists()
    assert not (wd / "scratch" / "assemble" / "SRR25474756").exists()  # reads removed


def test_reads_steps_assemble_a_public_run(run_repgenr, tmp_path: Path) -> None:
    """The stateless steps behind the Nextflow reads mode, on SRR25474756:
    assemble-run for the one run, then reads-gather into the genome contract
    (genome-qc needs a reference database, which the audit machine lacks)."""
    wd = tmp_path / "wd"
    listing = tmp_path / "runs.txt"
    listing.write_text("SRR25474756\n", encoding="utf-8")
    run_repgenr("reads", "-wd", wd, "--accession-file", listing)
    assemblies = tmp_path / "assemblies"
    versions = tmp_path / "versions.yml"
    run_repgenr(
        *DOCKER,
        "assemble-run",
        "--reads-tsv",
        wd / "reads.tsv",
        "--run",
        "SRR25474756",
        "-o",
        assemblies / "SRR25474756",
        "--assembler",
        "skesa",
        "-t",
        "4",
        "--memory-gb",
        "8",
        "--min-contig-length",
        "1000",
        "--versions-out",
        versions,
        timeout=3600,
    )
    assert (assemblies / "SRR25474756" / "assembly.ok").exists()
    assert "skesa:" in versions.read_text(encoding="utf-8")
    out = tmp_path / "out"
    run_repgenr(
        "reads-gather", "--reads-tsv", wd / "reads.tsv", "--assemblies", assemblies, "-o", out
    )
    rows = read_selection(out / "selection.tsv")
    assert [r.accession for r in rows] == ["SRR25474756"]
    assert (out / "genomes" / rows[0].filename).exists()
    stats = read_assembly_stats(out / "assembly_stats.tsv")[0]
    assert stats.assembler == "skesa" and stats.n50 > 20_000
    assert (out / "outgroup_accession.txt").read_text(encoding="utf-8") == ""
