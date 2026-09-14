"""Tools run through `--container docker`: pinned BioContainers under amd64
emulation, Wave-minted images for conda-spec adapters, the cache and
platform flags and their environment variables, the engine override, and the
resume rule that a container run never reuses a native result.

Skipped as a whole unless Docker is up, `linux/amd64` emulation works and the
Wave CLI is on PATH (Apple Silicon: Docker Desktop on the Virtualization
framework with Rosetta).
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest
from live_helpers import derep_partition, log_text, newick_leaves, truth_of, truth_partition

from repgenr.core.config import Config
from repgenr.core.contracts import MSA_FASTA, TREE_NWK

pytestmark = [pytest.mark.live, pytest.mark.container]

DOCKER = ["--container", "docker", "--platform", "linux/amd64"]


def _docker_ready() -> str | None:
    if shutil.which("docker") is None:
        return "docker binary not found"
    if subprocess.run(["docker", "info"], capture_output=True).returncode != 0:
        return "docker daemon not running"
    probe = subprocess.run(
        ["docker", "run", "--rm", "--platform", "linux/amd64", "alpine", "uname", "-m"],
        capture_output=True,
        text=True,
    )
    if probe.stdout.strip() != "x86_64":
        return f"linux/amd64 emulation unavailable ({probe.stdout.strip() or probe.stderr.strip()})"
    if shutil.which("wave") is None:
        return "wave CLI not on PATH"
    return None


@pytest.fixture(scope="module", autouse=True)
def docker_preflight() -> None:
    reason = _docker_ready()
    if reason:
        pytest.skip(reason)


@pytest.fixture
def clonal_wd(synthetic_set, ingested_workdir):
    genomes = synthetic_set("clonal", n=8, length=50_000)
    return genomes, ingested_workdir(genomes)


@pytest.fixture
def balanced_wd(synthetic_set, ingested_workdir):
    genomes = synthetic_set("balanced", n=5, length=20_000)
    return genomes, ingested_workdir(genomes)


# --- Wave-minted images (conda specs) ---------------------------------------------


def test_skder_in_a_wave_container_with_cache_and_env(
    run_repgenr, clonal_wd, tmp_path: Path
) -> None:
    genomes, wd = clonal_wd
    cache = tmp_path / "container-cache"
    run_repgenr(
        *DOCKER,
        "--wave",
        "--container-cache",
        cache,
        "dereplicate",
        "-wd",
        wd,
        "--tool",
        "skder",
        "-t",
        "2",
    )
    assert derep_partition(wd) == truth_partition(truth_of(genomes))
    assert "docker run" in log_text(wd)
    rec = Config.load(wd).stages["dereplicate"]
    assert rec.completed

    # The same request through the environment variables must skip: same
    # backend, platform and Wave choice give the same execution identity.
    again = run_repgenr(
        "dereplicate",
        "-wd",
        wd,
        "--tool",
        "skder",
        "-t",
        "2",
        env={
            "REPGENR_CONTAINER": "docker",
            "REPGENR_CONTAINER_PLATFORM": "linux/amd64",
            "REPGENR_WAVE": "1",
            "REPGENR_CONTAINER_CACHE": str(cache),
        },
    )
    assert "skipping" in again.stderr + again.stdout


def test_native_result_is_not_reused_by_a_container_run(run_repgenr, clonal_wd) -> None:
    _, wd = clonal_wd
    run_repgenr("dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2")
    stamp = Config.load(wd).stages["dereplicate"].completed
    out = run_repgenr(*DOCKER, "--wave", "dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2")
    assert "skipping" not in out.stderr + out.stdout, (
        "the container identity is part of the fingerprint"
    )
    assert Config.load(wd).stages["dereplicate"].completed != stamp


def test_drep_in_a_wave_container_with_virus_extra(run_repgenr, clonal_wd) -> None:
    genomes, wd = clonal_wd
    run_repgenr(
        *DOCKER,
        "--wave",
        "dereplicate",
        "-wd",
        wd,
        "--tool",
        "drep",
        "--virus",
        "-t",
        "2",
    )
    assert derep_partition(wd) == truth_partition(truth_of(genomes))
    assert Config.load(wd).stages["dereplicate"].params["virus"] is True


def test_simple_typer_multi_tool_image(run_repgenr, clonal_wd) -> None:
    _, wd = clonal_wd
    run_repgenr(*DOCKER, "--wave", "dereplicate", "-wd", wd, "--tool", "skder", "-t", "2")
    run_repgenr(
        *DOCKER, "--wave", "snptype", "-wd", wd, "--tool", "simple", "--all-genomes", "-t", "2"
    )
    core = wd / "snp" / "core_snp.fasta"
    assert sum(1 for line in core.open(encoding="utf-8") if line.startswith(">")) == 8


def test_sibeliaz_in_a_wave_container(run_repgenr, balanced_wd) -> None:
    genomes, wd = balanced_wd
    run_repgenr(
        *DOCKER, "--wave", "dereplicate", "-wd", wd, "--tool", "skder", "-sani", "0.999", "-t", "2"
    )
    run_repgenr(
        *DOCKER,
        "--wave",
        "phylo",
        "-wd",
        wd,
        "--aligner",
        "sibeliaz",
        "--treebuilder",
        "fasttree",
        "--no-outgroup",
        "--aligner-arg",
        "kmer=15",
        "-t",
        "2",
    )
    msa = wd / "align" / MSA_FASTA
    assert sum(1 for line in msa.open(encoding="utf-8") if line.startswith(">")) == 5
    assert "-k 15" in log_text(wd), "--aligner-arg reaches the sibeliaz command line"
    assert len(newick_leaves((wd / "tree" / TREE_NWK).read_text(encoding="utf-8"))) == 5


# --- pinned BioContainers --------------------------------------------------------


def test_progressivemauve_pinned_image(run_repgenr, balanced_wd) -> None:
    genomes, wd = balanced_wd
    run_repgenr(*DOCKER, "dereplicate", "-wd", wd, "--tool", "skder", "-sani", "0.999", "-t", "2")
    run_repgenr(
        *DOCKER,
        "phylo",
        "-wd",
        wd,
        "--aligner",
        "progressivemauve",
        "--treebuilder",
        "fasttree",
        "--no-outgroup",
        "-t",
        "2",
    )
    assert "quay.io/biocontainers/mauve" in log_text(wd)
    msa = wd / "align" / MSA_FASTA
    assert sum(1 for line in msa.open(encoding="utf-8") if line.startswith(">")) == 5


def test_cactus_pinned_image(run_repgenr, synthetic_set, ingested_workdir) -> None:
    genomes = synthetic_set("balanced", n=4, length=50_000)
    wd = ingested_workdir(genomes)
    run_repgenr(*DOCKER, "dereplicate", "-wd", wd, "--tool", "skder", "-sani", "0.999", "-t", "2")
    run_repgenr(
        *DOCKER,
        "phylo",
        "-wd",
        wd,
        "--aligner",
        "cactus",
        "--treebuilder",
        "fasttree",
        "--no-outgroup",
        "-t",
        "4",
        timeout=3600,
    )
    assert "cactus:v2.9.3" in log_text(wd)
    msa = wd / "align" / MSA_FASTA
    assert sum(1 for line in msa.open(encoding="utf-8") if line.startswith(">")) == 4


def test_glance_drep_compare(run_repgenr, clonal_wd) -> None:
    _, wd = clonal_wd
    run_repgenr(*DOCKER, "--wave", "glance", "-wd", wd, "--tool", "drep", "-t", "2", "--keep-files")
    assert (wd / "glance_clustering_dendrogram.pdf").is_file()


# --- engine override --------------------------------------------------------------


def test_container_engine_podman_is_reported_when_missing(run_repgenr, clonal_wd) -> None:
    _, wd = clonal_wd
    if shutil.which("podman"):
        pytest.skip("podman is installed; the negative path does not apply")
    out = run_repgenr(
        "--container",
        "docker",
        "--container-engine",
        "podman",
        "--wave",
        "dereplicate",
        "-wd",
        wd,
        "--tool",
        "skder",
        check=False,
    )
    assert out.returncode != 0 and "podman" in out.stderr + out.stdout


def test_pinned_biocontainer_without_wave(run_repgenr, clonal_wd) -> None:
    """Without --wave a conda-spec adapter runs in its pinned BioContainer, not on the host."""
    _genomes, wd = clonal_wd
    run_repgenr(*DOCKER, "dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2")
    text = log_text(wd)
    assert "quay.io/biocontainers/sourmash:" in text
    assert "running on the host despite --container" not in text


# --- pinned BioContainers, one stage per image, without Wave ---------------------

_PINNED_DEREP = [("galah", []), ("drep", ["--virus"]), ("sourmash", [])]


@pytest.mark.parametrize(("tool", "extra"), _PINNED_DEREP, ids=[t for t, _ in _PINNED_DEREP])
def test_pinned_dereplicator_image(run_repgenr, clonal_wd, tool, extra) -> None:
    _genomes, wd = clonal_wd
    run_repgenr(*DOCKER, "dereplicate", "-wd", wd, "--tool", tool, "-t", "2", *extra)
    text = log_text(wd)
    assert f"quay.io/biocontainers/{tool}:" in text
    assert "running on the host despite --container" not in text


_PINNED_TREES = ["mashtree", "sourmash", "iqtree", "fasttree", "raxmlng"]


@pytest.mark.parametrize("builder", _PINNED_TREES)
def test_pinned_tree_builder_image(run_repgenr, clonal_wd, builder) -> None:
    """Alignment-free builders run on the genomes; the ML builders on a ska2
    alignment of every genome (RAxML-NG needs at least four sequences)."""
    _genomes, wd = clonal_wd
    run_repgenr(*DOCKER, "dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2")
    source = [] if builder in ("mashtree", "sourmash") else ["--msa-source", "snptype"]
    source += [] if builder in ("mashtree", "sourmash") else ["--snptyper", "ska2", "--all-genomes"]
    run_repgenr(
        *DOCKER, "phylo", "-wd", wd, "--treebuilder", builder, "--no-outgroup", "-t", "2", *source
    )
    text = log_text(wd)
    image = {"raxmlng": "raxml-ng"}.get(builder, builder)
    assert f"quay.io/biocontainers/{image}:" in text
    assert "running on the host despite --container" not in text
    assert (wd / "tree" / TREE_NWK).exists()


def test_pinned_snippy_and_gubbins_images(run_repgenr, clonal_wd) -> None:
    _genomes, wd = clonal_wd
    run_repgenr(*DOCKER, "dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2")
    run_repgenr(
        *DOCKER,
        "snptype",
        "-wd",
        wd,
        "--tool",
        "snippy",
        "--all-genomes",
        "--mask",
        "gubbins",
        "-t",
        "2",
    )
    text = log_text(wd)
    assert "quay.io/biocontainers/snippy:" in text
    assert "quay.io/biocontainers/gubbins:" in text
    assert "running on the host despite --container" not in text


# --- assemblers: the stage through each pinned image on simulated reads -----------


def _simulated_reads_workdir(tmp_path: Path) -> Path:
    import random
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from readsim import simulate_paired_reads

    from repgenr.core.contracts import READS_TSV, ReadRow, write_reads

    rng = random.Random(3)
    genome = tmp_path / "syn.fasta"
    genome.write_text(
        ">syn\n" + "".join(rng.choice("ACGT") for _ in range(80_000)) + "\n", encoding="utf-8"
    )
    r1, r2 = simulate_paired_reads(genome, tmp_path / "reads", coverage=40)
    wd = tmp_path / "wd"
    wd.mkdir()
    write_reads(
        wd / READS_TSV,
        [
            ReadRow(
                "SRRSYN",
                "SAM1",
                "PRJ1",
                "Synthetic organism",
                "1",
                "ILLUMINA",
                "Illumina MiSeq",
                "PAIRED",
                80_000 * 40,
                21_000,
                "Synfam",
                "Syngen",
                "syn",
                (str(r1), str(r2)),
                (),
                (r1.stat().st_size, r2.stat().st_size),
            )
        ],
    )
    return wd


@pytest.mark.parametrize("tool", ["skesa", "shovill"])
def test_pinned_assembler_image_recovers_a_synthetic_genome(run_repgenr, tmp_path, tool) -> None:
    from repgenr.core.contracts import read_assembly_stats

    wd = _simulated_reads_workdir(tmp_path)
    run_repgenr(
        *DOCKER,
        "assemble",
        "-wd",
        wd,
        "--assembler",
        tool,
        "-t",
        "4",
        "--jobs",
        "1",
        "--memory-gb",
        "6",
    )
    assert f"quay.io/biocontainers/{tool}:" in log_text(wd)
    stats = read_assembly_stats(wd / "assembly_stats.tsv")[0]
    assert stats.assembler == tool and stats.total_length > 70_000 and stats.n_contigs < 20
    assert (wd / "genomes" / "Synfam_Syngen_syn_SRRSYN.fasta").exists()


def test_pinned_flye_image_recovers_a_synthetic_genome_from_long_reads(run_repgenr, tmp_path):
    import random
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from readsim import simulate_long_reads

    from repgenr.core.contracts import READS_TSV, ReadRow, read_assembly_stats, write_reads

    rng = random.Random(5)
    genome = tmp_path / "syn.fasta"
    genome.write_text(
        ">syn\n" + "".join(rng.choice("ACGT") for _ in range(120_000)) + "\n", encoding="utf-8"
    )
    reads = simulate_long_reads(genome, tmp_path / "reads", coverage=40, error_rate=0.02)
    wd = tmp_path / "wd"
    wd.mkdir()
    write_reads(
        wd / READS_TSV,
        [
            ReadRow(
                "ONTSYN",
                "SAM1",
                "PRJ1",
                "Synthetic organism",
                "1",
                "OXFORD_NANOPORE",
                "GridION",
                "SINGLE",
                120_000 * 40,
                400,
                "Synfam",
                "Syngen",
                "syn",
                (str(reads),),
                (),
                (reads.stat().st_size,),
            )
        ],
    )
    run_repgenr(*DOCKER, "assemble", "-wd", wd, "--assembler", "flye", "-t", "4")
    assert "quay.io/biocontainers/flye:" in log_text(wd)
    stats = read_assembly_stats(wd / "assembly_stats.tsv")[0]
    assert stats.assembler == "flye" and stats.total_length > 100_000 and stats.n_contigs < 10
