"""Live dereplication on a seeded clonal set: adapters and every flag.

The clonal set (n=12, 50 kb, clone fraction 0.5) has one near-identical block
of six and two balanced clusters of three at ~99.6% intra- and ~96%
inter-cluster ANI, so the default thresholds recover the generator's
partition and moving --secondary-ani across those bands changes the count.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from live_helpers import derep_partition, log_text, representatives, truth_of, truth_partition

from repgenr.core.config import Config

pytestmark = pytest.mark.live

ADAPTERS = [
    pytest.param("skder", marks=pytest.mark.requires_binary("skder")),
    pytest.param("sourmash", marks=pytest.mark.requires_binary("sourmash")),
    pytest.param("galah", marks=pytest.mark.requires_binary("galah")),
]


@pytest.fixture
def clonal(synthetic_set, ingested_workdir):
    genomes = synthetic_set("clonal", n=12, length=50_000)
    return genomes, ingested_workdir(genomes)


@pytest.mark.parametrize("tool", ADAPTERS)
def test_adapter_recovers_the_synthetic_partition(run_repgenr, clonal, tool: str) -> None:
    genomes, wd = clonal
    run_repgenr("dereplicate", "-wd", wd, "--tool", tool, "-t", "2")
    assert derep_partition(wd) == truth_partition(truth_of(genomes))
    rec = Config.load(wd).stages["dereplicate"]
    assert rec.tool == tool and rec.completed


@pytest.mark.requires_binary("sourmash")
def test_secondary_ani_sweep_changes_representative_count(run_repgenr, clonal) -> None:
    _, wd = clonal
    counts = {}
    for sani in ("0.95", "0.999"):
        run_repgenr("dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2", "-sani", sani)
        counts[sani] = len(representatives(wd))
    assert counts["0.95"] < counts["0.999"], counts
    assert counts["0.95"] == 1, "at 95% every cluster founder merges"


@pytest.mark.requires_binary("sourmash")
def test_process_size_runs_the_chunked_path(run_repgenr, clonal) -> None:
    genomes, wd = clonal
    run_repgenr(
        "dereplicate",
        "-wd",
        wd,
        "--tool",
        "sourmash",
        "-t",
        "2",
        "--process-size",
        "4",
        "--num-processes",
        "2",
        "-pani",
        "0.85",
        "-af",
        "0.4",
        "--pre-primary-ani",
        "0.85",
        "--pre-secondary-ani",
        "0.98",
    )
    chunk_dirs = sorted(p.name for p in (wd / "scratch" / "dereplicate" / "level0").glob("chunk*"))
    assert chunk_dirs == ["chunk0", "chunk1", "chunk2"], "12 genomes at size 4 give three chunks"
    assert "Chunking 12 genomes at size 4" in log_text(wd)
    assert derep_partition(wd) == truth_partition(truth_of(genomes))
    params = Config.load(wd).stages["dereplicate"].params
    assert (params["pre_primary_ani"], params["pre_secondary_ani"]) == (0.85, 0.98)
    assert (params["primary_ani"], params["aligned_fraction"]) == (0.85, 0.4)
    assert params["num_processes"] == 2
    assert "3 chunks, 2 parallel worker(s) at 1 threads each" in log_text(wd), (
        "-t 2 split across --num-processes 2"
    )


@pytest.mark.requires_binary("sourmash")
def test_target_reps_lands_on_the_requested_count(run_repgenr, clonal) -> None:
    _, wd = clonal
    run_repgenr("dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2", "--target-reps", "3")
    assert len(representatives(wd)) == 3
    assert "target-reps: chose secondary-ani=" in log_text(wd)


@pytest.mark.requires_binary("sourmash")
def test_reduce_species_keeps_one_representative_per_species(
    run_repgenr, synthetic_set, selection_for, ingested_workdir
) -> None:
    genomes = synthetic_set("clonal", n=12, length=50_000)
    truth = truth_of(genomes)
    species = {name: f"species_{cluster}" for name, cluster in truth["clusters"].items()}
    wd = ingested_workdir(genomes, selection=selection_for(genomes, species=species))
    # 0.999 splits the two balanced clusters into singletons; the reduction
    # then folds each species back to one representative.
    run_repgenr(
        "dereplicate",
        "-wd",
        wd,
        "--tool",
        "sourmash",
        "-t",
        "2",
        "-sani",
        "0.999",
        "--reduce",
        "species",
    )
    assert len(representatives(wd)) == 3
    assert "Taxonomy reduction (one per species)" in log_text(wd)


@pytest.mark.requires_binary("sourmash")
def test_keeper_quality_promotes_the_best_scored_member(
    run_repgenr, synthetic_set, selection_for, ingested_workdir
) -> None:
    genomes = synthetic_set("clonal", n=12, length=50_000)
    truth = truth_of(genomes)
    clone = sorted(n for n, c in truth["clusters"].items() if c == "clone")
    best = clone[-1]  # the last file in sort order is never the tool's own pick
    quality = {name: (90.0, 2.0) for name in truth["clusters"]}
    quality[best] = (99.9, 0.1)
    wd = ingested_workdir(genomes, selection=selection_for(genomes, quality=quality))

    run_repgenr("dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2", "--keeper", "quality")
    assert best in representatives(wd)
    run_repgenr("dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2", "--keeper", "tool")
    assert best not in representatives(wd), "the adapter's own pick is not the best-scored member"
    assert Config.load(wd).stages["dereplicate"].params["keeper"] == "tool"


@pytest.mark.requires_binary("sourmash")
def test_tool_arg_reaches_the_tool_command_line(run_repgenr, clonal) -> None:
    _, wd = clonal
    run_repgenr(
        "--verbose",
        "dereplicate",
        "-wd",
        wd,
        "--tool",
        "sourmash",
        "-t",
        "2",
        "--tool-arg",
        "ksize=21",
    )
    assert "-k 21" in log_text(wd), "sourmash pairwise must be called with the requested k"


@pytest.mark.requires_binary("skder")
def test_virus_flag_on_auto_tool_is_reported_when_ignored(run_repgenr, clonal) -> None:
    _, wd = clonal
    run_repgenr("dereplicate", "-wd", wd, "--tool", "auto", "--virus", "-t", "2")
    text = log_text(wd)
    assert "Auto-selected dereplicator" in text
    tool = Config.load(wd).stages["dereplicate"].tool
    from repgenr.dereplicators.base import registry

    if "virus" not in registry.get(tool).capabilities.accepted_extras:
        assert "ignores extra parameter(s): virus" in text


@pytest.mark.requires_binary("sourmash")
def test_allow_incomplete_gates_a_missing_genome(run_repgenr, clonal) -> None:
    _, wd = clonal
    victim = next(iter(sorted((wd / "genomes").glob("*.fasta"))))
    victim.unlink()
    refused = run_repgenr("dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2", check=False)
    assert refused.returncode != 0
    run_repgenr("dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2", "--allow-incomplete")
    assert len({m for group in derep_partition(wd) for m in group}) == 11


def _rep_files(wd: Path) -> set[str]:
    return {p.name for p in (wd / "derep" / "representatives").iterdir()}


@pytest.mark.requires_binary("sourmash")
def test_representatives_dir_matches_clusters(run_repgenr, clonal) -> None:
    _, wd = clonal
    run_repgenr("dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2")
    assert _rep_files(wd) == representatives(wd)
