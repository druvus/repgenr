"""End-to-end dereplicate stage test using an in-process fake adapter."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import (
    CLUSTERS_TSV,
    GENOME_STATUS_TSV,
    read_clusters,
)
from repgenr.core.plugins import ToolCapabilities
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepResult,
    registry,
)
from repgenr.stages.dereplicate import DereplicateParams, run


class _FakeDereplicator(Dereplicator):
    capabilities = ToolCapabilities(name="fake", supports_native_scaling=True)

    def preflight(self) -> dict[str, str]:
        return {"fake": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:
        genomes = list(genomes)
        rep = genomes[0]
        members = [g.name for g in genomes[1:]]
        return DerepResult(
            representatives=[rep],
            clusters={rep.name: members},
            genome_status={
                rep.name: STATUS_REPRESENTATIVE,
                **{m: STATUS_CONTAINED for m in members},
            },
        )


class _NonScalingDereplicator(_FakeDereplicator):
    """Same behaviour but flagged as not scaling natively, to exercise chunking."""

    capabilities = ToolCapabilities(name="chunky", supports_native_scaling=False)


class _RecordingDereplicator(_FakeDereplicator):
    """Native-scaling adapter that records the (size, secondary_ani) of each call."""

    capabilities = ToolCapabilities(name="recording", supports_native_scaling=True)
    calls: list[tuple[int, float]] = []

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:
        genomes = list(genomes)
        type(self).calls.append((len(genomes), params.secondary_ani))
        return super().dereplicate(genomes, out_dir, params, logger)


@pytest.fixture
def fake_tool() -> None:
    registry._load()
    registry._classes["fake"] = _FakeDereplicator
    registry._classes["chunky"] = _NonScalingDereplicator
    registry._classes["recording"] = _RecordingDereplicator
    _RecordingDereplicator.calls = []
    yield
    registry._classes.pop("fake", None)
    registry._classes.pop("chunky", None)
    registry._classes.pop("recording", None)


def test_dereplicate_writes_contract(workdir: Path, genome_files, fake_tool) -> None:
    ctx = WorkdirContext(workdir, create=True)
    result = run(ctx, DereplicateParams(tool="fake"))

    assert len(result.representatives) == 1

    rep_dir = ctx.representatives_dir
    assert rep_dir.is_dir()
    assert [p.name for p in rep_dir.iterdir()] == [genome_files[0].name]

    clusters = read_clusters(ctx.derep_dir / CLUSTERS_TSV)
    assert clusters[genome_files[0].name] == [genome_files[1].name, genome_files[2].name]

    assert (ctx.derep_dir / GENOME_STATUS_TSV).exists()
    assert (workdir / "repgenr.yaml").exists()
    assert ctx.config.stages["dereplicate"].tool == "fake"


def test_chunking_composes_membership(workdir: Path, genome_files, fake_tool) -> None:
    # process_size=2 with 3 genomes -> the lone last chunk is merged, so a
    # single chunk runs. Use a larger set to force a real two-stage pass.
    gdir = workdir / "genomes"
    extra = []
    for i in range(4, 9):
        name = f"Francisellaceae_francisella_tularensis_GCA_00000{i}.fasta"
        (gdir / name).write_text(f">s{i}\n{'ACGT' * 10}\n")
        extra.append(name)

    ctx = WorkdirContext(workdir, create=True)
    result = run(ctx, DereplicateParams(tool="chunky", process_size=2))
    # every original genome must be accounted for as rep or contained
    all_names = {p.name for p in gdir.iterdir()}
    accounted = set(result.genome_status)
    assert all_names <= accounted


def test_native_scaling_single_pass_by_default(workdir: Path, genome_files, fake_tool) -> None:
    # No process_size -> native-scaling tool runs in a single pass (one adapter call).
    ctx = WorkdirContext(workdir, create=True)
    result = run(ctx, DereplicateParams(tool="recording"))
    assert len(result.representatives) == 1
    assert len(_RecordingDereplicator.calls) == 1  # one pass over all genomes


def test_native_scaling_can_be_chunked(workdir: Path, genome_files, fake_tool) -> None:
    # With process_size set and exceeded, even a native-scaling tool is chunked
    # (escape hatch for very large sets). 6 genomes / size 2 -> 3 stage-1 chunks
    # + 1 stage-2 pass = 4 adapter calls.
    gdir = workdir / "genomes"
    for i in range(4, 7):
        (gdir / f"Francisellaceae_f_t_GCA_00000{i}.fasta").write_text(f">s{i}\n{'ACGT' * 10}\n")
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, DereplicateParams(tool="recording", process_size=2))
    calls = _RecordingDereplicator.calls
    assert len(calls) == 4  # 3 stage-1 chunks + 1 stage-2
    assert sorted(n for n, _ in calls[:-1]) == [2, 2, 2]  # each stage-1 chunk has 2 genomes
    assert calls[-1][0] == 3  # stage-2 over the 3 chunk representatives


def test_stage1_uses_pre_thresholds(workdir: Path, genome_files, fake_tool) -> None:
    # 3 genomes, process_size=2 -> trailing singleton merges into one chunk, so
    # only stage 1 runs; bump to >=4 genomes to force a real two-stage pass.
    gdir = workdir / "genomes"
    for i in range(4, 7):
        (gdir / f"Francisellaceae_f_t_GCA_00000{i}.fasta").write_text(f">s{i}\n{'ACGT' * 10}\n")

    ctx = WorkdirContext(workdir, create=True)
    run(ctx, DereplicateParams(
        tool="recording", process_size=2,
        secondary_ani=0.99, pre_secondary_ani=0.95,
    ))
    calls = _RecordingDereplicator.calls
    # stage-1 chunk calls use the looser pre threshold; the final stage-2 call
    # (on the union of stage-1 reps) uses the main threshold.
    stage1 = [s for s in calls[:-1]]
    assert stage1 and all(sec == 0.95 for _, sec in stage1)
    assert calls[-1][1] == 0.99
