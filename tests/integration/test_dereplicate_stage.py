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
    registry.register("fake", _FakeDereplicator, replace=True)
    registry.register("chunky", _NonScalingDereplicator, replace=True)
    registry.register("recording", _RecordingDereplicator, replace=True)
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


def test_keeper_quality_promotes_best_scoring_member(
    workdir: Path, genome_files, fake_tool
) -> None:
    # _FakeDereplicator clusters all three genomes under genome_files[0]. Seed the
    # manifest with quality so the tool's pick (rep) is the worst scorer, and
    # genome_files[2] the best -- the keeper should promote it to representative.
    ctx = WorkdirContext(workdir, create=True)
    from repgenr.core.contracts import accession_from_filename
    from repgenr.core.manifest import GenomeRecord

    quality = {
        genome_files[0].name: (80.0, 5.0),   # score 55.0 -- worst
        genome_files[1].name: (95.0, 1.0),   # score 90.0
        genome_files[2].name: (99.0, 0.2),   # score 98.0 -- best
    }
    for filename, (completeness, contamination) in quality.items():
        ctx.manifest.upsert(GenomeRecord(
            accession=accession_from_filename(filename), filename=filename,
            completeness=completeness, contamination=contamination,
        ))

    result = run(ctx, DereplicateParams(tool="fake", keeper="quality"))

    assert len(result.representatives) == 1
    assert result.representatives[0].name == genome_files[2].name

    rep_dir = ctx.representatives_dir
    assert [p.name for p in rep_dir.iterdir()] == [genome_files[2].name]

    assert ctx.config.stages["dereplicate"].params["keeper"] == "quality"
    assert ctx.config.stages["dereplicate"].params["keeper_swaps"] == 1


def test_keeper_tool_keeps_adapter_pick(workdir: Path, genome_files, fake_tool) -> None:
    ctx = WorkdirContext(workdir, create=True)
    from repgenr.core.contracts import accession_from_filename
    from repgenr.core.manifest import GenomeRecord

    quality = {
        genome_files[0].name: (80.0, 5.0),
        genome_files[1].name: (95.0, 1.0),
        genome_files[2].name: (99.0, 0.2),
    }
    for filename, (completeness, contamination) in quality.items():
        ctx.manifest.upsert(GenomeRecord(
            accession=accession_from_filename(filename), filename=filename,
            completeness=completeness, contamination=contamination,
        ))

    result = run(ctx, DereplicateParams(tool="fake", keeper="tool"))

    assert result.representatives[0].name == genome_files[0].name
    assert ctx.config.stages["dereplicate"].params["keeper"] == "tool"
    assert ctx.config.stages["dereplicate"].params["keeper_swaps"] == 0


def test_keeper_quality_without_manifest_quality_warns_and_records_fallback(
    workdir: Path, genome_files, fake_tool, caplog
) -> None:
    # No completeness/contamination in the manifest (e.g. an API selection that
    # could not fetch quality): the stage must say so loudly and record that the
    # adapter's own picks were kept, so "quality, 0 swaps" is not ambiguous.
    import logging

    ctx = WorkdirContext(workdir, create=True)
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.WARNING):
        run(ctx, DereplicateParams(tool="fake", keeper="quality"))

    params = ctx.config.stages["dereplicate"].params
    assert params["keeper"] == "quality"
    assert params["keeper_effective"] == "tool"
    warnings = [r.message for r in caplog.records if r.levelname == "WARNING"]
    assert any("quality" in m.lower() for m in warnings)


def test_keeper_quality_with_manifest_quality_records_quality(
    workdir: Path, genome_files, fake_tool
) -> None:
    ctx = WorkdirContext(workdir, create=True)
    from repgenr.core.contracts import accession_from_filename
    from repgenr.core.manifest import GenomeRecord

    for f in genome_files:
        ctx.manifest.upsert(GenomeRecord(
            accession=accession_from_filename(f.name), filename=f.name,
            completeness=99.0, contamination=0.5,
        ))
    run(ctx, DereplicateParams(tool="fake", keeper="quality"))
    assert ctx.config.stages["dereplicate"].params["keeper_effective"] == "quality"


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
