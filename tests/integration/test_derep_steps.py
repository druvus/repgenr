"""Discrete dereplicate-chunk / dereplicate-merge steps (scatter-gather substrate)."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import CLUSTERS_TSV, GENOME_STATUS_TSV, read_clusters
from repgenr.core.plugins import ToolCapabilities
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepResult,
    registry,
)
from repgenr.stages.derep_steps import (
    ChunkParams,
    MergeParams,
    dereplicate_chunk,
    dereplicate_merge,
)
from repgenr.stages.dereplicate import DereplicateParams
from repgenr.stages.dereplicate import run as dereplicate_run

_LOG = __import__("logging").getLogger("test")


class _Halver(Dereplicator):
    """Collapses each adjacent pair to one rep, so every pass ~halves the set."""

    capabilities = ToolCapabilities(name="halver", supports_native_scaling=True)

    def preflight(self) -> dict[str, str]:
        return {"halver": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        genomes = list(genomes)
        reps = genomes[::2]
        clusters: dict[str, list[str]] = {}
        status: dict[str, str] = {}
        for i, rep in enumerate(reps):
            members = [genomes[2 * i + 1].name] if 2 * i + 1 < len(genomes) else []
            clusters[rep.name] = members
            status[rep.name] = STATUS_REPRESENTATIVE
            for m in members:
                status[m] = STATUS_CONTAINED
        return DerepResult(representatives=reps, clusters=clusters, genome_status=status)


@pytest.fixture
def reg():
    registry._load()
    registry._classes["halver"] = _Halver
    yield
    registry._classes.pop("halver", None)


def _make_genomes(gdir: Path, n: int) -> list[Path]:
    gdir.mkdir(parents=True, exist_ok=True)
    out = []
    for i in range(n):
        p = gdir / f"Fam_g_s_GCA_{i:06d}.1.fasta"
        p.write_text(">x\nACGT\n")
        out.append(p)
    return out


def test_chunk_writes_a_valid_contract(tmp_path: Path, reg) -> None:
    genomes = _make_genomes(tmp_path / "genomes", 4)
    out = tmp_path / "chunk0"
    res = dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes, out_dir=out), _LOG
    )
    # halver of 4 -> reps [g0, g2]
    assert {r.name for r in res.representatives} == {genomes[0].name, genomes[2].name}
    # contract files + representative FASTAs are present on disk
    assert (out / CLUSTERS_TSV).exists()
    assert (out / GENOME_STATUS_TSV).exists()
    rep_files = {p.name for p in (out / "representatives").iterdir()}
    assert rep_files == {genomes[0].name, genomes[2].name}
    clusters = read_clusters(out / CLUSTERS_TSV)
    assert clusters[genomes[0].name] == [genomes[1].name]
    # tool intermediates are trimmed; the staged representatives survive
    assert not (out / "scratch").exists()
    assert (out / "representatives" / genomes[0].name).read_text() == ">x\nACGT\n"


def test_merge_composes_membership_over_chunks(tmp_path: Path, reg) -> None:
    genomes = _make_genomes(tmp_path / "genomes", 8)
    # two chunks of four
    c0 = dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes[:4], out_dir=tmp_path / "c0"), _LOG
    )
    c1 = dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes[4:], out_dir=tmp_path / "c1"), _LOG
    )
    assert len(c0.representatives) == 2 and len(c1.representatives) == 2

    final = dereplicate_merge(
        MergeParams(
            tool="halver",
            chunk_dirs=[tmp_path / "c0", tmp_path / "c1"],
            out_dir=tmp_path / "merged",
        ),
        _LOG,
    )
    # union [g0,g2,g4,g6] halved -> [g0, g4]; composition pulls in every original.
    assert {r.name for r in final.representatives} == {genomes[0].name, genomes[4].name}
    clusters = read_clusters(tmp_path / "merged" / CLUSTERS_TSV)
    assert sorted(clusters[genomes[0].name]) == sorted(
        [genomes[1].name, genomes[2].name, genomes[3].name]
    )
    assert sorted(clusters[genomes[4].name]) == sorted(
        [genomes[5].name, genomes[6].name, genomes[7].name]
    )
    # every original genome accounted for exactly once (8 = 2 reps + 6 members)
    all_members = [m for ms in clusters.values() for m in ms]
    assert len(all_members) + len(clusters) == 8


def test_discrete_matches_in_process_chunked(workdir: Path, reg) -> None:
    """Scatter-gather steps reproduce the in-process two-stage membership."""
    genomes = _make_genomes(workdir / "genomes", 8)

    # Reference: the shared-workdir stage with the same chunk size.
    ctx = WorkdirContext(workdir, create=True)
    ref = dereplicate_run(ctx, DereplicateParams(tool="halver", process_size=4))

    # Discrete: two chunks then a merge, with the same final thresholds.
    dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes[:4], out_dir=workdir / "c0"), _LOG
    )
    dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes[4:], out_dir=workdir / "c1"), _LOG
    )
    final = dereplicate_merge(
        MergeParams(
            tool="halver",
            chunk_dirs=[workdir / "c0", workdir / "c1"],
            out_dir=workdir / "merged",
        ),
        _LOG,
    )

    def membership(result: DerepResult) -> dict[str, list[str]]:
        return {rep: sorted(members) for rep, members in result.clusters.items()}

    assert {r.name for r in final.representatives} == {r.name for r in ref.representatives}
    assert membership(final) == membership(ref)


def test_chunk_rejects_missing_genome(tmp_path: Path, reg) -> None:
    from repgenr.core.errors import WorkdirError

    bogus = [tmp_path / "nope.fasta"]
    with pytest.raises(WorkdirError):
        dereplicate_chunk(
            ChunkParams(tool="halver", genomes=bogus, out_dir=tmp_path / "x"), _LOG
        )
