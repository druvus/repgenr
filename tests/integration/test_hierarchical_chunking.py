"""Hierarchical (recursive) chunked dereplication."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import GENOME_STATUS_TSV, read_genome_status
from repgenr.core.plugins import ToolCapabilities
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_FAIL_QC,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepResult,
    registry,
)
from repgenr.stages.dereplicate import DereplicateParams, run


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


class _NoRep(Dereplicator):
    """Keeps every genome (no reduction) -- exercises the termination guard."""

    capabilities = ToolCapabilities(name="norep2", supports_native_scaling=True)

    def preflight(self) -> dict[str, str]:
        return {"norep2": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        genomes = list(genomes)
        return DerepResult(
            representatives=list(genomes),
            clusters={g.name: [] for g in genomes},
            genome_status={g.name: STATUS_REPRESENTATIVE for g in genomes},
        )


class _QcHalver(_Halver):
    """Halver that rejects the genomes in ``fail_names`` before clustering.

    Mirrors dRep's checkM filter: a rejected genome gets a ``fail_qc`` status and
    belongs to no cluster.
    """

    capabilities = ToolCapabilities(name="qchalver", supports_native_scaling=True)
    fail_names: frozenset[str] = frozenset()

    def preflight(self) -> dict[str, str]:
        return {"qchalver": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        genomes = list(genomes)
        failed = [g for g in genomes if g.name in type(self).fail_names]
        kept = [g for g in genomes if g.name not in type(self).fail_names]
        result = super().dereplicate(kept, out_dir, params, logger)
        for g in failed:
            result.genome_status[g.name] = STATUS_FAIL_QC
        return result


def _make_genomes(workdir: Path, n: int) -> None:
    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    for i in range(n):
        (gdir / f"Fam_g_s_GCA_{i:06d}.1.fasta").write_text(">x\nACGT\n")


def _genome_name(i: int) -> str:
    return f"Fam_g_s_GCA_{i:06d}.1.fasta"


@pytest.fixture
def reg():
    registry._load()
    registry.register("halver", _Halver, replace=True)
    registry.register("norep2", _NoRep, replace=True)
    registry.register("qchalver", _QcHalver, replace=True)
    yield
    registry._classes.pop("halver", None)
    registry._classes.pop("norep2", None)
    registry._classes.pop("qchalver", None)
    _QcHalver.fail_names = frozenset()


def test_recursive_reduction_accounts_all(workdir: Path, reg) -> None:
    _make_genomes(workdir, 8)
    ctx = WorkdirContext(workdir, create=True)
    # 8 genomes, size 2 -> multi-level reduce-tree (depth 0 and 1)
    res = run(ctx, DereplicateParams(tool="halver", process_size=2))
    assert len(res.representatives) < 8           # genuinely reduced
    assert len(res.genome_status) == 8            # every original genome accounted for


def test_fail_qc_genome_survives_chunked_composition(workdir: Path, reg) -> None:
    """A genome the tool rejects on QC keeps its status through the chunked path.

    It belongs to no cluster, so composing the two stages from clusters alone
    drops it and the completeness check then aborts the stage.
    """
    _make_genomes(workdir, 8)
    bad = _genome_name(7)
    _QcHalver.fail_names = frozenset({bad})
    ctx = WorkdirContext(workdir, create=True)

    res = run(ctx, DereplicateParams(tool="qchalver", process_size=4))

    assert res.genome_status[bad] == STATUS_FAIL_QC
    assert len(res.genome_status) == 8
    on_disk = read_genome_status(ctx.derep_dir / GENOME_STATUS_TSV)
    assert on_disk[bad] == STATUS_FAIL_QC
    assert len(on_disk) == 8
    # a fail_qc genome is never a representative
    assert bad not in {r.name for r in res.representatives}


def test_non_shrinking_terminates(workdir: Path, reg) -> None:
    _make_genomes(workdir, 6)
    ctx = WorkdirContext(workdir, create=True)
    # union never shrinks; the guard must fall back to a single pass, not recurse forever
    res = run(ctx, DereplicateParams(tool="norep2", process_size=2))
    assert len(res.representatives) == 6
    assert len(res.genome_status) == 6
