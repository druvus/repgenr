"""Taxonomy-aware reduction (--reduce species|genus) after ANI dereplication."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import GENOME_STATUS_TSV, read_genome_status
from repgenr.core.manifest import GenomeRecord
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

# (accession, filename, genus, species)
_GENOMES = [
    ("GCF_000001.1", "Fam_aaa_sp1_GCF_000001.1.fasta", "aaa", "aaa-sp1"),
    ("GCF_000002.1", "Fam_aaa_sp1_GCF_000002.1.fasta", "aaa", "aaa-sp1"),  # same species
    ("GCF_000003.1", "Fam_aaa_sp2_GCF_000003.1.fasta", "aaa", "aaa-sp2"),
    ("GCF_000004.1", "Fam_bbb_sp3_GCF_000004.1.fasta", "bbb", "bbb-sp3"),
]


class _NoRep(Dereplicator):
    """Keeps every input genome as its own representative (no ANI clustering)."""

    capabilities = ToolCapabilities(name="norep", supports_native_scaling=True)

    def preflight(self) -> dict[str, str]:
        return {"norep": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        genomes = list(genomes)
        return DerepResult(
            representatives=list(genomes),
            clusters={g.name: [] for g in genomes},
            genome_status={g.name: STATUS_REPRESENTATIVE for g in genomes},
        )


class _QcNoRep(_NoRep):
    """NoRep that rejects one genome on QC (as dRep's checkM filter does)."""

    capabilities = ToolCapabilities(name="qcnorep", supports_native_scaling=True)
    fail_names: frozenset[str] = frozenset()

    def preflight(self) -> dict[str, str]:
        return {"qcnorep": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        genomes = list(genomes)
        failed = [g for g in genomes if g.name in type(self).fail_names]
        kept = [g for g in genomes if g.name not in type(self).fail_names]
        result = super().dereplicate(kept, out_dir, params, logger)
        for g in failed:
            result.genome_status[g.name] = STATUS_FAIL_QC
        return result


@pytest.fixture
def taxo_workdir(workdir: Path):
    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    for _acc, fn, _g, _s in _GENOMES:
        (gdir / fn).write_text(">x\nACGT\n")
    registry._load()
    registry.register("norep", _NoRep, replace=True)
    registry.register("qcnorep", _QcNoRep, replace=True)
    ctx = WorkdirContext(workdir, create=True)
    ctx.manifest.upsert_many(
        [GenomeRecord(accession=a, filename=fn, genus=g, species=s) for a, fn, g, s in _GENOMES]
    )
    yield ctx
    registry._classes.pop("norep", None)
    registry._classes.pop("qcnorep", None)
    _QcNoRep.fail_names = frozenset()


def test_reduce_none_keeps_all(taxo_workdir) -> None:
    res = run(taxo_workdir, DereplicateParams(tool="norep", reduce="none"))
    assert len(res.representatives) == 4


def test_reduce_species(taxo_workdir) -> None:
    res = run(taxo_workdir, DereplicateParams(tool="norep", reduce="species"))
    # aaa-sp1 (2) collapses to 1; aaa-sp2; bbb-sp3 -> 3 representatives
    assert len(res.representatives) == 3
    # every original genome is accounted for (rep or contained)
    assert len(res.genome_status) == 4


def test_reduce_genus(taxo_workdir) -> None:
    res = run(taxo_workdir, DereplicateParams(tool="norep", reduce="genus"))
    # genus aaa (3) -> 1, genus bbb (1) -> 1  => 2 representatives
    assert len(res.representatives) == 2
    assert len(res.genome_status) == 4


class _FixedClusters(Dereplicator):
    """Returns a canned DerepResult regardless of input; lets a test control
    cluster sizes directly, independent of ANI/quality logic."""

    capabilities = ToolCapabilities(name="fixedclusters", supports_native_scaling=True)
    representatives: list[str] = []
    clusters: dict[str, list[str]] = {}

    def preflight(self) -> dict[str, str]:
        return {"fixedclusters": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        by_name = {g.name: g for g in genomes}
        reps = [by_name[n] for n in type(self).representatives]
        status = {n: STATUS_REPRESENTATIVE for n in type(self).representatives}
        for members in type(self).clusters.values():
            for m in members:
                status[m] = STATUS_CONTAINED
        return DerepResult(
            representatives=reps, clusters=dict(type(self).clusters), genome_status=status
        )


# Two representatives (A, B) of the same species: A has a bigger cluster but
# worse assembly quality; B has no cluster members but the better quality.
# (accession, filename, genus, species, completeness, contamination)
_KEEPER_GENOMES = [
    ("GCF_100001.1", "Fam_aaa_sp1_GCF_100001.1.fasta", "aaa", "aaa-sp1", 70.0, 5.0),  # A: rep
    ("GCF_100002.1", "Fam_aaa_sp1_GCF_100002.1.fasta", "aaa", "aaa-sp1", 99.0, 0.1),  # B: rep
    ("GCF_100003.1", "Fam_aaa_sp1_GCF_100003.1.fasta", "aaa", "aaa-sp1", None, None),  # M1
    ("GCF_100004.1", "Fam_aaa_sp1_GCF_100004.1.fasta", "aaa", "aaa-sp1", None, None),  # M2
    ("GCF_100005.1", "Fam_aaa_sp1_GCF_100005.1.fasta", "aaa", "aaa-sp1", None, None),  # M3
]
_A, _B, _M1, _M2, _M3 = (fn for _acc, fn, *_rest in _KEEPER_GENOMES)


@pytest.fixture
def keeper_workdir(workdir: Path):
    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    for _acc, fn, *_rest in _KEEPER_GENOMES:
        (gdir / fn).write_text(">x\nACGT\n")
    registry._load()
    registry.register("fixedclusters", _FixedClusters, replace=True)
    _FixedClusters.representatives = [_A, _B]
    _FixedClusters.clusters = {_A: [_M1, _M2, _M3], _B: []}
    ctx = WorkdirContext(workdir, create=True)
    ctx.manifest.upsert_many(
        [
            GenomeRecord(
                accession=a,
                filename=fn,
                genus=g,
                species=s,
                completeness=comp,
                contamination=contam,
            )
            for a, fn, g, s, comp, contam in _KEEPER_GENOMES
        ]
    )
    yield ctx
    registry._classes.pop("fixedclusters", None)


def test_reduce_keeper_quality_prefers_better_score_over_cluster_size(keeper_workdir) -> None:
    """Under keeper='quality', the taxonomy-reduce keeper follows assembly
    quality even when the rival representative has the larger ANI cluster."""
    res = run(
        keeper_workdir,
        DereplicateParams(tool="fixedclusters", reduce="species", keeper="quality"),
    )
    assert [r.name for r in res.representatives] == [_B]
    assert sorted(res.clusters[_B]) == sorted([_A, _M1, _M2, _M3])


def test_reduce_keeper_tool_prefers_larger_cluster(keeper_workdir) -> None:
    """Under keeper='tool', the taxonomy-reduce keeper ignores quality and
    keeps the largest-cluster rule (pre-existing behavior)."""
    res = run(
        keeper_workdir,
        DereplicateParams(tool="fixedclusters", reduce="species", keeper="tool"),
    )
    assert [r.name for r in res.representatives] == [_A]
    assert sorted(res.clusters[_A]) == sorted([_B, _M1, _M2, _M3])


def test_reduce_keeps_fail_qc_genomes(taxo_workdir) -> None:
    """A QC-rejected genome belongs to no cluster; reduction must still keep it."""
    bad = _GENOMES[3][1]
    _QcNoRep.fail_names = frozenset({bad})

    res = run(taxo_workdir, DereplicateParams(tool="qcnorep", reduce="species"))

    # aaa-sp1 (2) collapses to 1, aaa-sp2 stays => 2 representatives
    assert len(res.representatives) == 2
    assert res.genome_status[bad] == STATUS_FAIL_QC
    assert len(res.genome_status) == 4
    on_disk = read_genome_status(taxo_workdir.derep_dir / GENOME_STATUS_TSV)
    assert on_disk[bad] == STATUS_FAIL_QC
    assert len(on_disk) == 4
