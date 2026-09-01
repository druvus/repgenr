"""Taxonomy-aware reduction (--reduce species|genus) after ANI dereplication."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import GENOME_STATUS_TSV, read_genome_status
from repgenr.core.manifest import GenomeRecord
from repgenr.core.plugins import ToolCapabilities
from repgenr.dereplicators.base import (
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
