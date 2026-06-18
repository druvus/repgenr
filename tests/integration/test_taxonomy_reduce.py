"""Taxonomy-aware reduction (--reduce species|genus) after ANI dereplication."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.manifest import GenomeRecord
from repgenr.core.plugins import ToolCapabilities
from repgenr.dereplicators.base import (
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


@pytest.fixture
def taxo_workdir(workdir: Path):
    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    for _acc, fn, _g, _s in _GENOMES:
        (gdir / fn).write_text(">x\nACGT\n")
    registry._load()
    registry._classes["norep"] = _NoRep
    ctx = WorkdirContext(workdir, create=True)
    ctx.manifest.upsert_many(
        [GenomeRecord(accession=a, filename=fn, genus=g, species=s) for a, fn, g, s in _GENOMES]
    )
    yield ctx
    registry._classes.pop("norep", None)


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
