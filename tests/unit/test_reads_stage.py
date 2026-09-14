"""The reads stage: select runs from ENA, label them, write reads.tsv."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from repgenr.core import ena
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import READS_TSV, read_reads
from repgenr.core.errors import UserInputError
from repgenr.stages import reads as reads_mod
from repgenr.stages.reads import ReadsParams, run

DATA = Path(__file__).parent / "data"


def _records(*names: str) -> list[dict]:
    out = []
    for name in names:
        out.extend(json.loads((DATA / name).read_text(encoding="utf-8")))
    return out


def _fake_entrez(taxids, logger):
    lineage = {
        "2097": ("Mycoplasmoidaceae", "Mycoplasmoides", "Mycoplasmoides genitalium"),
        "2094": ("Metamycoplasmataceae", "Mycoplasmopsis", "Mycoplasmopsis arginini"),
        "2102": ("Mycoplasmataceae", "Mycoplasma", "Mycoplasma mycoides"),
    }
    data = {}
    for taxid in taxids:
        fam, gen, sp = lineage[taxid]
        data[taxid] = {
            "taxid": taxid,
            "name": sp,
            "taxdata": {
                "family": {"taxid": "x", "name": fam, "level": "family"},
                "genus": {"taxid": "y", "name": gen, "level": "genus"},
                "species": {"taxid": taxid, "name": sp, "level": "species"},
            },
        }
    return data, [], {}


@pytest.fixture
def ena_fake(monkeypatch):
    """ENA answers from the frozen fixtures; Entrez from a canned lineage map."""
    seen: dict = {}

    def search_runs(query, **kw):
        seen["query"] = query
        if "tax_tree" in query:
            return _records("ena_read_run_taxon.json", "ena_read_run_mixed.json")
        return _records("ena_read_run_accessions.json")

    monkeypatch.setattr(ena, "search_runs", search_runs)
    monkeypatch.setattr(
        ena,
        "resolve_taxon",
        lambda name: ena.TaxonHit("2097", "Mycoplasmoides genitalium", "species"),
    )
    monkeypatch.setattr(reads_mod, "get_taxon_data_from_entrez", _fake_entrez)
    return seen


def test_taxon_query_writes_labelled_rows(workdir: Path, ena_fake) -> None:
    ctx = WorkdirContext(workdir, create=True)
    n = run(ctx, ReadsParams(target_species="Mycoplasma genitalium", one_per_sample=False))
    rows = read_reads(workdir / READS_TSV)
    assert n == len(rows) == 8  # 6 Illumina + 1 ONT + 1 PacBio without a mirror
    assert "tax_tree(2097)" in ena_fake["query"]
    ont = next(r for r in rows if r.run_accession == "SRR28800588")
    assert (ont.family, ont.genus, ont.species) == (
        "Mycoplasmoidaceae",
        "Mycoplasmoides",
        "genitalium",
    )
    record = ctx.config.stages["reads"]
    assert record.params["taxid"] == "2097" and record.params["selected"] == 8
    assert record.params["target_species"] == "Mycoplasma genitalium"


def test_rows_are_sorted_by_bases_and_capped(workdir: Path, ena_fake) -> None:
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, ReadsParams(target_species="x", one_per_sample=False, max_runs=3))
    rows = read_reads(workdir / READS_TSV)
    assert [r.bases for r in rows] == sorted((r.bases for r in rows), reverse=True)
    assert len(rows) == 3


def test_platform_filter_and_min_bases(workdir: Path, ena_fake) -> None:
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, ReadsParams(target_species="x", platform="ont", one_per_sample=False))
    assert [r.run_accession for r in read_reads(workdir / READS_TSV)] == ["SRR28800588"]
    run(ctx, ReadsParams(target_species="x", platform="illumina", min_bases=50_000_000))
    assert {r.platform for r in read_reads(workdir / READS_TSV)} == {"ILLUMINA"}
    assert all(r.bases >= 50_000_000 for r in read_reads(workdir / READS_TSV))


def test_one_per_sample_keeps_the_best_run_of_each_sample(workdir: Path, ena_fake) -> None:
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, ReadsParams(accessions=["PRJNA954307"]))
    rows = read_reads(workdir / READS_TSV)
    samples = [r.biosample for r in rows]
    assert len(samples) == len(set(samples))
    # SAMN34138608 has one run in the fixture; the largest run of the project leads
    assert rows[0].bases == max(r.bases for r in rows)
    assert 'study_accession="PRJNA954307"' in ena_fake["query"]


def test_accession_file_is_read_and_recorded(workdir: Path, ena_fake) -> None:
    listing = workdir.parent / "runs.txt"
    listing.write_text("SRR25474756\n# comment\n\nPRJNA954307\n", encoding="utf-8")
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, ReadsParams(accession_file=str(listing)))
    assert 'run_accession="SRR25474756"' in ena_fake["query"]
    assert ctx.config.stages["reads"].params["accession_file"] == str(listing)


def test_selection_is_required(workdir: Path, ena_fake) -> None:
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(UserInputError, match="taxon|accession"):
        run(ctx, ReadsParams())


def test_no_matching_runs_is_an_input_error(workdir: Path, monkeypatch, ena_fake) -> None:
    monkeypatch.setattr(ena, "search_runs", lambda query, **kw: [])
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(UserInputError, match="No sequencing runs"):
        run(ctx, ReadsParams(target_genus="Nothing"))
