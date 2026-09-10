"""`repgenr ingest` flags as a real process (no external tool needed)."""

from __future__ import annotations

import pytest

from repgenr.core.contracts import SELECTION_TSV, read_selection

pytestmark = pytest.mark.live


def test_outgroup_and_copy(synthetic_set, ingested_workdir) -> None:
    genomes = synthetic_set("balanced", n=4, length=20_000)
    files = sorted(genomes.glob("*.fasta"))
    og_accession = files[-1].stem.rsplit("_", 1)[-1]
    wd = ingested_workdir(genomes, outgroup=og_accession, copy=True)

    staged = sorted(p.name for p in (wd / "genomes").iterdir())
    assert staged == [f.name for f in files[:-1]], "the outgroup is kept out of the ingroup"
    assert all(not (wd / "genomes" / n).is_symlink() for n in staged), "--copy makes real files"
    assert (wd / "outgroup" / files[-1].name).is_file()
    assert (wd / "outgroup_accession.txt").read_text(encoding="utf-8").strip() == og_accession
    rows = read_selection(wd / SELECTION_TSV)
    assert [r.accession for r in rows if r.is_outgroup] == [og_accession]


def test_selection_table_drives_taxonomy_and_subset(
    synthetic_set, selection_for, ingested_workdir
) -> None:
    genomes = synthetic_set("balanced", n=4, length=20_000)
    files = sorted(genomes.glob("*.fasta"))
    species = {f.name: "picked" for f in files[:2]}
    sel = selection_for(genomes, species=species, path=genomes.parent / "sel.tsv")
    # Keep only the first two rows: the selection, not the directory, decides.
    lines = sel.read_text(encoding="utf-8").splitlines()
    sel.write_text("\n".join(lines[:3]) + "\n", encoding="utf-8")
    wd = ingested_workdir(genomes, selection=sel)
    assert sorted(p.name for p in (wd / "genomes").iterdir()) == [f.name for f in files[:2]]
    assert {r.species for r in read_selection(wd / SELECTION_TSV)} == {"picked"}
