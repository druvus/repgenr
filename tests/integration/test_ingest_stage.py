"""ingest stage: populate a workdir from a local directory of genome FASTAs.

The stage is the offline counterpart of metadata + genome: it writes
``genomes/``, ``selection.tsv``, the manifest and (optionally) the outgroup
pair, so dereplicate -> phylo -> tree2tax can run without any download.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import SelectionRow, read_selection, write_selection
from repgenr.core.errors import UserInputError
from repgenr.stages.ingest import IngestParams, run

_SEQ = ">s\n" + "ACGT" * 10 + "\n"


def _source(tmp_path: Path, names: list[str]) -> Path:
    src = tmp_path / "src"
    src.mkdir()
    for name in names:
        (src / name).write_text(_SEQ)
    return src


def test_ingest_symlinks_genomes_and_parses_canonical_names(tmp_path: Path, workdir: Path) -> None:
    src = _source(
        tmp_path,
        ["Fam_Gen_sp1_GCA_000001.1.fasta", "Fam_Gen_sp2_GCA_000002.1.fasta", "._junk.fasta"],
    )
    ctx = WorkdirContext(workdir, create=True)

    n = run(ctx, IngestParams(genomes_dir=str(src)))

    assert n == 2
    linked = sorted(p.name for p in ctx.genomes_dir.iterdir())
    assert linked == ["Fam_Gen_sp1_GCA_000001.1.fasta", "Fam_Gen_sp2_GCA_000002.1.fasta"]
    assert all((ctx.genomes_dir / name).is_symlink() for name in linked)
    rows = {r.accession: r for r in read_selection(workdir / "selection.tsv")}
    assert rows["GCA_000002.1"].species == "sp2"
    assert rows["GCA_000002.1"].filename == "Fam_Gen_sp2_GCA_000002.1.fasta"
    assert not rows["GCA_000002.1"].is_outgroup
    manifest_rows = {g.accession: g for g in ctx.manifest.all_genomes(include_outgroup=True)}
    assert manifest_rows["GCA_000001.1"].genus == "Gen"
    assert manifest_rows["GCA_000001.1"].source == "local"
    assert ctx.config.stages["ingest"].completed
    assert not (workdir / "outgroup_accession.txt").exists()


def test_ingest_non_canonical_name_uses_stem_as_accession(tmp_path: Path, workdir: Path) -> None:
    src = _source(tmp_path, ["sample1.fa"])
    ctx = WorkdirContext(workdir, create=True)

    run(ctx, IngestParams(genomes_dir=str(src)))

    (row,) = read_selection(workdir / "selection.tsv")
    assert row.accession == "sample1"
    assert (row.family, row.genus, row.species) == ("", "", "")
    assert (ctx.genomes_dir / "sample1.fa").exists()


def test_ingest_copy_mode_writes_real_files(tmp_path: Path, workdir: Path) -> None:
    src = _source(tmp_path, ["a.fasta"])
    ctx = WorkdirContext(workdir, create=True)

    run(ctx, IngestParams(genomes_dir=str(src), copy=True))

    target = ctx.genomes_dir / "a.fasta"
    assert target.is_file() and not target.is_symlink()
    assert target.read_text() == _SEQ


def test_ingest_selection_tsv_drives_set_and_quality(tmp_path: Path, workdir: Path) -> None:
    src = _source(tmp_path, ["a.fasta", "b.fasta", "c.fasta"])
    selection = tmp_path / "sel.tsv"
    write_selection(
        selection,
        [
            SelectionRow("A1", "F", "G", "s1", False, "a.fasta", 99.0, 0.5),
            SelectionRow("B1", "F", "G", "s1", False, "b.fasta", 80.0, 3.0),
        ],
    )
    ctx = WorkdirContext(workdir, create=True)

    n = run(ctx, IngestParams(genomes_dir=str(src), selection=str(selection)))

    assert n == 2
    assert sorted(p.name for p in ctx.genomes_dir.iterdir()) == ["a.fasta", "b.fasta"]
    rows = {r.accession: r for r in read_selection(workdir / "selection.tsv")}
    assert rows["A1"].completeness == 99.0 and rows["B1"].contamination == 3.0
    manifest = {g.accession: g for g in ctx.manifest.all_genomes(include_outgroup=True)}
    assert manifest["A1"].completeness == 99.0
    assert manifest["A1"].species == "s1"


def test_ingest_selection_row_without_file_is_an_error(tmp_path: Path, workdir: Path) -> None:
    src = _source(tmp_path, ["a.fasta"])
    selection = tmp_path / "sel.tsv"
    write_selection(selection, [SelectionRow("Z9", "F", "G", "s", False, "zzz.fasta")])
    ctx = WorkdirContext(workdir, create=True)

    with pytest.raises(UserInputError, match="zzz.fasta"):
        run(ctx, IngestParams(genomes_dir=str(src), selection=str(selection)))


def test_ingest_outgroup_from_selection_row(tmp_path: Path, workdir: Path) -> None:
    src = _source(tmp_path, ["a.fasta", "og.fasta"])
    selection = tmp_path / "sel.tsv"
    write_selection(
        selection,
        [
            SelectionRow("A1", "F", "G", "s", False, "a.fasta"),
            SelectionRow("OG1", "F", "H", "t", True, "og.fasta"),
        ],
    )
    ctx = WorkdirContext(workdir, create=True)

    run(ctx, IngestParams(genomes_dir=str(src), selection=str(selection)))

    assert [p.name for p in ctx.genomes_dir.iterdir()] == ["a.fasta"]
    assert (ctx.outgroup_dir / "og.fasta").exists()
    assert (workdir / "outgroup_accession.txt").read_text().strip() == "OG1"
    og = next(g for g in ctx.manifest.all_genomes(include_outgroup=True) if g.is_outgroup)
    assert og.accession == "OG1"


def test_ingest_outgroup_flag_names_a_source_genome(tmp_path: Path, workdir: Path) -> None:
    src = _source(tmp_path, ["Fam_Gen_sp1_GCA_000001.1.fasta", "Fam_Out_sp9_GCA_000009.1.fasta"])
    ctx = WorkdirContext(workdir, create=True)

    n = run(ctx, IngestParams(genomes_dir=str(src), outgroup="GCA_000009.1"))

    assert n == 1
    assert [p.name for p in ctx.genomes_dir.iterdir()] == ["Fam_Gen_sp1_GCA_000001.1.fasta"]
    assert (ctx.outgroup_dir / "Fam_Out_sp9_GCA_000009.1.fasta").exists()
    assert (workdir / "outgroup_accession.txt").read_text().strip() == "GCA_000009.1"


def test_ingest_outgroup_flag_accepts_external_fasta_path(tmp_path: Path, workdir: Path) -> None:
    src = _source(tmp_path, ["a.fasta"])
    external = tmp_path / "Out_grp_sp_X1.fasta"
    external.write_text(_SEQ)
    ctx = WorkdirContext(workdir, create=True)

    run(ctx, IngestParams(genomes_dir=str(src), outgroup=str(external)))

    assert (ctx.outgroup_dir / "Out_grp_sp_X1.fasta").exists()
    assert (workdir / "outgroup_accession.txt").read_text().strip() == "X1"


def test_ingest_unknown_outgroup_is_an_error(tmp_path: Path, workdir: Path) -> None:
    src = _source(tmp_path, ["a.fasta"])
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(UserInputError, match="nope"):
        run(ctx, IngestParams(genomes_dir=str(src), outgroup="nope"))


def test_ingest_empty_source_is_an_error(tmp_path: Path, workdir: Path) -> None:
    src = tmp_path / "empty"
    src.mkdir()
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(UserInputError, match="No genome FASTA"):
        run(ctx, IngestParams(genomes_dir=str(src)))


def test_ingest_rerun_prunes_stale_genomes_and_outgroup(tmp_path: Path, workdir: Path) -> None:
    src = _source(tmp_path, ["a.fasta", "b.fasta"])
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, IngestParams(genomes_dir=str(src), outgroup="b"))
    assert (ctx.outgroup_dir / "b.fasta").exists()

    os.remove(src / "b.fasta")
    (src / "c.fasta").write_text(_SEQ)
    run(ctx, IngestParams(genomes_dir=str(src)))

    assert sorted(p.name for p in ctx.genomes_dir.iterdir()) == ["a.fasta", "c.fasta"]
    assert not (workdir / "outgroup_accession.txt").exists()
    assert not any(ctx.outgroup_dir.iterdir()) if ctx.outgroup_dir.exists() else True
    accessions = {g.accession for g in ctx.manifest.all_genomes(include_outgroup=True)}
    assert accessions == {"a", "c"}
