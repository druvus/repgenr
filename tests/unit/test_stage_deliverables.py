"""STAGE_DELIVERABLES: the outputs whose absence makes a completed stage rerun."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from repgenr.cli import base as cli


def _layout(workdir: Path) -> SimpleNamespace:
    return SimpleNamespace(workdir=workdir, derep_dir=workdir / "derep")


def test_every_stage_with_inputs_declares_deliverables() -> None:
    """Each workdir stage lists the outputs whose absence forces a rerun."""
    assert set(cli.STAGE_INPUTS) <= set(cli.STAGE_DELIVERABLES)


def test_derep_unpack_checks_unpacked_dir_with_representatives(tmp_path: Path) -> None:
    ctx = _layout(tmp_path)
    params = SimpleNamespace(no_representant=False)
    assert cli.missing_deliverables(ctx, "derep_unpack", params) == [tmp_path / "derep/unpacked"]
    (tmp_path / "derep" / "unpacked").mkdir(parents=True)
    # An empty unpacked/ still counts as missing when representatives are unpacked.
    assert cli.missing_deliverables(ctx, "derep_unpack", params) == [tmp_path / "derep/unpacked"]


def test_derep_unpack_without_representant_accepts_empty_unpacked_dir(tmp_path: Path) -> None:
    """--no-representant with only singleton clusters leaves unpacked/ empty;
    that must not rerun the stage on every call."""
    (tmp_path / "derep" / "unpacked").mkdir(parents=True)
    params = SimpleNamespace(no_representant=True)
    assert cli.missing_deliverables(_layout(tmp_path), "derep_unpack", params) == []


def _genome_layout(workdir: Path) -> SimpleNamespace:
    return SimpleNamespace(
        workdir=workdir, genomes_dir=workdir / "genomes", outgroup_dir=workdir / "outgroup"
    )


def test_genome_deliverables_list_every_promised_file(tmp_path: Path) -> None:
    """A genome or outgroup deleted by hand reruns the stage (it fetches only
    the absent files); an accession NCBI returned nothing for is excused."""
    from repgenr.core.contracts import SelectionRow, write_selection
    from repgenr.core.manifest import MANIFEST_FILENAME

    ctx = _genome_layout(tmp_path)
    rows = [
        SelectionRow("GCF_1.1", "F", "G", "s", False, "F_G_s_GCF_1.1.fasta"),
        SelectionRow("GCF_2.1", "F", "G", "s", False, "F_G_s_GCF_2.1.fasta"),
        SelectionRow("GCF_3.1", "F", "G", "s", False, "F_G_s_GCF_3.1.fasta"),
        SelectionRow("GCF_9.1", "F", "H", "t", True, "F_H_t_GCF_9.1.fasta"),
    ]
    write_selection(tmp_path / "selection.tsv", rows)
    (tmp_path / MANIFEST_FILENAME).write_bytes(b"")
    (tmp_path / "missing_accessions.txt").write_text("GCF_3.1\n", encoding="utf-8")
    ctx.genomes_dir.mkdir()
    ctx.outgroup_dir.mkdir()
    (ctx.genomes_dir / rows[0].filename).write_text(">a\nACGT\n", encoding="utf-8")
    params = SimpleNamespace()
    assert cli.missing_deliverables(ctx, "genome", params) == [
        ctx.genomes_dir / rows[1].filename,
        ctx.outgroup_dir / rows[3].filename,
    ]
    (ctx.genomes_dir / rows[1].filename).write_text(">b\nACGT\n", encoding="utf-8")
    (ctx.outgroup_dir / rows[3].filename).write_text(">o\nACGT\n", encoding="utf-8")
    assert cli.missing_deliverables(ctx, "genome", params) == []
