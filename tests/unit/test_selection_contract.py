"""Selection contract (metadata -> genome data-channel hand-off) and genome-fetch."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.contracts import (
    SelectionRow,
    genome_filename,
    read_selection,
    write_selection,
)
from repgenr.stages import genome as genome_stage
from repgenr.stages.genome_steps import GenomeFetchParams, genome_fetch

_LOG = logging.getLogger("test")


def test_genome_filename_is_canonical() -> None:
    assert (
        genome_filename("Fam", "Gen", "spec", "GCF_000001.1")
        == "Fam_Gen_spec_GCF_000001.1.fasta"
    )


def test_selection_round_trip(tmp_path: Path) -> None:
    rows = [
        SelectionRow("GCF_1.1", "Fam", "Gen", "sp1", False, "Fam_Gen_sp1_GCF_1.1.fasta"),
        SelectionRow("GCF_2.1", "Fam", "Gen", "sp2", True, "Fam_Gen_sp2_GCF_2.1.fasta"),
    ]
    path = tmp_path / "selection.tsv"
    write_selection(path, rows)
    back = read_selection(path)
    assert back == rows
    # the outgroup flag survives the round trip
    assert [r.is_outgroup for r in back] == [False, True]


def test_genome_fetch_downloads_selected_and_outgroup(monkeypatch, tmp_path: Path) -> None:
    calls: list[tuple[list[str], Path]] = []

    def fake_download(accessions, filenames, dest_dir, scratch_dir, logger, keep_files=False):
        dest_dir.mkdir(parents=True, exist_ok=True)
        for acc in accessions:
            (dest_dir / filenames[acc]).write_text(">x\nACGT\n")
        calls.append((list(accessions), dest_dir))

    # genome_steps imports the symbol into its namespace, so patch it there.
    monkeypatch.setattr("repgenr.stages.genome_steps.download_accessions", fake_download)

    rows = [
        SelectionRow("GCF_1.1", "Fam", "Gen", "sp1", False, "Fam_Gen_sp1_GCF_1.1.fasta"),
        SelectionRow("GCF_2.1", "Fam", "Gen", "sp2", False, "Fam_Gen_sp2_GCF_2.1.fasta"),
        SelectionRow("GCF_9.1", "Fam", "Out", "grp", True, "Fam_Out_grp_GCF_9.1.fasta"),
    ]
    selection = tmp_path / "selection.tsv"
    write_selection(selection, rows)

    out = tmp_path / "out"
    n = genome_fetch(GenomeFetchParams(selection_tsv=selection, out_dir=out), _LOG)

    assert n == 2  # selected count (outgroup excluded)
    genomes = sorted(p.name for p in (out / "genomes").iterdir())
    assert genomes == ["Fam_Gen_sp1_GCF_1.1.fasta", "Fam_Gen_sp2_GCF_2.1.fasta"]
    assert (out / "outgroup" / "Fam_Out_grp_GCF_9.1.fasta").exists()
    # selected and outgroup went to separate destination directories
    dests = {dest.name for _, dest in calls}
    assert dests == {"genomes", "outgroup"}


def test_genome_fetch_rejects_missing_selection(tmp_path: Path) -> None:
    from repgenr.core.errors import WorkdirError

    with pytest.raises(WorkdirError):
        genome_fetch(
            GenomeFetchParams(selection_tsv=tmp_path / "nope.tsv", out_dir=tmp_path / "o"), _LOG
        )


def test_output_name_matches_canonical(tmp_path: Path) -> None:
    # the legacy genome stage and the canonical helper agree on filenames
    from repgenr.core.manifest import GenomeRecord

    rec = GenomeRecord(
        accession="GCF_5.1", source="gtdb", is_outgroup=False,
        family="Fam", genus="Gen", species="sp",
    )
    assert genome_stage._output_name(rec) == genome_filename("Fam", "Gen", "sp", "GCF_5.1")
