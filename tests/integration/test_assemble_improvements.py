"""The assemble stage: stored CheckM2 scores, the genus-rename flag, the
classifier memory budget, single-file PAIRED planning and the outputs of a
stage in which every run was excused."""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import (
    ASSEMBLY_STATS_TSV,
    EXCUSED_RUNS_TSV,
    READS_TSV,
    SELECTION_TSV,
    read_assembly_stats,
    read_excused_runs,
    read_selection,
    write_reads,
)
from repgenr.core.errors import WorkdirError
from repgenr.stages import assemble as stage
from repgenr.stages.assemble import AssembleParams, planned_layout, run
from repgenr.stages.assemble_qc import CHECKM2_CACHE

sys.path.insert(0, str(Path(__file__).resolve().parent))
from assemble_fakes import (  # noqa: E402  # noqa: E402
    FakeAssembler,
    FakeClassifier,
    fake_checkm2,
    read_row,
    register_fake_assembler,
    register_fake_classifier,
    unregister_fake_assembler,
    unregister_fake_classifier,
)

_LOG = logging.getLogger("test")


@pytest.fixture
def fakes():
    register_fake_assembler()
    register_fake_classifier()
    yield
    unregister_fake_assembler()
    unregister_fake_classifier()


def _db(path: Path, text: str = "db\n") -> str:
    path.write_text(text, encoding="utf-8")
    return str(path)


def _prepare(workdir: Path, rows) -> WorkdirContext:
    ctx = WorkdirContext(workdir, create=True)
    write_reads(workdir / READS_TSV, rows)
    return ctx


# --- stored CheckM2 scores -----------------------------------------------------------


def _checkm2(monkeypatch, quality, version="1.1.0") -> list:
    calls: list = []
    monkeypatch.setattr(stage, "preflight_checkm2", lambda: {"checkm2": version})
    monkeypatch.setattr(stage, "run_checkm2", fake_checkm2(quality, calls))
    return calls


def test_a_gate_change_reapplies_stored_checkm2_scores(
    workdir, tmp_path, fakes, monkeypatch, caplog
) -> None:
    calls = _checkm2(monkeypatch, {"SRR1.fasta": (98.0, 0.5), "SRR2.fasta": (70.0, 2.0)})
    db = _db(tmp_path / "checkm2.dmnd")
    ctx = _prepare(workdir, [read_row(tmp_path, "SRR1"), read_row(tmp_path, "SRR2")])
    run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=db))
    assert calls == [["SRR1.fasta", "SRR2.fasta"]]
    stored = json.loads((workdir / "assemblies" / "SRR1" / CHECKM2_CACHE).read_text())
    assert stored["completeness"] == 98.0 and stored["key"]["checkm2"] == "1.1.0"

    # Only the gate changes: the stored scores are applied, CheckM2 does not run.
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.INFO):
        run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=db, min_completeness=80.0))
    assert len(calls) == 1
    assert any("reusing the stored scores of 2 of 2" in r.getMessage() for r in caplog.records)
    rows = {r.accession: r for r in read_selection(workdir / SELECTION_TSV)}
    assert set(rows) == {"SRR1"} and rows["SRR1"].completeness == 98.0
    excused = {e.run_accession: e.reason for e in read_excused_runs(workdir / EXCUSED_RUNS_TSV)}
    assert excused["SRR2"].startswith("qc_failed: completeness 70.0 (min 80)")


def test_stored_checkm2_scores_follow_the_database_version_and_contigs(
    workdir, tmp_path, fakes, monkeypatch
) -> None:
    calls = _checkm2(monkeypatch, {"SRR1.fasta": (98.0, 0.5)})
    db = _db(tmp_path / "checkm2.dmnd")
    ctx = _prepare(workdir, [read_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=db))
    assert len(calls) == 1
    # Another database path scores again.
    other = _db(tmp_path / "other.dmnd")
    run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=other))
    assert len(calls) == 2
    # The same path replaced by a database of another size scores again.
    _db(tmp_path / "other.dmnd", "a newer release\n")
    run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=other))
    assert len(calls) == 3
    # Another CheckM2 version scores again.
    calls = _checkm2(monkeypatch, {"SRR1.fasta": (98.0, 0.5)}, version="1.2.0")
    run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=other))
    assert len(calls) == 1
    # Changed contigs score again.
    contigs = workdir / "assemblies" / "SRR1" / "contigs.fasta"
    contigs.write_text(contigs.read_text() + ">SRR1_contig2\n" + "ACGT" * 200 + "\n")
    run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=other))
    assert len(calls) == 2


def test_a_run_without_a_score_is_not_stored(workdir, tmp_path, fakes, monkeypatch) -> None:
    calls = _checkm2(monkeypatch, {})
    db = _db(tmp_path / "checkm2.dmnd")
    ctx = _prepare(workdir, [read_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=db))
    assert not (workdir / "assemblies" / "SRR1" / CHECKM2_CACHE).exists()
    run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=db, max_contamination=5.0))
    assert len(calls) == 2


# --- genus renames ---------------------------------------------------------------------


def test_a_renamed_genus_with_the_same_epithet_is_flagged_genus_renamed(
    workdir, tmp_path, fakes, caplog
) -> None:
    FakeClassifier.lineages = {
        "SRR1.fasta": (
            "d__Bacteria;p__Bacillota;c__Bacilli;o__Mycoplasmatales;f__Metamycoplasmataceae;"
            "g__Metamycoplasma;s__Metamycoplasma arginini"
        ),
        "SRR2.fasta": "d__Bacteria;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis",
    }
    myco = dict(
        organism="Mycoplasmopsis arginini",
        family="Metamycoplasmataceae",
        genus="Mycoplasmopsis",
        species="arginini",
    )
    ctx = _prepare(workdir, [read_row(tmp_path, "SRR1", **myco), read_row(tmp_path, "SRR2")])
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.INFO):
        run(
            ctx,
            AssembleParams(
                assembler="fakeasm",
                classifier="fakecls",
                gtdb_sketch=_db(tmp_path / "db.sig.zip"),
            ),
        )
    stats = {s.run_accession: s for s in read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)}
    assert stats["SRR1"].taxonomy_flag == "genus_renamed"
    assert stats["SRR1"].label_source == "metadata"
    assert stats["SRR1"].filename == "Metamycoplasmataceae_Mycoplasmopsis_arginini_SRR1.fasta"
    assert stats["SRR2"].taxonomy_flag == "classifier_disagrees"
    record = ctx.config.stages["assemble"]
    assert record.params["n_disagree"] == 1 and record.params["n_genus_renamed"] == 1
    renamed = [r for r in caplog.records if "genus_renamed" in r.getMessage()]
    assert renamed and renamed[0].levelno == logging.INFO
    assert "g__Metamycoplasma" in renamed[0].getMessage()


@pytest.mark.parametrize(
    ("submitted", "gtdb", "same"),
    [
        ("arginini", "arginini", True),
        ("coli", "coli-A", True),
        ("unknown", "unknown", False),
        ("", "arginini", False),
        ("subtilis", "tularensis", False),
    ],
)
def test_epithet_comparison(submitted, gtdb, same) -> None:
    assert stage._same_epithet(submitted, gtdb) is same


# --- classifier memory budget ---------------------------------------------------------


def test_the_memory_budget_reaches_the_classifier(workdir, tmp_path, fakes) -> None:
    FakeClassifier.lineages = {}
    ctx = _prepare(workdir, [read_row(tmp_path, "SRR1")])
    run(
        ctx,
        AssembleParams(
            assembler="fakeasm",
            classifier="fakecls",
            gtdb_sketch=_db(tmp_path / "db.sig.zip"),
            memory_gb=3,
        ),
    )
    assert FakeClassifier.last_params is not None
    assert FakeClassifier.last_params.memory_gb == 3


# --- single-file PAIRED runs -----------------------------------------------------------


def _one_file_paired(tmp_path: Path, run_accession: str, **over):
    row = read_row(tmp_path, run_accession, layout="SINGLE", **over)
    from dataclasses import replace

    return replace(row, layout="PAIRED")


def test_planned_layout(tmp_path) -> None:
    from dataclasses import replace

    paired = read_row(tmp_path, "X", layout="SINGLE")
    one_file = replace(paired, layout="PAIRED")
    assert planned_layout(one_file) == "SINGLE"
    assert planned_layout(replace(one_file, platform="OXFORD_NANOPORE")) == "PAIRED"
    two = replace(one_file, fastq_urls=("a_1.fastq.gz", "a_2.fastq.gz"))
    assert planned_layout(two) == "PAIRED"
    assert planned_layout(replace(one_file, fastq_urls=())) == "PAIRED"


def test_a_pair_only_assembler_excuses_a_one_file_paired_run_before_download(
    workdir, tmp_path, fakes, monkeypatch, caplog
) -> None:
    # The file does not exist: a fetch attempt would excuse the run as
    # download_failed, so the reason shows the run was excused when planned.
    monkeypatch.setattr(FakeAssembler, "layouts", frozenset({"PAIRED"}))
    one = _one_file_paired(
        tmp_path,
        "ONEFILE",
        fastq_urls=(str(tmp_path / "absent" / "ONEFILE.fastq.gz"),),
    )
    ctx = _prepare(workdir, [read_row(tmp_path, "SRR1"), one])
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.INFO):
        n = run(ctx, AssembleParams(assembler="fakeasm"))
    assert n == 1
    excused = {e.run_accession: e.reason for e in read_excused_runs(workdir / EXCUSED_RUNS_TSV)}
    assert excused == {"ONEFILE": "unsupported_layout"}
    assert "ONEFILE" not in FakeAssembler.calls
    assert any("planned as single-end" in r.getMessage() for r in caplog.records)


def test_auto_assembles_a_one_file_paired_run_as_single_end(workdir, tmp_path, fakes) -> None:
    ctx = _prepare(workdir, [_one_file_paired(tmp_path, "ONEFILE")])
    assert run(ctx, AssembleParams(assembler="fakeasm")) == 1
    assert FakeAssembler.layouts_seen == {"ONEFILE": "SINGLE"}


# --- every run excused -----------------------------------------------------------------


def test_a_total_failure_clears_the_previous_genome_set(
    workdir, tmp_path, fakes, monkeypatch
) -> None:
    from repgenr.stages.dereplicate import DereplicateParams
    from repgenr.stages.dereplicate import precheck as derep_precheck

    _checkm2(monkeypatch, {"SRR1.fasta": (60.0, 1.0)})
    og = tmp_path / "Fam_Gen_sp_GCF_000009.1.fasta"
    og.write_text(">og\nACGT\n", encoding="utf-8")
    db = _db(tmp_path / "checkm2.dmnd")
    ctx = _prepare(workdir, [read_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=db, outgroup=str(og)))
    assert len(list(ctx.genomes_dir.iterdir())) == 1

    with pytest.raises(WorkdirError, match="previous genome set was cleared") as info:
        run(
            ctx,
            AssembleParams(
                assembler="fakeasm", checkm2_db=db, outgroup=str(og), min_completeness=90.0
            ),
        )
    assert info.value.exit_code == 3
    assert ctx.genomes_dir.is_dir() and not any(ctx.genomes_dir.iterdir())
    assert read_selection(workdir / SELECTION_TSV) == []
    assert not (workdir / ASSEMBLY_STATS_TSV).exists()
    assert not ctx.outgroup_dir.exists()
    assert not (workdir / "outgroup_accession.txt").exists()
    assert ctx.manifest.all_genomes(include_outgroup=True) == []
    assert [e.run_accession for e in read_excused_runs(workdir / EXCUSED_RUNS_TSV)] == ["SRR1"]
    # The finished assembly stays for a later call with a looser gate.
    assert (workdir / "assemblies" / "SRR1" / "assembly.ok").exists()
    # dereplicate refuses rather than running on the old set.
    with pytest.raises(WorkdirError) as derep:
        derep_precheck(ctx, DereplicateParams())
    assert derep.value.exit_code == 3

    # A looser gate restores the set without assembling again.
    run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=db, outgroup=str(og)))
    assert FakeAssembler.calls == ["SRR1"]
    assert len(read_selection(workdir / SELECTION_TSV)) == 2


def test_a_total_failure_under_append_leaves_the_selection(
    workdir, tmp_path, fakes, monkeypatch
) -> None:
    from repgenr.core.contracts import SelectionRow, write_selection

    ctx = WorkdirContext(workdir, create=True)
    ctx.genomes_dir.mkdir(parents=True)
    kept = SelectionRow("GCF_1", "Fam", "Gen", "sp", False, "Fam_Gen_sp_GCF_1.fasta")
    (ctx.genomes_dir / kept.filename).write_text(">g\nACGT\n", encoding="utf-8")
    write_selection(workdir / SELECTION_TSV, [kept])
    FakeAssembler.fail_runs = frozenset({"SRR1"})
    write_reads(workdir / READS_TSV, [read_row(tmp_path, "SRR1")])
    with pytest.raises(WorkdirError, match="existing selection is unchanged"):
        run(ctx, AssembleParams(assembler="fakeasm", append=True))
    assert [r.accession for r in read_selection(workdir / SELECTION_TSV)] == ["GCF_1"]
    assert (ctx.genomes_dir / kept.filename).exists()


def test_a_total_failure_through_the_cli_marks_the_stage_and_dereplicate_exits_3(
    workdir, tmp_path, fakes, monkeypatch
) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    _checkm2(monkeypatch, {"SRR1.fasta": (60.0, 1.0)})
    db = _db(tmp_path / "checkm2.dmnd")
    WorkdirContext(workdir, create=True)
    write_reads(workdir / READS_TSV, [read_row(tmp_path, "SRR1")])
    runner = CliRunner()
    base = ["assemble", "-wd", str(workdir), "--assembler", "fakeasm", "--checkm2-db", db]
    first = runner.invoke(app, base)
    assert first.exit_code == 0, first.output
    failed = runner.invoke(app, [*base, "--min-completeness", "90"])
    assert failed.exit_code == 3, failed.output
    status = runner.invoke(app, ["status", "-wd", str(workdir)])
    line = next(ln for ln in status.output.splitlines() if "assemble" in ln)
    assert "[interrupted]" in line, status.output
    derep = runner.invoke(app, ["dereplicate", "-wd", str(workdir)])
    assert derep.exit_code == 3, derep.output


def test_a_first_call_in_which_every_run_fails_says_no_set_was_written(
    workdir, tmp_path, fakes
) -> None:
    FakeAssembler.fail_runs = frozenset({"SRR1"})
    ctx = _prepare(workdir, [read_row(tmp_path, "SRR1")])
    with pytest.raises(WorkdirError, match="no genome set was written"):
        run(ctx, AssembleParams(assembler="fakeasm"))
    assert read_selection(workdir / SELECTION_TSV) == []
    assert not any(ctx.genomes_dir.iterdir())
