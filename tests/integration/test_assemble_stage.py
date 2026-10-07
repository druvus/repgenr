"""The assemble stage: fetch reads, assemble, label, and write the genome contract."""

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
    ReadRow,
    SelectionRow,
    read_assembly_stats,
    read_excused_runs,
    read_selection,
    write_reads,
    write_selection,
)
from repgenr.core.errors import WorkdirError
from repgenr.core.integrity import check_genome_completeness
from repgenr.stages.assemble import AssembleParams, run

sys.path.insert(0, str(Path(__file__).resolve().parent))
from assemble_fakes import FakeAssembler as _FakeAssembler  # noqa: E402
from assemble_fakes import FakeClassifier as _FakeClassifier  # noqa: E402
from assemble_fakes import fake_checkm2 as _fake_checkm2  # noqa: E402
from assemble_fakes import read_row as _row  # noqa: E402
from assemble_fakes import (  # noqa: E402
    register_fake_assembler,
    register_fake_classifier,
    unregister_fake_assembler,
    unregister_fake_classifier,
)

_LOG = logging.getLogger("test")


@pytest.fixture
def fake_assembler():
    register_fake_assembler()
    yield
    unregister_fake_assembler()


def _db(path: Path) -> str:
    """A stand-in database file: the stage refuses database paths that do not exist."""
    path.write_text("db\n", encoding="utf-8")
    return str(path)


def _prepare(workdir: Path, rows: list[ReadRow]) -> WorkdirContext:
    ctx = WorkdirContext(workdir, create=True)
    write_reads(workdir / READS_TSV, rows)
    return ctx


def test_assemble_writes_the_genome_contract(workdir: Path, tmp_path: Path, fake_assembler) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")])
    n = run(ctx, AssembleParams(assembler="fakeasm", threads=4, jobs=2))
    assert n == 2
    names = sorted(p.name for p in ctx.genomes_dir.iterdir())
    assert names == [
        "Francisellaceae_Francisella_tularensis_SRR1.fasta",
        "Francisellaceae_Francisella_tularensis_SRR2.fasta",
    ]
    text = (ctx.genomes_dir / names[0]).read_text(encoding="utf-8")
    assert text.startswith(">SRR1_contig1\n") and "tiny" not in text  # filtered and renamed
    rows = read_selection(workdir / SELECTION_TSV)
    assert {r.accession for r in rows} == {"SRR1", "SRR2"}
    assert {r.filename for r in rows} == set(names)
    stats = {s.run_accession: s for s in read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)}
    assert stats["SRR1"].n_contigs == 1 and stats["SRR1"].total_length == 1200
    assert stats["SRR1"].assembler == "fakeasm" and stats["SRR1"].est_coverage == 1000.0
    assert stats["SRR1"].ncbi_taxonomy == "Francisellaceae;Francisella;tularensis"
    manifest = {g.accession: g for g in ctx.manifest.all_genomes()}
    assert manifest["SRR1"].source == "sra" and manifest["SRR1"].species == "tularensis"
    record = ctx.config.stages["assemble"]
    assert record.params["n_assembled"] == 2 and record.params["n_excused"] == 0
    assert record.tool_versions == {"fakeasm": "1.0"}
    assert check_genome_completeness(ctx.genomes_dir, workdir, logger=_LOG) == []
    assert not (ctx.scratch_dir / "assemble" / "SRR1").exists()  # reads and scratch removed


def test_a_second_run_skips_finished_assemblies(
    workdir: Path, tmp_path: Path, fake_assembler
) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm"))
    assert _FakeAssembler.calls == ["SRR1"]
    marker = workdir / "assemblies" / "SRR1" / "assembly.ok"
    assert json.loads(marker.read_text(encoding="utf-8"))["assembler"] == "fakeasm"
    run(ctx, AssembleParams(assembler="fakeasm"))
    assert _FakeAssembler.calls == ["SRR1"]  # not assembled again
    assert len(read_selection(workdir / SELECTION_TSV)) == 1


def test_failures_are_excused_and_the_rest_proceed(
    workdir: Path, tmp_path: Path, fake_assembler
) -> None:
    rows = [
        _row(tmp_path, "SRR1"),
        _row(tmp_path, "NOMIRROR", fastq_urls=(), fastq_md5=(), fastq_bytes=()),
        _row(tmp_path, "BADSUM", fastq_md5=("0" * 32, "0" * 32)),
        _row(tmp_path, "CRASH"),
        _row(tmp_path, "IONT", platform="ION_TORRENT"),
    ]
    _FakeAssembler.fail_runs = frozenset({"CRASH"})
    ctx = _prepare(workdir, rows)
    n = run(ctx, AssembleParams(assembler="auto"))
    assert n == 1
    excused = {
        e.run_accession: (e.step, e.reason) for e in read_excused_runs(workdir / EXCUSED_RUNS_TSV)
    }
    assert excused["NOMIRROR"] == ("fetch", "no_fastq_mirror")
    assert excused["BADSUM"][0] == "fetch" and "download_failed" in excused["BADSUM"][1]
    assert excused["CRASH"][0] == "assemble" and "assembly_failed" in excused["CRASH"][1]
    assert excused["IONT"] == ("assemble", "unsupported_platform")
    assert [r.accession for r in read_selection(workdir / SELECTION_TSV)] == ["SRR1"]
    assert check_genome_completeness(ctx.genomes_dir, workdir, logger=_LOG) == []
    assert ctx.config.stages["assemble"].params["n_excused"] == 4


def test_an_excused_tool_failure_stays_on_one_tsv_line(
    workdir: Path, tmp_path: Path, fake_assembler
) -> None:
    # A ToolExecutionError carries the tool's output tail on further lines;
    # excused_runs.tsv must still hold one line per run for awk/cut readers.
    _FakeAssembler.fail_runs = frozenset({"CRASH"})
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "CRASH")])
    run(ctx, AssembleParams(assembler="fakeasm"))
    lines = (workdir / EXCUSED_RUNS_TSV).read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2
    run_accession, step, reason = lines[1].split("\t")
    assert (run_accession, step) == ("CRASH", "assemble")
    assert reason.startswith("assembly_failed:") and "boom" in reason


def test_auto_excuses_runs_whose_assembler_is_not_installed(
    workdir: Path, tmp_path: Path, fake_assembler, monkeypatch, caplog
) -> None:
    # Only an Illumina assembler is installed: the Illumina run assembles, the
    # ONT run (which flye would take) is excused with an accurate reason.
    from repgenr.assemblers import base as assemblers_base

    monkeypatch.setattr(_FakeAssembler, "read_types", frozenset({"ILLUMINA"}))
    monkeypatch.setattr(assemblers_base, "tool_available", lambda caps: caps.name == "fakeasm")
    rows = [
        _row(tmp_path, "SRR1"),
        _row(tmp_path, "ONT1", "OXFORD_NANOPORE", "SINGLE"),
        _row(tmp_path, "IONT", platform="ION_TORRENT"),
    ]
    ctx = _prepare(workdir, rows)
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.WARNING):
        n = run(ctx, AssembleParams(assembler="auto"))
    assert n == 1
    excused = {
        e.run_accession: (e.step, e.reason) for e in read_excused_runs(workdir / EXCUSED_RUNS_TSV)
    }
    assert excused == {
        "ONT1": ("assemble", "assembler_not_installed"),
        "IONT": ("assemble", "unsupported_platform"),
    }
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert any("OXFORD_NANOPORE" in w and "flye" in w and "--container" in w for w in warnings)
    assert check_genome_completeness(ctx.genomes_dir, workdir, logger=_LOG) == []


def test_auto_with_no_installed_assembler_is_a_missing_tool(
    workdir: Path, tmp_path: Path, monkeypatch
) -> None:
    # Nothing can be assembled because no assembler is installed: exit 4,
    # naming the adapters per platform, not exit 3 with every run excused.
    from repgenr.assemblers import base as assemblers_base
    from repgenr.core.errors import MissingBinaryError

    monkeypatch.setattr(assemblers_base, "tool_available", lambda caps: False)
    rows = [_row(tmp_path, "SRR1"), _row(tmp_path, "ONT1", "OXFORD_NANOPORE", "SINGLE")]
    ctx = _prepare(workdir, rows)
    with pytest.raises(MissingBinaryError, match="skesa") as info:
        run(ctx, AssembleParams(assembler="auto"))
    assert info.value.exit_code == 4 and "flye" in str(info.value)


def test_outgroup_fasta_is_staged(workdir: Path, tmp_path: Path, fake_assembler) -> None:
    og = tmp_path / "Fam_Gen_sp_GCF_000009.1.fasta"
    og.write_text(">og\nACGT\n", encoding="utf-8")
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", outgroup=str(og)))
    assert [p.name for p in ctx.outgroup_dir.iterdir()] == [og.name]
    assert (workdir / "outgroup_accession.txt").read_text().strip() == "GCF_000009.1"
    rows = read_selection(workdir / SELECTION_TSV)
    assert {r.accession for r in rows if r.is_outgroup} == {"GCF_000009.1"}


def test_missing_reads_table_is_a_workdir_error(workdir: Path, fake_assembler) -> None:
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(WorkdirError, match="reads"):
        run(ctx, AssembleParams(assembler="fakeasm"))


def test_threads_are_split_across_jobs(workdir: Path, tmp_path: Path, fake_assembler) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")])
    run(ctx, AssembleParams(assembler="fakeasm", threads=8, jobs=2))
    marker = json.loads(
        (workdir / "assemblies" / "SRR1" / "assembly.ok").read_text(encoding="utf-8")
    )
    assert marker["tool_stats"]["threads"] == 4


def test_excused_runs_are_written_even_when_nothing_assembles(workdir, tmp_path, fake_assembler):
    _FakeAssembler.fail_runs = frozenset({"SRR1"})
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    with pytest.raises(WorkdirError, match="None of the 1 runs"):
        run(ctx, AssembleParams(assembler="fakeasm"))
    assert read_excused_runs(workdir / EXCUSED_RUNS_TSV)[0].step == "assemble"


def test_jobs_default_is_one_when_long_reads_are_pending(workdir, tmp_path, fake_assembler, caplog):
    """Memory bounds concurrent assemblies; a long-read run gets the machine alone."""
    rows = [
        _row(tmp_path, "SRR1"),
        _row(
            tmp_path,
            "ONT1",
            platform="OXFORD_NANOPORE",
            layout="SINGLE",
            instrument_model="GridION",
        ),
    ]
    ctx = _prepare(workdir, rows)
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.INFO):
        run(ctx, AssembleParams(assembler="fakeasm", threads=8))
    assert any("with 1 concurrent job" in r.message for r in caplog.records)
    assert ctx.config.stages["assemble"].params["jobs"] is None  # the request, not the resolution


def test_jobs_default_is_two_for_short_reads(workdir, tmp_path, fake_assembler, caplog):
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")])
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.INFO):
        run(ctx, AssembleParams(assembler="fakeasm", threads=8))
    assert any("with 2 concurrent jobs" in r.message for r in caplog.records)


# --- quality and classification -----------------------------------------------------


@pytest.fixture
def fake_classifier():
    register_fake_classifier()
    yield
    unregister_fake_classifier()


def test_checkm2_quality_gates_and_feeds_the_selection(
    workdir, tmp_path, fake_assembler, monkeypatch
) -> None:
    from repgenr.stages import assemble as stage

    monkeypatch.setattr(stage, "preflight_checkm2", lambda: {"checkm2": "1.1.0"})
    monkeypatch.setattr(
        stage,
        "run_checkm2",
        _fake_checkm2(
            {
                "SRR1.fasta": (98.5, 0.4),
                "SRR2.fasta": (40.0, 15.0),
            }
        ),
    )
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")])
    n = run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=_db(tmp_path / "db")))
    assert n == 1
    rows = read_selection(workdir / SELECTION_TSV)
    assert [(r.accession, r.completeness, r.contamination) for r in rows] == [("SRR1", 98.5, 0.4)]
    assert ctx.manifest.quality()["Francisellaceae_Francisella_tularensis_SRR1.fasta"] == (
        98.5,
        0.4,
    )
    excused = {e.run_accession: e for e in read_excused_runs(workdir / EXCUSED_RUNS_TSV)}
    assert excused["SRR2"].step == "qc" and "qc_failed" in excused["SRR2"].reason
    stats = {s.run_accession: s for s in read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)}
    assert stats["SRR1"].completeness == 98.5 and "SRR2" not in stats
    assert ctx.config.stages["assemble"].params["checkm2_db"] == str(tmp_path / "db")


def test_checkm2_gate_and_missing_results_are_warned_about(
    workdir, tmp_path, fake_assembler, monkeypatch, caplog
) -> None:
    """Both reach the console under --quiet: an excused assembly and one kept unscored."""
    from repgenr.stages import assemble as stage

    monkeypatch.setattr(stage, "preflight_checkm2", lambda: {"checkm2": "1.1.0"})
    monkeypatch.setattr(stage, "run_checkm2", _fake_checkm2({"SRR2.fasta": (40.0, 15.0)}))
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")])
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.WARNING):
        run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=_db(tmp_path / "db")))
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert any(
        "SRR1" in w and "no quality" in w and "kept without quality values" in w for w in warnings
    )
    assert any("SRR2" in w and "qc_failed" in w and "completeness 40.0" in w for w in warnings)


def test_classifier_agreement_names_the_genome_with_gtdb_tokens(
    workdir, tmp_path, fake_assembler, fake_classifier
) -> None:
    _FakeClassifier.lineages = {
        "SRR1.fasta": (
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Francisellales;"
            "f__Francisellaceae;g__Francisella;s__Francisella tularensis_A"
        ),
        "SRR2.fasta": (
            "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;"
            "g__Bacillus;s__Bacillus subtilis"
        ),
    }
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")])
    run(
        ctx,
        AssembleParams(
            assembler="fakeasm",
            classifier="fakecls",
            gtdb_sketch=_db(tmp_path / "gtdb-rs226-reps.k31-sc10k.sig.zip"),
            gtdb_lineages=_db(tmp_path / "lineages.csv"),
        ),
    )
    rows = {r.accession: r for r in read_selection(workdir / SELECTION_TSV)}
    # agreement at genus: the GTDB tokens name the file
    assert rows["SRR1"].filename == "Francisellaceae_Francisella_tularensis-A_SRR1.fasta"
    assert rows["SRR1"].species == "tularensis-A"
    assert (ctx.genomes_dir / rows["SRR1"].filename).exists()
    # disagreement: NCBI tokens stay and the genome is flagged
    assert rows["SRR2"].filename == "Francisellaceae_Francisella_tularensis_SRR2.fasta"
    stats = {s.run_accession: s for s in read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)}
    assert stats["SRR1"].label_source == "classifier" and stats["SRR1"].taxonomy_flag == ""
    assert stats["SRR2"].label_source == "metadata"
    assert stats["SRR2"].taxonomy_flag == "classifier_disagrees"
    assert stats["SRR2"].gtdb_taxonomy.endswith("s__Bacillus subtilis")
    record = ctx.config.stages["assemble"]
    assert record.params["n_disagree"] == 1
    assert record.tool_versions["gtdb_sketch"] == "gtdb-rs226-reps.k31-sc10k"


# --- appending into an existing working directory -----------------------------------


def _gtdb_workdir(workdir: Path) -> WorkdirContext:
    """A workdir as metadata + genome leave it: two GTDB genomes and an outgroup."""
    from repgenr.core.manifest import GenomeRecord

    ctx = WorkdirContext(workdir, create=True)
    ctx.genomes_dir.mkdir(parents=True)
    rows = [
        SelectionRow("GCF_1", "Fam", "Gen", "sp", False, "Fam_Gen_sp_GCF_1.fasta", 99.0, 0.1),
        SelectionRow("GCF_2", "Fam", "Gen", "sp", False, "Fam_Gen_sp_GCF_2.fasta", 98.0, 0.2),
        SelectionRow("GCF_9", "Fam", "Out", "og", True, "Fam_Out_og_GCF_9.fasta"),
    ]
    for r in rows[:2]:
        (ctx.genomes_dir / r.filename).write_text(">g\nACGT\n", encoding="utf-8")
    ctx.outgroup_dir.mkdir(parents=True)
    (ctx.outgroup_dir / rows[2].filename).write_text(">o\nACGT\n", encoding="utf-8")
    (workdir / "outgroup_accession.txt").write_text("GCF_9\n", encoding="utf-8")
    write_selection(workdir / SELECTION_TSV, rows)
    ctx.manifest.replace_genomes(
        [
            GenomeRecord(
                r.accession,
                r.filename,
                "gtdb",
                r.family,
                r.genus,
                r.species,
                r.is_outgroup,
                completeness=r.completeness,
                contamination=r.contamination,
            )
            for r in rows
        ]
    )
    return ctx


def test_append_adds_assemblies_to_a_gtdb_workdir(workdir, tmp_path, fake_assembler) -> None:
    ctx = _gtdb_workdir(workdir)
    write_reads(workdir / READS_TSV, [_row(tmp_path, "SRR1")])
    n = run(ctx, AssembleParams(assembler="fakeasm", append=True))
    assert n == 1
    rows = {r.accession: r for r in read_selection(workdir / SELECTION_TSV)}
    assert set(rows) == {"GCF_1", "GCF_2", "GCF_9", "SRR1"}
    assert rows["GCF_1"].completeness == 99.0  # existing rows untouched
    assert rows["GCF_9"].is_outgroup and rows["SRR1"].filename.endswith("_SRR1.fasta")
    names = {p.name for p in ctx.genomes_dir.iterdir()}
    assert names == {"Fam_Gen_sp_GCF_1.fasta", "Fam_Gen_sp_GCF_2.fasta", rows["SRR1"].filename}
    assert (workdir / "outgroup_accession.txt").read_text().strip() == "GCF_9"
    sources = {g.accession: g.source for g in ctx.manifest.all_genomes(include_outgroup=True)}
    assert sources == {"GCF_1": "gtdb", "GCF_2": "gtdb", "GCF_9": "gtdb", "SRR1": "sra"}
    assert check_genome_completeness(ctx.genomes_dir, workdir, logger=_LOG) == []
    assert ctx.config.stages["assemble"].params["append"] is True


def test_append_replaces_its_own_earlier_rows(workdir, tmp_path, fake_assembler) -> None:
    ctx = _gtdb_workdir(workdir)
    write_reads(workdir / READS_TSV, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", append=True))
    write_reads(workdir / READS_TSV, [_row(tmp_path, "SRR1", species="holarctica")])
    for marker in (workdir / "assemblies").rglob("assembly.ok"):
        marker.unlink()  # force a re-assembly under the new label
    run(ctx, AssembleParams(assembler="fakeasm", append=True))
    rows = read_selection(workdir / SELECTION_TSV)
    assert [r.accession for r in rows].count("SRR1") == 1
    assert next(r for r in rows if r.accession == "SRR1").species == "holarctica"
    assert not (ctx.genomes_dir / "Francisellaceae_Francisella_tularensis_SRR1.fasta").exists()


def test_append_needs_an_existing_selection(workdir, tmp_path, fake_assembler) -> None:
    from repgenr.core.errors import UserInputError

    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    with pytest.raises(UserInputError, match="--append"):
        run(ctx, AssembleParams(assembler="fakeasm", append=True))


def test_databases_from_the_environment_are_recorded(
    workdir, tmp_path, fake_assembler, fake_classifier, monkeypatch
) -> None:
    sketch = tmp_path / "gtdb-rs226-reps.k31-sc10k.sig.zip"
    lineages = tmp_path / "lineages.csv"
    _db(sketch), _db(lineages)
    monkeypatch.setenv("REPGENR_GTDB_SKETCH", str(sketch))
    monkeypatch.setenv("REPGENR_GTDB_LINEAGES", str(lineages))
    _FakeClassifier.lineages = {
        "SRR1.fasta": "d__Bacteria;f__Francisellaceae;g__Francisella;s__Francisella tularensis"
    }
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", classifier="fakecls"))
    params = ctx.config.stages["assemble"].params
    assert params["gtdb_sketch"] == str(sketch) and params["gtdb_lineages"] == str(lineages)


def test_finished_runs_do_not_need_the_assembler_again(workdir, tmp_path, fake_assembler) -> None:
    """Re-running for QC or classification must not demand the assembler binary."""
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm"))
    calls: list[str] = []
    original = _FakeAssembler.preflight

    def spy(self):
        calls.append("preflight")
        return original(self)

    _FakeAssembler.preflight = spy  # type: ignore[method-assign]
    try:
        run(ctx, AssembleParams(assembler="fakeasm", max_contamination=5.0))
    finally:
        _FakeAssembler.preflight = original  # type: ignore[method-assign]
    assert calls == []


def test_disagreement_warning_names_the_compared_genera(
    workdir, tmp_path, fake_assembler, fake_classifier, caplog
) -> None:
    _FakeClassifier.lineages = {
        "SRR1.fasta": "d__Bacteria;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis"
    }
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    ctx.logger.addHandler(caplog.handler)
    run(
        ctx,
        AssembleParams(
            assembler="fakeasm", classifier="fakecls", gtdb_sketch=_db(tmp_path / "db.sig.zip")
        ),
    )
    messages = [r.getMessage() for r in caplog.records if "classifier_disagrees" in r.getMessage()]
    assert messages, [r.getMessage() for r in caplog.records]
    assert "genus Francisella" in messages[0] and "g__Bacillus" in messages[0]


# --- polishing -----------------------------------------------------------------------


@pytest.fixture
def fake_polisher():
    from assemble_fakes import register_fake_polisher, unregister_fake_polisher

    register_fake_polisher()
    yield
    unregister_fake_polisher()


def _ont_row(tmp_path, run="ONT1"):
    return _row(
        tmp_path, run, platform="OXFORD_NANOPORE", layout="SINGLE", instrument_model="GridION"
    )


def test_long_read_assemblies_are_polished_and_the_marker_says_so(
    workdir, tmp_path, fake_assembler, fake_polisher
) -> None:
    from assemble_fakes import FakePolisher

    ctx = _prepare(workdir, [_ont_row(tmp_path), _row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", polisher="fakepol", polish_rounds=2, threads=4))
    assert FakePolisher.calls == ["ONT1"]  # the Illumina run is not polished
    genome = next(p for p in ctx.genomes_dir.iterdir() if "ONT1" in p.name)
    assert "AAAA" in genome.read_text(encoding="utf-8")  # the polished sequence was kept
    marker = json.loads((workdir / "assemblies" / "ONT1" / "assembly.ok").read_text())
    assert marker["polisher"] == "fakepol" and marker["polish_rounds"] == 2
    stats = {s.run_accession: s for s in read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)}
    assert stats["ONT1"].polisher == "fakepol" and stats["ONT1"].total_length == 1204
    assert stats["SRR1"].polisher == "" and stats["SRR1"].total_length == 1200
    record = ctx.config.stages["assemble"]
    assert record.params["polishers_used"] == ["fakepol"]
    assert record.tool_versions["fakepol"] == "0.1"


def test_polisher_none_skips_polishing(workdir, tmp_path, fake_assembler, fake_polisher) -> None:
    from assemble_fakes import FakePolisher

    ctx = _prepare(workdir, [_ont_row(tmp_path)])
    run(ctx, AssembleParams(assembler="fakeasm", polisher="none"))
    assert FakePolisher.calls == []
    stats = read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)[0]
    assert stats.polisher == "" and stats.total_length == 1200


def test_a_failed_polish_excuses_the_run(workdir, tmp_path, fake_assembler, fake_polisher) -> None:
    from assemble_fakes import FakePolisher

    FakePolisher.fail_runs = frozenset({"ONT1"})
    ctx = _prepare(workdir, [_ont_row(tmp_path), _row(tmp_path, "SRR1")])
    n = run(ctx, AssembleParams(assembler="fakeasm", polisher="fakepol"))
    assert n == 1
    excused = read_excused_runs(workdir / EXCUSED_RUNS_TSV)
    assert excused[0].run_accession == "ONT1" and "polish_failed" in excused[0].reason


def test_auto_polisher_warns_when_the_accepting_polisher_is_not_installed(
    workdir, tmp_path, fake_assembler, monkeypatch, caplog
) -> None:
    from repgenr.polishers import base as polishers_base

    monkeypatch.setattr(polishers_base, "tool_available", lambda caps: False)
    ctx = _prepare(
        workdir, [_ont_row(tmp_path), _ont_row(tmp_path, "ONT2"), _row(tmp_path, "SRR1")]
    )
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.WARNING):
        n = run(ctx, AssembleParams(assembler="fakeasm", polisher="auto"))
    assert n == 3  # the runs are assembled unpolished, not excused
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    polish = [w for w in warnings if "polish" in w]
    assert len(polish) == 1, warnings
    assert "2 OXFORD_NANOPORE" in polish[0] and "medaka" in polish[0]
    assert "--container" in polish[0] and "--polisher none" in polish[0]
    stats = {s.run_accession: s for s in read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)}
    assert stats["ONT1"].polisher == ""


# --- reuse of finished runs ------------------------------------------------------------


def test_a_lower_contig_floor_assembles_finished_runs_again(
    workdir, tmp_path, fake_assembler, caplog
) -> None:
    """The finished contigs were filtered at the old floor; a lower one needs the raw
    assembly again, so the marker is not reused."""
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", min_contig_length=500))
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.INFO):
        run(ctx, AssembleParams(assembler="fakeasm", min_contig_length=2))
    assert _FakeAssembler.calls == ["SRR1", "SRR1"]
    genome = next(ctx.genomes_dir.iterdir()).read_text(encoding="utf-8")
    assert genome.count(">") == 2  # the 4 bp contig now passes
    assert any("min_contig_length 500 -> 2" in r.getMessage() for r in caplog.records)


def test_a_higher_contig_floor_refilters_finished_runs(workdir, tmp_path, fake_assembler) -> None:
    """Raising the floor needs no new assembly: the kept contigs are filtered again."""
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")])
    run(ctx, AssembleParams(assembler="fakeasm", min_contig_length=2))
    with pytest.raises(WorkdirError, match="None of the 2 runs"):
        run(ctx, AssembleParams(assembler="fakeasm", min_contig_length=5000))
    assert sorted(_FakeAssembler.calls) == ["SRR1", "SRR2"]  # not assembled again
    excused = read_excused_runs(workdir / EXCUSED_RUNS_TSV)
    assert [e.reason for e in excused] == ["assembly_failed: no contig of 5000 bp or more"] * 2
    # The refused floor left the finished contigs in place, so 1000 refilters them too.
    run(ctx, AssembleParams(assembler="fakeasm", min_contig_length=1000))
    assert sorted(_FakeAssembler.calls) == ["SRR1", "SRR2"]
    marker = json.loads((workdir / "assemblies" / "SRR1" / "assembly.ok").read_text())
    assert marker["settings"]["min_contig_length"] == 1000


def test_a_different_assembler_or_tool_arg_assembles_again(
    workdir, tmp_path, fake_assembler, fake_polisher
) -> None:
    from assemble_fakes import FakePolisher

    ctx = _prepare(workdir, [_ont_row(tmp_path)])
    run(ctx, AssembleParams(assembler="fakeasm", polisher="fakepol"))
    run(ctx, AssembleParams(assembler="fakeasm", polisher="fakepol", polish_rounds=3))
    assert FakePolisher.calls == ["ONT1", "ONT1"]
    run(ctx, AssembleParams(assembler="fakeasm", polisher="none"))
    assert _FakeAssembler.calls == ["ONT1", "ONT1", "ONT1"]
    # A tool argument no adapter of the run reads (a classifier's) changes nothing.
    run(ctx, AssembleParams(assembler="fakeasm", polisher="none", extra={"ksize": "21"}))
    assert _FakeAssembler.calls == ["ONT1", "ONT1", "ONT1"]
    stats = read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)[0]
    assert stats.polisher == "" and stats.total_length == 1200


def test_a_resumed_run_keeps_the_tool_versions_and_names_the_assembler(
    workdir, tmp_path, fake_assembler, fake_polisher
) -> None:
    ctx = _prepare(workdir, [_ont_row(tmp_path), _row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", polisher="fakepol"))
    first = ctx.config.stages["assemble"]
    assert first.tool == "fakeasm" and first.params["assembler"] == "fakeasm"
    run(ctx, AssembleParams(assembler="fakeasm", polisher="fakepol", max_contamination=5.0))
    record = ctx.config.stages["assemble"]
    assert record.tool_versions == {"fakeasm": "1.0", "fakepol": "0.1"}
    assert record.tool == "fakeasm"


def test_auto_records_the_assemblers_used(workdir, tmp_path, fake_assembler, monkeypatch) -> None:
    from repgenr.assemblers import base as assemblers_base

    monkeypatch.setattr(assemblers_base, "tool_available", lambda caps: caps.name == "fakeasm")
    monkeypatch.setattr(assemblers_base, "_PREFERENCE", ("fakeasm",))
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="auto"))
    record = ctx.config.stages["assemble"]
    assert record.tool == "fakeasm" and record.params["assembler"] == "auto"


def test_an_unreadable_marker_is_assembled_again(workdir, tmp_path, fake_assembler) -> None:
    """A marker cut short by a kill is not a finished run (no JSON traceback)."""
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm"))
    (workdir / "assemblies" / "SRR1" / "assembly.ok").write_text('{"assem', encoding="utf-8")
    run(ctx, AssembleParams(assembler="fakeasm", max_contamination=5.0))
    assert _FakeAssembler.calls == ["SRR1", "SRR1"]


def test_a_run_without_long_enough_contigs_leaves_no_reads(
    workdir, tmp_path, fake_assembler
) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")])
    with pytest.raises(WorkdirError):
        run(ctx, AssembleParams(assembler="fakeasm", min_contig_length=5000))
    assert not (ctx.scratch_dir / "assemble" / "SRR1").exists()
    assert not (workdir / "assemblies" / "SRR1" / "contigs.fasta").exists()


def test_a_rerun_with_nothing_to_fetch_needs_no_free_disk(
    workdir, tmp_path, fake_assembler, monkeypatch
) -> None:
    """A QC or classification rerun over finished runs downloads nothing."""
    import shutil as _shutil
    from collections import namedtuple

    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm"))
    usage = namedtuple("usage", "total used free")
    monkeypatch.setattr(_shutil, "disk_usage", lambda path: usage(1, 1, 1))
    assert run(ctx, AssembleParams(assembler="fakeasm", max_contamination=5.0)) == 1


# --- quality inputs are checked before any assembly ------------------------------------


def test_a_missing_checkm2_database_is_refused_before_assembling(
    workdir, tmp_path, fake_assembler
) -> None:
    from repgenr.core.errors import UserInputError

    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    with pytest.raises(UserInputError, match="--checkm2-db"):
        run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=str(tmp_path / "nope.dmnd")))
    assert _FakeAssembler.calls == []


def test_a_missing_checkm2_binary_is_found_before_assembling(
    workdir, tmp_path, fake_assembler, monkeypatch
) -> None:
    from repgenr.core.errors import MissingBinaryError
    from repgenr.stages import assemble as stage

    def absent():
        raise MissingBinaryError("checkm2 not found")

    monkeypatch.setattr(stage, "preflight_checkm2", absent)
    db = tmp_path / "db.dmnd"
    db.write_text("x", encoding="utf-8")
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    with pytest.raises(MissingBinaryError):
        run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=str(db)))
    assert _FakeAssembler.calls == []


def test_sourmash_without_lineages_is_refused_before_assembling(
    workdir, tmp_path, fake_assembler, monkeypatch
) -> None:
    from repgenr.core.errors import UserInputError

    monkeypatch.delenv("REPGENR_GTDB_LINEAGES", raising=False)
    sketch = tmp_path / "gtdb.sig.zip"
    sketch.write_text("x", encoding="utf-8")
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    with pytest.raises(UserInputError, match="lineages"):
        run(ctx, AssembleParams(assembler="fakeasm", gtdb_sketch=str(sketch)))
    with pytest.raises(UserInputError, match="--gtdb-lineages"):
        run(
            ctx,
            AssembleParams(
                assembler="fakeasm", gtdb_sketch=str(sketch), gtdb_lineages=str(tmp_path / "no")
            ),
        )
    assert _FakeAssembler.calls == []


def test_verified_reads_of_an_interrupted_run_are_not_fetched_again(
    workdir, tmp_path, fake_assembler, monkeypatch
) -> None:
    """A kill during assembly leaves the checksummed FASTQ files in scratch; the
    rerun keeps them and fetches only a file that is missing or fails its checksum."""
    from repgenr.stages import assemble as stage

    row = _row(tmp_path, "SRR1")
    ctx = _prepare(workdir, [row])
    run_scratch = ctx.scratch_dir / "assemble" / "SRR1"
    run_scratch.mkdir(parents=True)
    kept = run_scratch / Path(row.fastq_urls[0]).name
    kept.write_bytes(Path(row.fastq_urls[0]).read_bytes())
    (run_scratch / Path(row.fastq_urls[1]).name).write_bytes(b"truncated")
    (run_scratch / "asm").mkdir()
    (run_scratch / "asm" / "partial.fa").write_text(">x\nA\n", encoding="utf-8")
    copied: list[str] = []
    original = stage.shutil.copy2
    monkeypatch.setattr(
        stage.shutil, "copy2", lambda src, dst: (copied.append(Path(src).name), original(src, dst))
    )
    assert run(ctx, AssembleParams(assembler="fakeasm")) == 1
    assert copied == [Path(row.fastq_urls[1]).name]


def test_scratch_clearing_tolerates_files_that_vanish(tmp_path, monkeypatch) -> None:
    """On exFAT, macOS removes the AppleDouble twin '._asm' along with 'asm'."""
    from repgenr.stages import assemble as stage

    row = _row(tmp_path, "SRR1")
    scratch = tmp_path / "scratch" / "SRR1"
    (scratch / "asm").mkdir(parents=True)
    (scratch / "._asm").write_bytes(b"x")
    original = stage.remove_tree

    def remove_with_twin(path):
        original(path)
        (path.parent / f"._{path.name}").unlink(missing_ok=True)

    monkeypatch.setattr(stage, "remove_tree", remove_with_twin)
    entries = sorted(scratch.iterdir(), key=lambda p: p.name != "asm")  # 'asm' first
    monkeypatch.setattr(type(scratch), "iterdir", lambda self: iter(entries))
    stage._clear_scratch(row, scratch)
    assert not (scratch / "asm").exists() and not (scratch / "._asm").exists()


def test_failed_downloads_are_named_with_how_to_retry(
    workdir, tmp_path, fake_assembler, caplog
) -> None:
    """A repeat with the same settings skips the stage, so the log says how to retry."""
    rows = [_row(tmp_path, "SRR1"), _row(tmp_path, "BADSUM", fastq_md5=("0" * 32, "0" * 32))]
    ctx = _prepare(workdir, rows)
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.WARNING):
        run(ctx, AssembleParams(assembler="fakeasm"))
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert any("1 run(s) could not be fetched" in w and "--force" in w for w in warnings)


def test_a_marker_without_settings_still_gets_a_higher_floor(
    workdir, tmp_path, fake_assembler
) -> None:
    """Markers written before the settings were recorded are reused, at the requested floor."""
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", min_contig_length=2))
    marker_path = workdir / "assemblies" / "SRR1" / "assembly.ok"
    marker = json.loads(marker_path.read_text())
    del marker["settings"]
    marker_path.write_text(json.dumps(marker), encoding="utf-8")
    run(ctx, AssembleParams(assembler="fakeasm", min_contig_length=500))
    assert _FakeAssembler.calls == ["SRR1"]
    assert next(ctx.genomes_dir.iterdir()).read_text(encoding="utf-8").count(">") == 1
    assert "settings" not in json.loads(marker_path.read_text())


def test_a_kill_during_a_refilter_never_leaves_a_marker_over_other_contigs(
    workdir, tmp_path, fake_assembler, monkeypatch
) -> None:
    """The marker goes first: a kill after the contigs are replaced leaves no marker,
    so the run is assembled again instead of reused with the wrong statistics."""
    from repgenr.stages import assemble as stage

    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", min_contig_length=2))
    marker = workdir / "assemblies" / "SRR1" / "assembly.ok"
    seen: list[bool] = []
    original = stage._write_marker

    def killed(path, data):
        seen.append(path.exists())
        raise KeyboardInterrupt  # the kill, after the contigs were replaced

    monkeypatch.setattr(stage, "_write_marker", killed)
    with pytest.raises(KeyboardInterrupt):
        run(ctx, AssembleParams(assembler="fakeasm", min_contig_length=500))
    assert seen == [False] and not marker.exists()
    monkeypatch.setattr(stage, "_write_marker", original)
    run(ctx, AssembleParams(assembler="fakeasm", min_contig_length=500))
    assert _FakeAssembler.calls == ["SRR1", "SRR1"]
