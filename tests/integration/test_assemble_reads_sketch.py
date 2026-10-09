"""The reads sketch of the assemble stage: one sourmash call per run on its FASTQ
files, beside the assembler, with a fake sourmash (tests/conftest.py) and a fake
assembler (assemble_fakes.py). No real tool runs."""

from __future__ import annotations

import json
import logging
import sys
import threading
from pathlib import Path

import pytest
from typer.testing import CliRunner

from repgenr.assemblers.base import registry as asm_registry
from repgenr.cli.main import app
from repgenr.core import http, sketches
from repgenr.core.config import Config
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import (
    ASSEMBLY_STATS_TSV,
    EXCUSED_RUNS_TSV,
    READS_TSV,
    ReadRow,
    read_assembly_stats,
    read_excused_runs,
    write_reads,
)
from repgenr.core.errors import MissingBinaryError
from repgenr.stages import assemble as stage
from repgenr.stages.assemble import AssembleParams, run
from repgenr.stages.assemble_steps import AssembleRunParams, assemble_run

sys.path.insert(0, str(Path(__file__).resolve().parent))
from assemble_fakes import FakeAssembler as _FakeAssembler  # noqa: E402
from assemble_fakes import fake_checkm2 as _fake_checkm2  # noqa: E402
from assemble_fakes import read_row as _row  # noqa: E402
from assemble_fakes import register_fake_assembler, unregister_fake_assembler  # noqa: E402

_runner = CliRunner()
_PARAMS = "k=21,k=31,k=51,scaled=1000,abund"


@pytest.fixture
def fake_assembler():
    register_fake_assembler()
    yield
    unregister_fake_assembler()


def _prepare(workdir: Path, rows: list[ReadRow], caplog=None) -> WorkdirContext:  # noqa: ANN001
    ctx = WorkdirContext(workdir, create=True)
    write_reads(workdir / READS_TSV, rows)
    if caplog is not None:
        ctx.logger.addHandler(caplog.handler)
    return ctx


def _marker(workdir: Path, run_acc: str) -> dict:
    return json.loads((workdir / "assemblies" / run_acc / "assembly.ok").read_text("utf-8"))


def _sketch(workdir: Path, run_acc: str) -> Path:
    return workdir / "assemblies" / run_acc / "reads.sig.zip"


def _hidden(directory: Path) -> list[str]:
    return sorted(p.name for p in directory.iterdir() if p.name.startswith("."))


def _stats(workdir: Path) -> dict[str, bool]:
    return {
        s.run_accession: s.reads_sketch for s in read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)
    }


class _SyncAssembler(_FakeAssembler):
    """The fake assembler, waiting at a barrier the reads sketch also waits at."""

    barrier: threading.Barrier | None = None

    def assemble(self, reads, out_dir, params, logger):  # noqa: ANN001
        if type(self).barrier is not None:
            type(self).barrier.wait(timeout=10)
        return super().assemble(reads, out_dir, params, logger)


def test_reads_are_sketched_beside_the_assembler(
    workdir, tmp_path, fake_assembler, fake_sourmash, register_tool, caplog
) -> None:
    register_tool(asm_registry, "syncasm", _SyncAssembler)
    # Both must reach the barrier for either to go on: the sketch and the
    # assembler run at the same time, or the barrier times out.
    _SyncAssembler.barrier = threading.Barrier(2)
    fake_sourmash.reads_hook = lambda name: _SyncAssembler.barrier.wait(timeout=10)
    try:
        ctx = _prepare(workdir, [_row(tmp_path, "SRR1")], caplog)
        with caplog.at_level(logging.INFO):
            run(ctx, AssembleParams(assembler="syncasm", threads=4, jobs=1))
    finally:
        _SyncAssembler.barrier = None

    # One call on both files of the pair, named by run accession.
    assert [name for name, _ in fake_sourmash.reads_calls] == ["SRR1"]
    inputs = fake_sourmash.reads_calls[0][1]
    assert [Path(f).name for f in inputs] == ["SRR1_1.fastq.gz", "SRR1_2.fastq.gz"]
    sketch = _sketch(workdir, "SRR1")
    assert sketch.read_text("utf-8").splitlines()[0] == "SRR1"
    marker = _marker(workdir, "SRR1")
    assert marker["reads_sketch"] == "reads.sig.zip"
    assert marker["reads_sketch_params"] == _PARAMS
    assert marker["reads_sketch_version"] == "4.9.4"
    assert _stats(workdir) == {"SRR1": True}
    record = ctx.config.stages["assemble"]
    assert record.params["n_reads_sketches"] == 1 and record.params["reads_sketch"] is None
    assert record.tool_versions["sourmash"] == "4.9.4"
    assert "SRR1: reads sketch written in" in caplog.text
    # Never part of the genome sketch contract.
    names = sorted(p.name for p in (workdir / "sketches").iterdir())
    assert all("reads" not in n for n in names) and "SRR1" not in fake_sourmash.calls
    assert _hidden(sketch.parent) == []
    # The reads are removed after the sketch ended.
    assert not (ctx.scratch_dir / "assemble" / "SRR1").exists()


def test_the_assembler_gives_one_thread_to_the_sketch(
    workdir, tmp_path, fake_assembler, fake_sourmash
) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")])
    run(ctx, AssembleParams(assembler="fakeasm", threads=8, jobs=2))
    assert _marker(workdir, "SRR1")["tool_stats"]["threads"] == 3  # 8 // 2 - 1
    assert _marker(workdir, "SRR2")["tool_stats"]["threads"] == 3


def test_with_one_thread_each_the_sketch_runs_alongside(
    workdir, tmp_path, fake_assembler, fake_sourmash
) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", threads=1, jobs=1))
    assert _marker(workdir, "SRR1")["tool_stats"]["threads"] == 1
    assert _sketch(workdir, "SRR1").is_file()


def test_no_reads_sketch_keeps_every_thread_for_the_assembler(
    workdir, tmp_path, fake_assembler, fake_sourmash
) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", threads=4, jobs=1, reads_sketch=False))
    assert fake_sourmash.reads_calls == []
    assert not _sketch(workdir, "SRR1").exists()
    marker = _marker(workdir, "SRR1")
    assert marker["tool_stats"]["threads"] == 4
    assert marker["reads_sketch"] is None
    assert marker["reads_sketch_params"] is None and marker["reads_sketch_version"] is None
    assert _stats(workdir) == {"SRR1": False}


def test_without_sourmash_the_default_logs_once_and_sketches_nothing(
    workdir, tmp_path, fake_assembler, caplog
) -> None:
    # tests/conftest.py reports sourmash as unavailable.
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")], caplog)
    with caplog.at_level(logging.INFO):
        run(ctx, AssembleParams(assembler="fakeasm", threads=4, jobs=2))
    assert caplog.text.count("reads sketches not written") == 1
    assert _marker(workdir, "SRR1")["reads_sketch"] is None
    assert _marker(workdir, "SRR1")["tool_stats"]["threads"] == 2
    assert _stats(workdir) == {"SRR1": False, "SRR2": False}


def test_explicit_reads_sketch_without_sourmash_exits_before_any_download(
    workdir, tmp_path, fake_assembler, monkeypatch
) -> None:
    def missing(caps):  # noqa: ANN001, ANN202
        raise MissingBinaryError("Required binary 'sourmash' not found on PATH.")

    def no_download(*args, **kwargs):  # noqa: ANN002, ANN003, ANN202
        raise AssertionError("downloaded before the preflight")

    monkeypatch.setattr(sketches, "preflight", missing)
    monkeypatch.setattr(http, "download", no_download)
    row = _row(
        tmp_path,
        "SRR1",
        fastq_urls=("https://example.org/SRR1_1.fastq.gz", "https://example.org/SRR1_2.fastq.gz"),
    )
    _prepare(workdir, [row]).close()
    result = _runner.invoke(
        app, ["assemble", "-wd", str(workdir), "--assembler", "fakeasm", "--reads-sketch"]
    )
    assert result.exit_code == MissingBinaryError.exit_code == 4, result.output
    assert _FakeAssembler.calls == []
    assert not (workdir / "assemblies" / "SRR1").exists()
    assert "assemble" not in Config.load(workdir).stages


def test_a_failed_sketch_is_a_warning_and_leaves_no_file(
    workdir, tmp_path, fake_assembler, fake_sourmash, caplog
) -> None:
    fake_sourmash.fail = {"SRR1"}  # writes a partial file, then fails
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")], caplog)
    with caplog.at_level(logging.INFO):
        n = run(ctx, AssembleParams(assembler="fakeasm", threads=4, jobs=1))
    assert n == 1
    assert "SRR1: reads sketch not written" in caplog.text
    assert any(
        r.levelno == logging.WARNING and "reads sketch not written" in r.getMessage()
        for r in caplog.records
    )
    run_dir = workdir / "assemblies" / "SRR1"
    assert not _sketch(workdir, "SRR1").exists()
    assert _hidden(run_dir) == []  # the temporary file went with the failure
    assert _marker(workdir, "SRR1")["reads_sketch"] is None
    assert _stats(workdir) == {"SRR1": False}
    assert not (workdir / EXCUSED_RUNS_TSV).exists()


def test_an_explicit_flag_does_not_turn_a_failed_sketch_into_a_failure(
    workdir, tmp_path, fake_assembler, fake_sourmash
) -> None:
    fake_sourmash.fail = {"SRR1"}
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    assert run(ctx, AssembleParams(assembler="fakeasm", reads_sketch=True)) == 1
    assert _marker(workdir, "SRR1")["reads_sketch"] is None


def test_an_excused_run_keeps_no_reads_sketch(
    workdir, tmp_path, fake_assembler, fake_sourmash
) -> None:
    _FakeAssembler.fail_runs = frozenset({"CRASH"})
    rows = [
        _row(tmp_path, "SRR1"),
        _row(tmp_path, "CRASH"),
        _row(tmp_path, "BADSUM", fastq_md5=("0" * 32, "0" * 32)),
    ]
    ctx = _prepare(workdir, rows)
    run(ctx, AssembleParams(assembler="fakeasm", threads=2, jobs=1))
    excused = {e.run_accession for e in read_excused_runs(workdir / EXCUSED_RUNS_TSV)}
    assert excused == {"CRASH", "BADSUM"}
    # The crashed run was sketched while it was assembled; the sketch is removed.
    assert "CRASH" in [name for name, _ in fake_sourmash.reads_calls]
    assert not _sketch(workdir, "CRASH").exists()
    assert not _sketch(workdir, "BADSUM").exists()  # never fetched, never sketched
    assert "BADSUM" not in [name for name, _ in fake_sourmash.reads_calls]
    assert _sketch(workdir, "SRR1").is_file()


def test_a_run_excused_by_quality_loses_its_reads_sketch(
    workdir, tmp_path, fake_assembler, fake_sourmash, monkeypatch
) -> None:
    db = tmp_path / "checkm2.dmnd"
    db.write_text("db\n", encoding="utf-8")
    monkeypatch.setattr(stage, "preflight_checkm2", lambda: {"checkm2": "1.0"})
    monkeypatch.setattr(
        stage, "run_checkm2", _fake_checkm2({"SRR1.fasta": (99.0, 0.5), "LOW.fasta": (20.0, 1.0)})
    )
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "LOW")])
    run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=str(db)))
    assert _sketch(workdir, "SRR1").is_file()
    assert not _sketch(workdir, "LOW").exists()
    # The marker stays (contigs and scores are kept for a later gate) but names no sketch.
    assert _marker(workdir, "LOW")["reads_sketch"] is None
    assert _stats(workdir) == {"SRR1": True}


def test_a_finished_run_is_reused_and_not_sketched_again(
    workdir, tmp_path, fake_assembler, fake_sourmash
) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm"))
    assert len(fake_sourmash.reads_calls) == 1
    run(ctx, AssembleParams(assembler="fakeasm"))
    assert _FakeAssembler.calls == ["SRR1"]
    assert len(fake_sourmash.reads_calls) == 1
    assert _stats(workdir) == {"SRR1": True}


def test_kept_reads_fill_a_missing_sketch_without_assembling_again(
    workdir, tmp_path, fake_assembler, fake_sourmash, caplog
) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")], caplog)
    run(ctx, AssembleParams(assembler="fakeasm", keep_reads=True, reads_sketch=False))
    marker = _marker(workdir, "SRR1")
    # A marker written before reads sketches existed has no such fields.
    for key in ("reads_sketch", "reads_sketch_params", "reads_sketch_version"):
        marker.pop(key)
    (workdir / "assemblies" / "SRR1" / "assembly.ok").write_text(json.dumps(marker), "utf-8")
    assert (ctx.scratch_dir / "assemble" / "SRR1" / "SRR1_1.fastq.gz").is_file()

    with caplog.at_level(logging.INFO):
        run(ctx, AssembleParams(assembler="fakeasm"))
    assert _FakeAssembler.calls == ["SRR1"]  # not assembled again
    assert [name for name, _ in fake_sourmash.reads_calls] == ["SRR1"]
    assert _sketch(workdir, "SRR1").is_file()
    marker = _marker(workdir, "SRR1")
    assert marker["reads_sketch"] == "reads.sig.zip"
    assert marker["reads_sketch_params"] == _PARAMS
    assert marker["stats"]["n_contigs"] == 1  # the rest of the marker is unchanged
    assert _stats(workdir) == {"SRR1": True}
    assert "Sketching the kept reads of 1 finished run(s)" in caplog.text


def test_without_kept_reads_a_missing_sketch_stays_empty(
    workdir, tmp_path, fake_assembler, fake_sourmash, caplog
) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")], caplog)
    run(ctx, AssembleParams(assembler="fakeasm", reads_sketch=False))
    with caplog.at_level(logging.INFO):
        run(ctx, AssembleParams(assembler="fakeasm"))
    assert sorted(_FakeAssembler.calls) == ["SRR1", "SRR2"]
    assert fake_sourmash.reads_calls == []
    assert _marker(workdir, "SRR1")["reads_sketch"] is None
    assert caplog.text.count("their FASTQ files were not kept") == 1
    assert "2 finished run(s) have no reads sketch" in caplog.text


def test_a_kept_file_that_does_not_match_its_checksum_is_not_sketched(
    workdir, tmp_path, fake_assembler, fake_sourmash
) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", keep_reads=True, reads_sketch=False))
    (ctx.scratch_dir / "assemble" / "SRR1" / "SRR1_2.fastq.gz").write_bytes(b"changed")
    run(ctx, AssembleParams(assembler="fakeasm"))
    assert fake_sourmash.reads_calls == []
    assert _FakeAssembler.calls == ["SRR1"]


def test_a_reassembled_run_drops_its_earlier_sketch(
    workdir, tmp_path, fake_assembler, fake_sourmash
) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm"))
    assert _sketch(workdir, "SRR1").is_file()
    # A lower contig floor assembles again; without a sketch this time.
    run(ctx, AssembleParams(assembler="fakeasm", min_contig_length=100, reads_sketch=False))
    assert _FakeAssembler.calls == ["SRR1", "SRR1"]
    assert not _sketch(workdir, "SRR1").exists()
    assert _marker(workdir, "SRR1")["reads_sketch"] is None


def test_assemble_run_sketches_into_its_output_directory(
    tmp_path, fake_assembler, fake_sourmash
) -> None:
    reads_tsv = tmp_path / READS_TSV
    write_reads(reads_tsv, [_row(tmp_path, "SRR1"), _row(tmp_path, "CRASH")])
    _FakeAssembler.fail_runs = frozenset({"CRASH"})
    log = logging.getLogger("test")
    out = tmp_path / "runs" / "SRR1"
    assert assemble_run(
        AssembleRunParams(reads_tsv=reads_tsv, run="SRR1", out_dir=out, assembler="fakeasm"),
        log,
    )
    assert (out / "reads.sig.zip").is_file()
    assert json.loads((out / "assembly.ok").read_text("utf-8"))["reads_sketch_version"] == "4.9.4"
    crash = tmp_path / "runs" / "CRASH"
    assert not assemble_run(
        AssembleRunParams(reads_tsv=reads_tsv, run="CRASH", out_dir=crash, assembler="fakeasm"),
        log,
    )
    assert not (crash / "reads.sig.zip").exists()
    off = tmp_path / "runs_off" / "SRR1"
    assemble_run(
        AssembleRunParams(
            reads_tsv=reads_tsv, run="SRR1", out_dir=off, assembler="fakeasm", reads_sketch=False
        ),
        log,
    )
    assert not (off / "reads.sig.zip").exists()


def test_assemble_run_explicit_flag_without_sourmash_is_refused(
    tmp_path, fake_assembler, monkeypatch
) -> None:
    def missing(caps):  # noqa: ANN001, ANN202
        raise MissingBinaryError("Required binary 'sourmash' not found on PATH.")

    monkeypatch.setattr(sketches, "preflight", missing)
    reads_tsv = tmp_path / READS_TSV
    write_reads(reads_tsv, [_row(tmp_path, "SRR1")])
    with pytest.raises(MissingBinaryError):
        assemble_run(
            AssembleRunParams(
                reads_tsv=reads_tsv,
                run="SRR1",
                out_dir=tmp_path / "SRR1",
                assembler="fakeasm",
                reads_sketch=True,
            ),
            logging.getLogger("test"),
        )
    assert _FakeAssembler.calls == []


def test_a_fresh_download_is_checked_while_it_is_written(
    workdir, tmp_path, fake_assembler, monkeypatch
) -> None:
    """The fetch hands the md5 to the download and does not hash the file again."""
    local = _row(tmp_path, "SRR1")
    seen: list[tuple[str, str | None]] = []

    def download(url, dest, *, logger=None, md5=None):  # noqa: ANN001, ANN202
        seen.append((Path(url).name, md5))
        dest.write_bytes(Path(tmp_path / Path(url).name).read_bytes())
        return dest

    def no_verify(path, md5):  # noqa: ANN001, ANN202
        raise AssertionError(f"{path.name} hashed again after the download")

    monkeypatch.setattr(http, "download", download)
    monkeypatch.setattr(http, "verify_md5", no_verify)
    row = _row(
        tmp_path,
        "SRR1",
        fastq_urls=tuple(f"https://example.org/{Path(u).name}" for u in local.fastq_urls),
    )
    ctx = _prepare(workdir, [row])
    assert run(ctx, AssembleParams(assembler="fakeasm", reads_sketch=False)) == 1
    assert seen == [
        ("SRR1_1.fastq.gz", local.fastq_md5[0]),
        ("SRR1_2.fastq.gz", local.fastq_md5[1]),
    ]


def test_a_stats_table_without_the_reads_sketch_column_still_reads(tmp_path) -> None:
    old = tmp_path / ASSEMBLY_STATS_TSV
    header = (
        "run_accession\tfilename\tassembler\tn_contigs\ttotal_length\tn50\tlargest_contig\t"
        "est_coverage\tcompleteness\tcontamination\tncbi_taxonomy\tgtdb_taxonomy\t"
        "label_source\ttaxonomy_flag\tpolisher\n"
    )
    row = "SRR1\tF_G_s_SRR1.fasta\tskesa\t3\t1200\t500\t600\t10.00\t\t\tF;G;s\t\tmetadata\t\t\n"
    old.write_text(header + row, encoding="utf-8")
    [parsed] = read_assembly_stats(old)
    assert parsed.run_accession == "SRR1" and parsed.reads_sketch is False
    # Older still: no polisher column either.
    old.write_text(
        header.replace("\tpolisher", "") + row.rsplit("\t", 1)[0] + "\n", encoding="utf-8"
    )
    [parsed] = read_assembly_stats(old)
    assert parsed.polisher == "" and parsed.reads_sketch is False


def test_the_reads_sketch_command_takes_every_file_of_the_run(tmp_path) -> None:
    files = [tmp_path / "a_1.fastq.gz", tmp_path / "a_2.fastq.gz"]
    out = tmp_path / "reads.sig.zip"
    argv = [str(a) for a in sketches.reads_sketch_command(files, "SRR1", out)]
    assert argv == [
        "sourmash",
        "sketch",
        "dna",
        "-p",
        _PARAMS,
        "--name",
        "SRR1",
        "-o",
        str(out),
        *map(str, files),
    ]
    assert sketches.READS_SKETCH_PARAMS == sketches.SKETCH_PARAMS + ",abund"


def test_sketch_reads_needs_files(tmp_path) -> None:
    from repgenr.core.errors import WorkdirError

    with pytest.raises(WorkdirError, match="no FASTQ files"):
        sketches.sketch_reads([], "SRR1", tmp_path / "reads.sig.zip", logging.getLogger("t"))
