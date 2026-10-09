"""assemble --screen-reads: the taxon, fraction and duplicate checks on each
run's reads sketch, before an assembler runs. sourmash is faked throughout:
the reads sketch by FakeSourmash (tests/conftest.py), gather and tax genome by
tables written here, the sketch comparison by a lookup."""

from __future__ import annotations

import json
import logging
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from repgenr.classifiers import sourmash as sm
from repgenr.cli.main import app
from repgenr.core import http
from repgenr.core.config import Config
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import (
    ASSEMBLY_STATS_TSV,
    EXCUSED_RUNS_TSV,
    READS_TSV,
    AssemblyStatsRow,
    ReadRow,
    read_assembly_stats,
    read_excused_runs,
    read_selection,
    write_assembly_stats,
    write_reads,
)
from repgenr.core.errors import WorkdirError
from repgenr.stages import assemble as stage
from repgenr.stages import assemble_screen as screen
from repgenr.stages.assemble import AssembleParams, run
from repgenr.stages.assemble_steps import (
    AssembleRunParams,
    ReadsGatherParams,
    assemble_run,
    reads_gather,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from assemble_fakes import FakeAssembler as _FakeAssembler  # noqa: E402
from assemble_fakes import read_row as _row  # noqa: E402
from assemble_fakes import register_fake_assembler, unregister_fake_assembler  # noqa: E402

_runner = CliRunner()
_LOG = logging.getLogger("test")

_GAMMA = "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria"
FT = f"{_GAMMA};o__Francisellales;f__Francisellaceae;g__Francisella;s__Francisella tularensis"
EC = f"{_GAMMA};o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli"
MA = (
    "d__Bacteria;p__Bacillota;c__Bacilli;o__Mycoplasmatales;f__Metamycoplasmataceae;"
    "g__Metamycoplasma;s__Metamycoplasma arginini"
)

_TAX_HEADER = (
    "query_name,status,rank,fraction,lineage,query_md5,query_filename,"
    "f_weighted_at_rank,bp_match_at_rank,query_ani_at_rank\n"
)


class FakeGather:
    """``sourmash gather`` and ``tax genome`` from canned tables, by run accession.

    ``hits[run]`` is (lineage, weighted fraction) or None for no match. The
    query of a gather is the fake reads sketch, whose first line is the run.
    ``events`` is shared with the other fakes to check the order of the steps.
    """

    def __init__(self, events: list[tuple[str, str]]) -> None:
        self.hits: dict[str, tuple[str, float] | None] = {}
        self.events = events
        self.gathers: list[str] = []
        self.fail: set[str] = set()
        self._lock = threading.Lock()

    def run_tool(self, caps, command, *, logger, **kwargs) -> int:  # noqa: ANN001
        from repgenr.core.errors import ToolExecutionError

        argv = [str(c) for c in command]
        if argv[:2] == ["sourmash", "gather"]:
            run = Path(argv[2]).read_text(encoding="utf-8").splitlines()[0]
            with self._lock:
                self.gathers.append(run)
                self.events.append(("gather", run))
            assert argv[argv.index("-k") + 1] == "31"
            if run in self.fail:
                raise ToolExecutionError(argv, 1, output="gather died", tool="sourmash")
            out = Path(argv[argv.index("-o") + 1])
            hit = self.hits.get(run)
            lines = ["intersect_bp,f_unique_weighted,name,query_name"]
            if hit is not None:
                lines.append(f"500000,{hit[1]},GCF_000000001.1 {hit[0].split('s__')[-1]},{run}")
            out.write_text("\n".join(lines) + "\n", encoding="utf-8")
            return 0
        assert argv[:3] == ["sourmash", "tax", "genome"]
        assert argv[argv.index("--containment-threshold") + 1] == "0"
        gather_csv = Path(argv[argv.index("--gather-csv") + 1])
        run = gather_csv.read_text(encoding="utf-8").splitlines()[1].split(",")[-1]
        lineage, fraction = self.hits[run]  # type: ignore[misc]
        base = Path(argv[argv.index("--output-base") + 1])
        Path(str(base) + ".classifications.csv").write_text(
            _TAX_HEADER + f'{run},match,species,0.04,"{lineage}",x,y,{fraction},500000,0.9\n',
            encoding="utf-8",
        )
        return 0


@pytest.fixture
def events() -> list[tuple[str, str]]:
    return []


@pytest.fixture
def fake_gather(monkeypatch, events) -> FakeGather:
    fake = FakeGather(events)
    monkeypatch.setattr(screen, "run_tool", fake.run_tool)
    monkeypatch.setattr(sm, "run_tool", fake.run_tool)
    return fake


class AniTable(dict):
    """ANI estimates by (query run, reference run); unlisted pairs are 0.9.
    ``compared`` lists the pairs the duplicate check compared."""

    def __init__(self) -> None:
        super().__init__()
        self.compared: list[tuple[str, str]] = []


@pytest.fixture
def ani(monkeypatch) -> AniTable:
    """The sketch comparison by lookup; the loaded "hashes" of a fake sketch
    are its run accession."""
    table = AniTable()

    def load(path: Path, ksize: int = 31):  # noqa: ANN202
        return path.read_text(encoding="utf-8").splitlines()[0]

    def containment(query, reference, ksize=31):  # noqa: ANN001, ANN202
        table.compared.append((query, reference))
        return table.get((query, reference), 0.9)

    monkeypatch.setattr(screen, "load_reads_hashes", load)
    monkeypatch.setattr(screen, "containment_ani", containment)
    return table


@pytest.fixture
def assembler(events):
    register_fake_assembler()

    original = _FakeAssembler.assemble

    def assemble(self, reads, out_dir, params, logger):  # noqa: ANN001, ANN202
        events.append(("assemble", reads.run_accession))
        return original(self, reads, out_dir, params, logger)

    _FakeAssembler.assemble = assemble  # type: ignore[method-assign]
    yield
    _FakeAssembler.assemble = original  # type: ignore[method-assign]
    unregister_fake_assembler()


@pytest.fixture
def gtdb(tmp_path) -> tuple[str, str]:
    db = tmp_path / "db" / "gtdb-rs226-reps.k31.sig.zip"
    db.parent.mkdir()
    db.write_bytes(b"zip")
    lineages = tmp_path / "db" / "lineages.csv"
    lineages.write_text("ident,superkingdom\n", encoding="utf-8")
    return str(db), str(lineages)


@pytest.fixture
def setup(workdir, tmp_path, fake_sourmash, fake_gather, ani, assembler, gtdb, events, monkeypatch):
    """A workdir factory and the fakes, with fetch and sketch recorded in ``events``."""
    fetch = stage._fetch

    def recording_fetch(row, run_scratch, logger):  # noqa: ANN001, ANN202
        events.append(("fetch", row.run_accession))
        return fetch(row, run_scratch, logger)

    monkeypatch.setattr(stage, "_fetch", recording_fetch)
    fake_sourmash.reads_hook = lambda name: events.append(("sketch", name))

    class Setup:
        def __init__(self) -> None:
            self.gather = fake_gather
            self.ani = ani
            self.events = events
            self.sourmash = fake_sourmash

        def prepare(self, rows: list[ReadRow]) -> WorkdirContext:
            ctx = WorkdirContext(workdir, create=True)
            write_reads(workdir / READS_TSV, rows)
            return ctx

        def params(self, **over) -> AssembleParams:  # noqa: ANN003
            base = dict(
                assembler="fakeasm",
                threads=2,
                jobs=1,
                screen_reads=True,
                gtdb_sketch=gtdb[0],
                gtdb_lineages=gtdb[1],
                classifier="none",
            )
            base.update(over)
            return AssembleParams(**base)

        def row(self, run: str, **over) -> ReadRow:  # noqa: ANN003
            return _row(tmp_path, run, **over)

    return Setup()


def _excused(workdir: Path) -> dict[str, tuple[str, str]]:
    path = workdir / EXCUSED_RUNS_TSV
    if not path.exists():
        return {}
    return {e.run_accession: (e.step, e.reason) for e in read_excused_runs(path)}


def _screen_rows(workdir: Path) -> dict[str, dict[str, str]]:
    lines = (workdir / screen.SCREEN_READS_TSV).read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    assert header == [
        "run_accession",
        "top_match",
        "top_genus",
        "fraction",
        "duplicate_of",
        "decision",
        "reason",
    ]
    return {c[0]: dict(zip(header, c, strict=True)) for c in (ln.split("\t") for ln in lines[1:])}


def _selected(workdir: Path) -> list[str]:
    return [r.accession for r in read_selection(workdir / "selection.tsv")]


# --- the three checks -------------------------------------------------------------------


def test_a_matching_run_passes_in_the_order_fetch_sketch_screen_assemble(
    workdir, setup, caplog
) -> None:
    setup.gather.hits["SRR1"] = (FT, 0.93)
    ctx = setup.prepare([setup.row("SRR1")])
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.INFO):
        assert run(ctx, setup.params()) == 1
    assert setup.events == [
        ("fetch", "SRR1"),
        ("sketch", "SRR1"),
        ("gather", "SRR1"),
        ("assemble", "SRR1"),
    ]
    rows = _screen_rows(workdir)
    assert rows["SRR1"]["decision"] == "pass" and rows["SRR1"]["top_genus"] == "Francisella"
    assert rows["SRR1"]["fraction"] == "0.9300"
    assert rows["SRR1"]["top_match"] == "GCF_000000001.1 Francisella tularensis"
    stats = read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)
    assert [s.screen for s in stats] == ["pass"]
    header = (workdir / ASSEMBLY_STATS_TSV).read_text(encoding="utf-8").splitlines()[0]
    assert header.split("\t")[-1] == "screen"
    marker = json.loads((workdir / "assemblies" / "SRR1" / "assembly.ok").read_text("utf-8"))
    assert marker["screen"]["decision"] == "pass"
    assert marker["screen"]["settings"]["min_fraction"] == 0.5
    assert marker["screen"]["settings"]["dup_ani"] == 0.999
    assert (workdir / "assemblies" / "SRR1" / "screen.json").is_file()
    assert "SRR1: screen passed" in caplog.text
    params = ctx.config.stages["assemble"].params
    assert params["screen_reads"] is True and params["n_screened_out"] == 0


def test_a_renamed_genus_of_the_same_family_and_epithet_passes(workdir, setup) -> None:
    row = setup.row(
        "SRR1",
        organism="Mycoplasmopsis arginini",
        family="Metamycoplasmataceae",
        genus="Mycoplasmopsis",
        species="arginini",
    )
    setup.gather.hits["SRR1"] = (MA, 0.58)
    ctx = setup.prepare([row])
    assert run(ctx, setup.params()) == 1
    assert _screen_rows(workdir)["SRR1"]["reason"] == "genus_renamed"
    assert _screen_rows(workdir)["SRR1"]["top_genus"] == "Metamycoplasma"


def test_a_foreign_genus_is_excused_as_taxon_mismatch_before_any_assembler(workdir, setup) -> None:
    setup.gather.hits.update({"SRR1": (FT, 0.9), "SRR2": (EC, 0.95)})
    ctx = setup.prepare([setup.row("SRR1"), setup.row("SRR2")])
    assert run(ctx, setup.params()) == 1
    assert _FakeAssembler.calls == ["SRR1"]
    step, reason = _excused(workdir)["SRR2"]
    assert step == "screen" and reason.startswith("taxon_mismatch: submitted genus Francisella")
    assert "GTDB genus Escherichia" in reason
    # The sketch and the gather output stay for inspection.
    run_dir = workdir / "assemblies" / "SRR2"
    assert (run_dir / "reads.sig.zip").is_file()
    assert (run_dir / "screen" / "gather.csv").is_file()
    assert not (run_dir / "assembly.ok").exists()
    assert json.loads((run_dir / "screen.json").read_text("utf-8"))["decision"] == "excused"
    assert not (ctx.scratch_dir / "assemble" / "SRR2").exists()  # reads removed
    rows = _screen_rows(workdir)
    assert rows["SRR2"]["decision"] == "excused" and rows["SRR2"]["top_genus"] == "Escherichia"
    assert ctx.config.stages["assemble"].params["n_screened_out"] == 1
    assert [s.run_accession for s in read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)] == ["SRR1"]


def test_no_gtdb_match_is_a_taxon_mismatch(workdir, setup) -> None:
    setup.gather.hits.update({"SRR1": (FT, 0.9), "SRR2": None})
    ctx = setup.prepare([setup.row("SRR1"), setup.row("SRR2")])
    run(ctx, setup.params())
    assert _excused(workdir)["SRR2"][1].startswith("taxon_mismatch: no GTDB match")


def test_a_low_top_species_fraction_is_excused_as_host_dominated(workdir, setup) -> None:
    setup.gather.hits.update({"SRR1": (FT, 0.9), "SRR2": (FT, 0.21)})
    ctx = setup.prepare([setup.row("SRR1"), setup.row("SRR2")])
    run(ctx, setup.params())
    step, reason = _excused(workdir)["SRR2"]
    assert step == "screen"
    assert reason == "host_dominated: weighted fraction of the top species 0.210 (min 0.5)"
    # The threshold is a flag.
    assert run(ctx, setup.params(screen_min_fraction=0.2)) == 2


def test_a_duplicate_of_an_earlier_run_of_the_same_biosample_is_excused(workdir, setup) -> None:
    for r in ("SRR1", "SRR2", "SRR3"):
        setup.gather.hits[r] = (FT, 0.9)
    setup.ani[("SRR2", "SRR1")] = 0.9996
    rows = [
        setup.row("SRR1", biosample="SAMX", taxid="1"),
        setup.row("SRR2", biosample="SAMX", taxid="2"),  # another taxid, the same sample
        setup.row("SRR3", biosample="SAMY", taxid="3"),  # neither: never compared
    ]
    ctx = setup.prepare(rows)
    assert run(ctx, setup.params()) == 2
    step, reason = _excused(workdir)["SRR2"]
    assert step == "screen"
    assert reason == "duplicate_isolate: contained in SRR1 at ANI 0.9996 (min 0.999)"
    assert _screen_rows(workdir)["SRR2"]["duplicate_of"] == "SRR1"
    assert ("SRR2", "SRR1") in setup.ani.compared
    assert not any("SRR3" in pair for pair in setup.ani.compared)
    assert _FakeAssembler.calls == ["SRR1", "SRR3"]


def test_duplicates_are_compared_within_a_taxid_and_kept_below_the_threshold(
    workdir, setup
) -> None:
    for r in ("SRR1", "SRR2", "SRR3"):
        setup.gather.hits[r] = (FT, 0.9)
    setup.ani[("SRR2", "SRR1")] = 0.9995  # same taxid
    setup.ani[("SRR3", "SRR1")] = 0.9999  # another taxid and sample: not compared
    rows = [
        setup.row("SRR1", taxid="263"),
        setup.row("SRR2", taxid="263"),
        setup.row("SRR3", taxid="119857"),
    ]
    ctx = setup.prepare(rows)
    assert run(ctx, setup.params()) == 2
    assert _excused(workdir) == {
        "SRR2": ("screen", "duplicate_isolate: contained in SRR1 at ANI 0.9995 (min 0.999)")
    }
    # A lower ANI threshold is a flag; at 0.99995 SRR2 is kept.
    stage_dir = workdir / "assemblies"
    assert (stage_dir / "SRR2" / "reads.sig.zip").is_file()
    assert run(ctx, setup.params(screen_dup_ani=0.99995)) == 3


def test_the_duplicate_check_follows_reads_order_whatever_finishes_first(workdir, setup) -> None:
    """SRR1 (larger, first in reads.tsv) is sketched slowly; SRR2 waits for its
    decision instead of being accepted first and excusing SRR1."""
    for r in ("SRR1", "SRR2"):
        setup.gather.hits[r] = (FT, 0.9)
    setup.ani[("SRR2", "SRR1")] = 0.9999
    setup.ani[("SRR1", "SRR2")] = 0.9999

    def slow(name: str) -> None:
        setup.events.append(("sketch", name))
        if name == "SRR1":
            time.sleep(0.5)

    setup.sourmash.reads_hook = slow
    ctx = setup.prepare([setup.row("SRR1"), setup.row("SRR2")])
    assert run(ctx, setup.params(jobs=2, threads=4)) == 1
    assert _selected(workdir) == ["SRR1"]
    assert _excused(workdir)["SRR2"][1].startswith("duplicate_isolate: contained in SRR1")


def test_a_failed_gather_is_screen_failed_and_not_judged(workdir, setup) -> None:
    setup.gather.hits.update({"SRR1": (FT, 0.9), "SRR2": (FT, 0.9)})
    ctx = setup.prepare([setup.row("SRR1"), setup.row("SRR2")])
    assert run(ctx, setup.params(screen_reads=False)) == 2
    # Assembled before; now every run's screen fails: the genome set is kept.
    setup.gather.fail = {"SRR1", "SRR2"}
    with pytest.raises(WorkdirError, match="from an earlier assemble call was kept"):
        run(ctx, setup.params())
    assert {r[1].split(":")[0] for r in _excused(workdir).values()} == {"screen_failed"}
    assert _selected(workdir) == ["SRR1", "SRR2"]


def test_screened_out_runs_are_judged_and_clear_an_earlier_set(workdir, setup) -> None:
    setup.gather.hits.update({"SRR1": (EC, 0.9), "SRR2": (EC, 0.9)})
    ctx = setup.prepare([setup.row("SRR1"), setup.row("SRR2")])
    assert run(ctx, setup.params(screen_reads=False)) == 2
    setup.events.clear()
    with pytest.raises(WorkdirError, match="previous genome set was cleared"):
        run(ctx, setup.params())
    # Screened from the kept reads sketches: nothing fetched again.
    assert [e for e in setup.events if e[0] == "fetch"] == []
    assert {r[0] for r in _excused(workdir).values()} == {"screen"}
    # The contigs and markers stay for a later call.
    assert (workdir / "assemblies" / "SRR1" / "assembly.ok").is_file()


# --- resume ---------------------------------------------------------------------------------


def test_a_repeat_with_the_same_settings_screens_nothing_again(workdir, setup) -> None:
    setup.gather.hits.update({"SRR1": (FT, 0.9), "SRR2": (EC, 0.9), "SRR3": (FT, 0.1)})
    ctx = setup.prepare([setup.row(r) for r in ("SRR1", "SRR2", "SRR3")])
    assert run(ctx, setup.params()) == 1
    first = _excused(workdir)
    setup.events.clear()
    setup.gather.gathers.clear()
    assert run(ctx, setup.params()) == 1
    assert setup.events == []  # no fetch, sketch, gather or assembly
    assert _excused(workdir) == first
    assert set(_screen_rows(workdir)) == {"SRR1", "SRR2", "SRR3"}


def test_a_changed_threshold_rescreens_from_the_kept_sketch(workdir, setup) -> None:
    setup.gather.hits.update({"SRR1": (FT, 0.9), "SRR2": (FT, 0.3)})
    ctx = setup.prepare([setup.row("SRR1"), setup.row("SRR2")])
    assert run(ctx, setup.params()) == 1
    # A stricter fraction: both are screened again from their sketches, no fetch;
    # the finished SRR1 now fails as well and keeps its contigs.
    setup.events.clear()
    with pytest.raises(WorkdirError):
        run(ctx, setup.params(screen_min_fraction=0.95))
    assert sorted(setup.events) == [("gather", "SRR1"), ("gather", "SRR2")]
    assert _excused(workdir)["SRR1"][1].startswith("host_dominated")
    assert (workdir / "assemblies" / "SRR1" / "contigs.fasta").is_file()
    # The same stricter fraction again: both excuses are reused, nothing gathered.
    setup.events.clear()
    with pytest.raises(WorkdirError):
        run(ctx, setup.params(screen_min_fraction=0.95))
    assert setup.events == []
    # A laxer one: SRR1 passes from its sketch, SRR2 passes from its sketch and
    # only then is fetched and assembled, with no second gather or sketch.
    setup.events.clear()
    assert run(ctx, setup.params(screen_min_fraction=0.2)) == 2
    assert sorted(setup.events) == [
        ("assemble", "SRR2"),
        ("fetch", "SRR2"),
        ("gather", "SRR1"),
        ("gather", "SRR2"),
    ]
    stats = read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)
    assert [s.screen for s in stats] == ["pass", "pass"]


def test_a_finished_run_assembled_unscreened_is_screened_from_its_sketch(workdir, setup) -> None:
    setup.gather.hits.update({"SRR1": (FT, 0.9)})
    ctx = setup.prepare([setup.row("SRR1")])
    run(ctx, setup.params(screen_reads=False))
    assert read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)[0].screen == ""
    setup.events.clear()
    run(ctx, setup.params())
    assert setup.events == [("gather", "SRR1")]
    assert read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)[0].screen == "pass"
    marker = json.loads((workdir / "assemblies" / "SRR1" / "assembly.ok").read_text("utf-8"))
    assert marker["screen"]["decision"] == "pass"


def test_a_finished_run_without_a_reads_sketch_is_not_screened(workdir, setup, caplog) -> None:
    ctx = setup.prepare([setup.row("SRR1")])
    run(ctx, setup.params(screen_reads=False, reads_sketch=False))
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.INFO):
        assert run(ctx, setup.params()) == 1
    assert "1 finished run(s) have no reads sketch and are not screened" in caplog.text
    assert setup.gather.gathers == []


# --- --max-runs ------------------------------------------------------------------------------


def test_a_screened_out_run_frees_its_place_in_the_next_call(workdir, setup) -> None:
    setup.gather.hits.update({"SRR1": (EC, 0.9), "SRR2": (FT, 0.9), "SRR3": (FT, 0.9)})
    ctx = setup.prepare([setup.row(r) for r in ("SRR1", "SRR2", "SRR3")])
    assert run(ctx, setup.params(max_runs=2)) == 1
    assert _excused(workdir) == {
        "SRR1": ("screen", _excused(workdir)["SRR1"][1]),
        "SRR3": ("assemble", "deferred"),
    }
    # The same N again: SRR1 is excused from its decision before the places
    # are counted, so SRR3 is assembled.
    setup.events.clear()
    assert run(ctx, setup.params(max_runs=2)) == 2
    assert ("fetch", "SRR1") not in setup.events and ("fetch", "SRR3") in setup.events
    assert _selected(workdir) == ["SRR2", "SRR3"]


# --- refusals --------------------------------------------------------------------------------


def test_screen_reads_without_a_gtdb_sketch_exits_2_before_any_download(
    workdir, tmp_path, fake_sourmash, monkeypatch
) -> None:
    def no_download(*args, **kwargs):  # noqa: ANN002, ANN003, ANN202
        raise AssertionError("downloaded before the check")

    monkeypatch.setattr(http, "download", no_download)
    monkeypatch.delenv(stage.GTDB_SKETCH_ENV, raising=False)
    monkeypatch.delenv(stage.GTDB_LINEAGES_ENV, raising=False)
    row = _row(tmp_path, "SRR1", fastq_urls=("https://example.org/SRR1_1.fastq.gz",))
    ctx = WorkdirContext(workdir, create=True)
    write_reads(workdir / READS_TSV, [row])
    ctx.close()
    result = _runner.invoke(app, ["assemble", "-wd", str(workdir), "--screen-reads"])
    assert result.exit_code == 2, result.output
    assert "--screen-reads needs the GTDB sketch" in result.output
    assert not (workdir / "assemblies" / "SRR1").exists()
    assert "assemble" not in Config.load(workdir).stages


def test_screen_reads_with_no_reads_sketch_is_refused(workdir, setup) -> None:
    from repgenr.core.errors import UserInputError

    ctx = setup.prepare([setup.row("SRR1")])
    with pytest.raises(UserInputError, match="--no-reads-sketch"):
        run(ctx, setup.params(reads_sketch=False))


def test_screen_settings_are_part_of_the_resume_fingerprint() -> None:
    from repgenr.cli.base import _stage_fingerprint

    def fp(**over) -> str:  # noqa: ANN003
        return _stage_fingerprint("assemble", AssembleParams(**over), {}, {})

    fingerprints = {
        fp(),
        fp(screen_reads=True),
        fp(screen_reads=True, screen_min_fraction=0.4),
        fp(screen_reads=True, screen_dup_ani=0.99),
    }
    assert len(fingerprints) == 4


# --- assembly_stats.tsv -----------------------------------------------------------------------


def test_an_assembly_stats_table_without_the_screen_column_is_read(tmp_path) -> None:
    path = tmp_path / ASSEMBLY_STATS_TSV
    write_assembly_stats(path, [AssemblyStatsRow("SRR1", "f.fasta", "skesa", 1, 10, 10, 10)])
    lines = path.read_text(encoding="utf-8").splitlines()
    old = ["\t".join(ln.split("\t")[:-1]) for ln in lines]  # drop the screen column
    path.write_text("\n".join(old) + "\n", encoding="utf-8")
    [row] = read_assembly_stats(path)
    assert row.screen == "" and row.run_accession == "SRR1"


# --- the data-channel step ---------------------------------------------------------------------


def test_assemble_run_screens_without_the_duplicate_check(tmp_path, setup, gtdb) -> None:
    setup.gather.hits.update({"SRR1": (FT, 0.9), "SRR2": (EC, 0.9)})
    reads_tsv = tmp_path / "reads.tsv"
    write_reads(reads_tsv, [setup.row("SRR1"), setup.row("SRR2")])
    setup.ani[("SRR2", "SRR1")] = 1.0
    runs = tmp_path / "runs"
    for r in ("SRR1", "SRR2"):
        assemble_run(
            AssembleRunParams(
                reads_tsv=reads_tsv,
                run=r,
                out_dir=runs / r,
                assembler="fakeasm",
                screen_reads=True,
                gtdb_sketch=gtdb[0],
                gtdb_lineages=gtdb[1],
            ),
            _LOG,
        )
    assert setup.ani.compared == []
    [excuse] = read_excused_runs(runs / "SRR2" / EXCUSED_RUNS_TSV)
    assert excuse.step == "screen" and excuse.reason.startswith("taxon_mismatch")
    assert (runs / "SRR2" / "reads.sig.zip").is_file()
    assert (runs / "SRR2" / "screen.json").is_file()
    out = tmp_path / "out"
    assert reads_gather(ReadsGatherParams(reads_tsv, runs, out), _LOG) == 1
    rows = _screen_rows(out)
    assert rows["SRR1"]["decision"] == "pass" and rows["SRR2"]["decision"] == "excused"
    assert [s.screen for s in read_assembly_stats(out / ASSEMBLY_STATS_TSV)] == ["pass"]


# --- the sketch comparison ---------------------------------------------------------------------


def _write_sig_zip(path: Path, mins: list[int], abund: list[int], ksize: int = 31) -> None:
    import gzip
    import zipfile

    sig = [
        {
            "name": path.stem,
            "signatures": [
                {"ksize": 21, "molecule": "DNA", "mins": [1], "abundances": [5]},
                {"ksize": ksize, "molecule": "DNA", "mins": mins, "abundances": abund},
            ],
        }
    ]
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("signatures/a.sig.gz", gzip.compress(json.dumps(sig).encode()))
        zf.writestr("SOURMASH-MANIFEST.csv", "# SOURMASH-MANIFEST-VERSION: 1.0\n")


def test_solid_cutoff_separates_error_kmers_from_genomic_ones() -> None:
    abund = np.array([1] * 9000 + [2] * 900 + [3] * 100 + [250] * 600)
    assert screen.solid_cutoff(abund) == 25
    assert screen.solid_cutoff(np.array([1, 1, 2, 30, 30, 31])) == 3
    assert screen.solid_cutoff(np.array([], dtype=np.int64)) == 2


def test_containment_ani_of_real_format_sketches(tmp_path) -> None:
    genome = list(range(1000, 1600))
    a, b, c = tmp_path / "a.sig.zip", tmp_path / "b.sig.zip", tmp_path / "c.sig.zip"
    # Two runs of one isolate: the genomic hashes at coverage, different errors.
    _write_sig_zip(a, genome + list(range(5000, 9000)), [80] * 600 + [1] * 4000)
    _write_sig_zip(b, genome + list(range(9000, 12000)), [40] * 600 + [1] * 3000)
    # Another strain: 10 percent of the genomic hashes differ.
    _write_sig_zip(c, genome[:540] + list(range(2000, 2060)), [60] * 600)
    ha, hb, hc = (screen.load_reads_hashes(p) for p in (a, b, c))
    assert ha is not None and hb is not None and hc is not None
    assert ha.solid.size == 600  # the errors are not taken as genomic
    assert screen.containment_ani(ha, hb) == pytest.approx(1.0)
    assert screen.containment_ani(hc, ha) == pytest.approx(0.9 ** (1 / 31))
    assert screen.containment_ani(hc, ha) < 0.999
    assert screen.load_reads_hashes(a, ksize=51) is None
