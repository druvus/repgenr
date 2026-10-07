"""glance stage: dRep-compare ANI overview.

dRep is not available in CI, so the adapter's runner and preflight are
monkeypatched: the fake dRep writes the dendrogram and an Mdb.csv into the work
directory, and the test asserts glance copies the dendrogram out, renders the
boxplot/histogram from Mdb.csv, and honours --keep-files for the scratch dir.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: no display in CI

import repgenr.dereplicators.drep as drep_mod  # noqa: E402
from repgenr.core.context import WorkdirContext  # noqa: E402
from repgenr.stages.glance import GlanceParams  # noqa: E402
from repgenr.stages.glance import run as glance_run  # noqa: E402

_MDB = (
    "genome1,genome2,similarity\n"
    "a.fasta,b.fasta,0.95\n"
    "b.fasta,a.fasta,0.95\n"
    "a.fasta,a.fasta,1.00\n"  # self-comparison, must be skipped
)


def _fake_drep(caps, command, *, logger, **kwargs) -> None:
    """Stand in for `dRep compare`: write the outputs glance reads back."""
    glance_wd = Path(command[-1])
    (glance_wd / "figures").mkdir(parents=True)
    (glance_wd / "figures" / "Primary_clustering_dendrogram.pdf").write_text("%PDF fake")
    (glance_wd / "data_tables").mkdir(parents=True)
    (glance_wd / "data_tables" / "Mdb.csv").write_text(_MDB)


def _setup(workdir: Path) -> WorkdirContext:
    ctx = WorkdirContext(workdir, create=True)
    ctx.genomes_dir.mkdir(parents=True)
    for name in ("a.fasta", "b.fasta"):
        (ctx.genomes_dir / name).write_text(">x\nACGT\n")
    return ctx


def test_glance_happy_path(workdir: Path, monkeypatch) -> None:
    ctx = _setup(workdir)
    monkeypatch.setattr(drep_mod.DrepDereplicator, "preflight", lambda self: {})
    monkeypatch.setattr(drep_mod, "run_tool", _fake_drep)

    out_pdf = glance_run(ctx, GlanceParams(threads=2))

    assert out_pdf.exists()  # dendrogram copied out
    assert (ctx.workdir / "glance_MASH_ANI_similarity_boxplot.png").exists()
    assert (ctx.workdir / "glance_MASH_ANI_similarity_histogram.png").exists()
    assert not (ctx.workdir / "glance_wd").exists()  # scratch cleaned by default


def test_glance_keep_files(workdir: Path, monkeypatch) -> None:
    ctx = _setup(workdir)
    monkeypatch.setattr(drep_mod.DrepDereplicator, "preflight", lambda self: {})
    monkeypatch.setattr(drep_mod, "run_tool", _fake_drep)

    glance_run(ctx, GlanceParams(threads=2, keep_files=True))
    assert (ctx.workdir / "glance_wd").exists()  # scratch retained


def test_glance_passes_fofn(workdir: Path, monkeypatch) -> None:
    ctx = _setup(workdir)
    captured: dict = {}

    def fake(caps, command, *, logger, **kwargs):
        captured["cmd"] = [str(c) for c in command]
        return _fake_drep(caps, command, logger=logger, **kwargs)

    monkeypatch.setattr(drep_mod.DrepDereplicator, "preflight", lambda self: {})
    monkeypatch.setattr(drep_mod, "run_tool", fake)
    glance_run(ctx, GlanceParams())

    parts = captured["cmd"]
    gidx = parts.index("-g")
    assert parts[gidx + 1].endswith(".fofn")
    assert parts[gidx + 2] == "--processors"


def test_glance_records_tool_version(workdir: Path, monkeypatch) -> None:
    # The stage record carries the resolved dRep version, like dereplicate's.
    ctx = _setup(workdir)
    monkeypatch.setattr(drep_mod.DrepDereplicator, "preflight", lambda self: {"dRep": "3.7.1"})
    monkeypatch.setattr(drep_mod, "run_tool", _fake_drep)

    glance_run(ctx, GlanceParams(threads=2))

    record = ctx.config.stages["glance"]
    assert record.tool == "drep"
    assert record.tool_versions == {"dRep": "3.7.1"}


def test_glance_reports_missing_genomes_before_the_tool_check(workdir: Path, monkeypatch) -> None:
    # A workdir without genomes is a workdir error (exit 3) whether or not
    # dRep is installed; the tool check must not mask it.
    import pytest

    from repgenr.core.errors import MissingBinaryError, WorkdirError

    def _absent(self):
        raise MissingBinaryError("dRep: not found on PATH")

    monkeypatch.setattr(drep_mod.DrepDereplicator, "preflight", _absent)
    ctx = WorkdirContext(workdir / "absent")
    with pytest.raises(WorkdirError):
        glance_run(ctx, GlanceParams(threads=2))


def test_glance_removes_plots_of_an_earlier_run(workdir: Path, monkeypatch) -> None:
    ctx = _setup(workdir)
    monkeypatch.setattr(drep_mod.DrepDereplicator, "preflight", lambda self: {})
    monkeypatch.setattr(drep_mod, "run_tool", _fake_drep)
    glance_run(ctx, GlanceParams(threads=2))
    boxplot = ctx.workdir / "glance_MASH_ANI_similarity_boxplot.png"
    histogram = ctx.workdir / "glance_MASH_ANI_similarity_histogram.png"
    assert boxplot.exists() and histogram.exists()
    # No similarity falls in this range, so no plot is drawn this time.
    glance_run(ctx, GlanceParams(threads=2, plot_min=0.99, plot_max=0.999))
    assert not boxplot.exists() and not histogram.exists()


def test_a_drep_failure_prints_one_line_and_logs_the_tail(workdir: Path, monkeypatch) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app
    from repgenr.core.errors import ToolExecutionError

    traceback = "\n".join(
        ["Traceback (most recent call last):"]
        + [f'  File "drep/x.py", line {i}, in f' for i in range(40)]
        + ["OSError: [Errno 22] Invalid argument: '._a.fasta'"]
    )

    def _failing(caps, command, *, logger, **kwargs):
        raise ToolExecutionError(list(map(str, command)), 1, output=traceback, tool="drep")

    _setup(workdir)
    monkeypatch.setattr(drep_mod.DrepDereplicator, "preflight", lambda self: {})
    monkeypatch.setattr(drep_mod, "run_tool", _failing)
    result = CliRunner().invoke(app, ["glance", "-wd", str(workdir)])
    assert result.exit_code == 6, result.output
    assert "Traceback" not in result.output
    errors = [line for line in result.output.splitlines() if "ERROR" in line]
    assert len(errors) == 1 and "drep failed (exit 1)" in errors[0]
    log = (workdir / "repgenr.log").read_text(encoding="utf-8")
    assert "OSError: [Errno 22]" in log and "dRep compare" in log


def test_glance_counts_each_genome_pair_once(tmp_path: Path) -> None:
    # dRep's Mdb.csv lists every pair in both orders; the plots count a pair once.
    from repgenr.stages.glance import _pair_similarities

    mdb = tmp_path / "Mdb.csv"
    mdb.write_text(
        "genome1,genome2,dist,similarity\n"
        "a.fasta,a.fasta,0.0,1.0\n"
        "a.fasta,b.fasta,0.05,0.95\n"
        "b.fasta,a.fasta,0.05,0.95\n"
        "a.fasta,c.fasta,0.10,0.90\n"
        "c.fasta,a.fasta,0.10,0.90\n"
        "b.fasta,c.fasta,0.02,0.98\n"
        "c.fasta,b.fasta,0.02,0.98\n"
        "c.fasta,c.fasta,0.0,1.0\n"
    )
    assert sorted(_pair_similarities(mdb, 0.0, 1.0)) == [0.90, 0.95, 0.98]
    assert sorted(_pair_similarities(mdb, 0.92, 1.0)) == [0.95, 0.98]


def test_glance_histogram_axes_name_ani_and_pair_counts(workdir: Path, monkeypatch) -> None:
    # The histogram's x axis carries the ANI values and its y axis the pair
    # counts; the box plot's single box has no meaningless "1" tick.
    from matplotlib.figure import Figure

    labels: dict[str, tuple[str, str, list[str]]] = {}

    def _capture(self, fname, *args, **kwargs):
        ax = self.axes[0]
        ticks = [t.get_text() for t in ax.get_xticklabels()]
        labels[Path(fname).name] = (ax.get_xlabel(), ax.get_ylabel(), ticks)

    ctx = _setup(workdir)
    monkeypatch.setattr(drep_mod.DrepDereplicator, "preflight", lambda self: {})
    monkeypatch.setattr(drep_mod, "run_tool", _fake_drep)
    monkeypatch.setattr(Figure, "savefig", _capture)
    glance_run(ctx, GlanceParams(threads=2))

    hist_x, hist_y, _ = labels["glance_MASH_ANI_similarity_histogram.png"]
    assert hist_x == "MASH ANI" and hist_y == "Genome pairs"
    box_x, box_y, box_ticks = labels["glance_MASH_ANI_similarity_boxplot.png"]
    assert box_y == "MASH ANI" and box_ticks == [""] * len(box_ticks)


def test_glance_with_one_genome_exits_3_before_the_tool_runs(workdir: Path, monkeypatch) -> None:
    # dRep compare fails inside scipy on a single genome (empty distance
    # matrix); glance names the cause instead of reporting a tool failure.
    import pytest

    from repgenr.core.errors import WorkdirError

    ctx = WorkdirContext(workdir, create=True)
    ctx.genomes_dir.mkdir(parents=True)
    (ctx.genomes_dir / "a.fasta").write_text(">x\nACGT\n")
    called: list[object] = []
    monkeypatch.setattr(drep_mod.DrepDereplicator, "preflight", lambda self: {})
    monkeypatch.setattr(drep_mod, "run_tool", lambda *a, **k: called.append(a))
    with pytest.raises(WorkdirError, match="at least two genomes"):
        glance_run(ctx, GlanceParams(threads=2))
    assert not called


def test_glance_rejects_inverted_or_out_of_range_plot_bounds(workdir: Path, monkeypatch) -> None:
    # Bounds outside 0-1 (for example a percentage) or --plot-min above
    # --plot-max can never select a value; they were accepted and the run
    # removed the plots of the previous run.
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    _setup(workdir)
    called: list[object] = []
    monkeypatch.setattr(drep_mod.DrepDereplicator, "preflight", lambda self: {})
    monkeypatch.setattr(drep_mod, "run_tool", lambda *a, **k: called.append(a))
    for bounds in (
        ["--plot-min", "0.99", "--plot-max", "0.5"],
        ["--plot-max", "99"],
        ["--plot-min", "-0.1"],
    ):
        result = CliRunner().invoke(app, ["glance", "-wd", str(workdir), *bounds])
        assert result.exit_code == 2, (bounds, result.output)
        assert "--plot-m" in result.output
    assert not called


def test_glance_help_names_the_bound_units_and_the_kept_directory() -> None:
    import re

    from typer.testing import CliRunner

    from repgenr.cli.main import app

    result = CliRunner().invoke(app, ["glance", "--help"], terminal_width=200)
    text = " ".join(re.sub(r"\x1b\[[0-9;]*[A-Za-z]|[│╭╮╰╯─]", " ", result.output).split())
    assert "Mash ANI values plotted, as a fraction from 0 to 1" in text
    assert "Keep glance_wd/" in text


def test_glance_warns_when_the_tool_returns_no_dendrogram(workdir: Path, monkeypatch) -> None:
    def _no_dendrogram(caps, command, *, logger, **kwargs) -> None:
        glance_wd = Path(command[-1])
        (glance_wd / "data_tables").mkdir(parents=True)
        (glance_wd / "data_tables" / "Mdb.csv").write_text(_MDB)

    ctx = _setup(workdir)
    monkeypatch.setattr(drep_mod.DrepDereplicator, "preflight", lambda self: {})
    monkeypatch.setattr(drep_mod, "run_tool", _no_dendrogram)
    glance_run(ctx, GlanceParams(threads=2))
    assert not (ctx.workdir / "glance_clustering_dendrogram.pdf").exists()
    log = (ctx.workdir / "repgenr.log").read_text(encoding="utf-8")
    assert "WARNING The comparison returned no dendrogram" in log


def test_a_deleted_dendrogram_is_rebuilt_without_force(workdir: Path, monkeypatch) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    _setup(workdir)
    monkeypatch.setattr(drep_mod.DrepDereplicator, "preflight", lambda self: {})
    monkeypatch.setattr(drep_mod, "run_tool", _fake_drep)
    runner = CliRunner()
    assert runner.invoke(app, ["glance", "-wd", str(workdir), "-t", "2"]).exit_code == 0
    pdf = workdir / "glance_clustering_dendrogram.pdf"
    pdf.unlink()
    result = runner.invoke(app, ["glance", "-wd", str(workdir), "-t", "2"])
    assert result.exit_code == 0, result.output
    assert pdf.exists()
