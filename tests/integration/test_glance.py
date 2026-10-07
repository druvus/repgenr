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
    "b.fasta,a.fasta,0.80\n"
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
