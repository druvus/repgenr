"""glance stage: dRep-compare ANI overview.

dRep is not available in CI, so the runner and the binary preflight are
monkeypatched: the fake dRep writes the dendrogram and an Mdb.csv into the work
directory, and the test asserts glance copies the dendrogram out, renders the
boxplot/histogram from Mdb.csv, and honours --keep-files for the scratch dir.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: no display in CI

from repgenr.core.context import WorkdirContext  # noqa: E402
from repgenr.stages import glance as glance_mod  # noqa: E402
from repgenr.stages.glance import GlanceParams  # noqa: E402
from repgenr.stages.glance import run as glance_run  # noqa: E402

_MDB = (
    "genome1,genome2,similarity\n"
    "a.fasta,b.fasta,0.95\n"
    "b.fasta,a.fasta,0.80\n"
    "a.fasta,a.fasta,1.00\n"  # self-comparison, must be skipped
)


def _fake_drep(command, *, logger, **kwargs) -> None:
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
    monkeypatch.setattr(glance_mod, "check_binaries", lambda specs: None)
    monkeypatch.setattr(glance_mod, "run_cmd", _fake_drep)

    out_pdf = glance_run(ctx, GlanceParams(threads=2))

    assert out_pdf.exists()  # dendrogram copied out
    assert (ctx.workdir / "glance_MASH_ANI_similarity_boxplot.png").exists()
    assert (ctx.workdir / "glance_MASH_ANI_similarity_histogram.png").exists()
    assert not (ctx.workdir / "glance_wd").exists()  # scratch cleaned by default


def test_glance_keep_files(workdir: Path, monkeypatch) -> None:
    ctx = _setup(workdir)
    monkeypatch.setattr(glance_mod, "check_binaries", lambda specs: None)
    monkeypatch.setattr(glance_mod, "run_cmd", _fake_drep)

    glance_run(ctx, GlanceParams(threads=2, keep_files=True))
    assert (ctx.workdir / "glance_wd").exists()  # scratch retained
