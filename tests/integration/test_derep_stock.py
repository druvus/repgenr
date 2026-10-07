"""derep_stock stage: store/load named dereplication runs.

Pure file-management logic (no external binaries): exercises list/pack/unpack/
delete against a temp workdir holding a derep/ contract, asserting a pack ->
unpack round-trip restores the representatives and the flat contract files, and
that the error paths raise.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import write_clusters, write_genome_status
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.stages.derep_stock import DerepStockParams
from repgenr.stages.derep_stock import run as derep_stock_run

_GENOMES = [
    "Fam_Gen_sp_GCA_000001.1.fasta",
    "Fam_Gen_sp_GCA_000002.1.fasta",
    "Fam_Gen_sp_GCA_000003.1.fasta",
]
_REPS = _GENOMES[:2]  # two of the three are representatives


def _setup_contract(workdir: Path) -> WorkdirContext:
    ctx = WorkdirContext(workdir, create=True)
    ctx.genomes_dir.mkdir(parents=True)
    for name in _GENOMES:
        (ctx.genomes_dir / name).write_text(">x\nACGT\n")
    ctx.representatives_dir.mkdir(parents=True)
    for name in _REPS:
        (ctx.representatives_dir / name).write_text(">x\nACGT\n")
    write_clusters(ctx.derep_dir / "clusters.tsv", {_REPS[0]: [_GENOMES[2]], _REPS[1]: []})
    write_genome_status(
        ctx.derep_dir / "genome_status.tsv",
        {_REPS[0]: "representative", _REPS[1]: "representative", _GENOMES[2]: "contained"},
    )
    return ctx


def test_pack_unpack_round_trip(workdir: Path) -> None:
    ctx = _setup_contract(workdir)
    store = ctx.derep_dir / "stock"

    derep_stock_run(ctx, DerepStockParams(action="pack", name="run1"))
    packed = store / "run1"
    assert (packed / "clusters.tsv").exists()
    assert (packed / "genome_status.tsv").exists()
    assert {p.name for p in (packed / "representatives").iterdir()} == set(_REPS)

    # Wipe the live representatives + flat files, then unpack to restore them.
    for f in ("clusters.tsv", "genome_status.tsv"):
        (ctx.derep_dir / f).unlink()
    for rep in list(ctx.representatives_dir.iterdir()):
        rep.unlink()

    derep_stock_run(ctx, DerepStockParams(action="unpack", name="run1"))
    assert {p.name for p in ctx.representatives_dir.iterdir()} == set(_REPS)
    assert (ctx.derep_dir / "clusters.tsv").exists()
    assert (ctx.derep_dir / "genome_status.tsv").exists()


def test_list_and_delete(workdir: Path) -> None:
    ctx = _setup_contract(workdir)
    derep_stock_run(ctx, DerepStockParams(action="pack", name="run1"))
    # list does not raise whether or not runs exist
    derep_stock_run(ctx, DerepStockParams(action="list"))
    derep_stock_run(ctx, DerepStockParams(action="delete", name="run1"))
    assert not (ctx.derep_dir / "stock" / "run1").exists()
    derep_stock_run(ctx, DerepStockParams(action="list"))  # empty store, still fine


def test_error_paths(workdir: Path) -> None:
    ctx = _setup_contract(workdir)
    with pytest.raises(UserInputError):
        derep_stock_run(ctx, DerepStockParams(action="pack", name=None))  # name required
    with pytest.raises(UserInputError):
        derep_stock_run(ctx, DerepStockParams(action="unpack", name="missing"))
    with pytest.raises(WorkdirError, match="'missing'; stored runs: none"):
        derep_stock_run(ctx, DerepStockParams(action="delete", name="missing"))
    with pytest.raises(UserInputError):
        derep_stock_run(ctx, DerepStockParams(action="bogus", name="run1"))


@pytest.mark.parametrize(
    "bad_name",
    [
        "../escape",
        "/abs/path",
        "a/b",
        "..",
        ".hidden",
        "-x",
        "with space",
        pytest.param("x" * 300, id="too-long"),
    ],
)
def test_traversal_names_rejected(workdir: Path, bad_name: str) -> None:
    # --name becomes a directory under the stock store and is passed to rmtree
    # on delete, so separators and dot-prefixes must be rejected outright.
    ctx = _setup_contract(workdir)
    for action in ("pack", "unpack", "delete"):
        with pytest.raises(UserInputError, match="name"):
            derep_stock_run(ctx, DerepStockParams(action=action, name=bad_name))


def test_unpack_restores_summary_record_and_manifest(workdir: Path) -> None:
    """Unpacking a stored run leaves the workdir describing that run, not the
    dereplication that produced the record before it."""
    from repgenr.core.manifest import GenomeRecord

    ctx = _setup_contract(workdir)
    (ctx.derep_dir / "cluster_summary.tsv").write_text("representative\tn_members\n", "utf-8")
    ctx.manifest.replace_genomes(
        [
            GenomeRecord(accession=f"GCA_00000{i}.1", filename=name)
            for i, name in enumerate(_GENOMES, start=1)
        ]
    )
    ctx.config.record_stage(
        "dereplicate", tool="skder", params={"tool": "skder"}, completed="t0", fingerprint="fp0"
    )
    ctx.save_config()

    derep_stock_run(ctx, DerepStockParams(action="pack", name="run1"))
    (ctx.derep_dir / "cluster_summary.tsv").unlink()
    ctx.manifest.set_derep_status_many([("GCA_000001.1", "contained", "GCA_000003.1")])
    derep_stock_run(ctx, DerepStockParams(action="unpack", name="run1"))

    assert (ctx.derep_dir / "cluster_summary.tsv").exists()
    record = ctx.config.stages["dereplicate"]
    assert record.fingerprint is None  # the next `dereplicate` must not skip
    assert record.completed not in (None, "t0")
    assert record.params.get("stock") == "run1"
    status = {g.accession: g.derep_status for g in ctx.manifest.all_genomes()}
    assert status["GCA_000001.1"] == "representative"
    assert status["GCA_000003.1"] == "contained"


@pytest.mark.parametrize("missing", ["clusters.tsv", "representatives"])
def test_pack_without_dereplication_outputs_raises(workdir: Path, missing: str) -> None:
    # Packing a workdir that holds no dereplication must not store an empty run.
    import shutil

    ctx = _setup_contract(workdir)
    target = ctx.derep_dir / missing
    if target.is_dir():
        shutil.rmtree(target)
    else:
        target.unlink()
    with pytest.raises(WorkdirError, match=missing):
        derep_stock_run(ctx, DerepStockParams(action="pack", name="run1"))
    assert not (ctx.derep_dir / "stock" / "run1").exists()


@pytest.mark.parametrize("damage", ["stored_representatives", "stored_clusters", "genome"])
def test_unpack_of_incomplete_run_raises_before_changing_the_workdir(
    workdir: Path, damage: str
) -> None:
    # A stored run that cannot be restored in full must leave the current
    # dereplication in place and report the problem as a workdir error.
    import shutil

    ctx = _setup_contract(workdir)
    derep_stock_run(ctx, DerepStockParams(action="pack", name="run1"))
    packed = ctx.derep_dir / "stock" / "run1"
    if damage == "stored_representatives":
        shutil.rmtree(packed / "representatives")
    elif damage == "stored_clusters":
        (packed / "clusters.tsv").unlink()
    else:
        (ctx.genomes_dir / _REPS[0]).unlink()
    clusters_before = (ctx.derep_dir / "clusters.tsv").read_text()

    with pytest.raises(WorkdirError):
        derep_stock_run(ctx, DerepStockParams(action="unpack", name="run1"))
    assert {p.name for p in ctx.representatives_dir.iterdir()} == set(_REPS)
    assert (ctx.derep_dir / "clusters.tsv").read_text() == clusters_before


def test_unpack_of_a_run_without_summary_rebuilds_the_summary(workdir: Path) -> None:
    from repgenr.core.contracts import read_cluster_summary

    ctx = _setup_contract(workdir)
    derep_stock_run(ctx, DerepStockParams(action="pack", name="old"))  # no summary stored
    # A later dereplication with a different clustering writes its own summary.
    (ctx.derep_dir / "cluster_summary.tsv").write_text(
        "representative\tn_members\nstale.fasta\t9\n", "utf-8"
    )
    derep_stock_run(ctx, DerepStockParams(action="unpack", name="old"))
    rows = read_cluster_summary(ctx.derep_dir / "cluster_summary.tsv")
    assert {r.representative for r in rows} == set(_REPS)


def test_unpack_rebuilds_the_summary_from_the_live_manifest(workdir: Path) -> None:
    # The stored summary carries the quality of pack time. After unpack the
    # live derep/cluster_summary.tsv must reflect the manifest as it is now,
    # as `cluster-summary` would write it (it skips, since its inputs match).
    from repgenr.core.contracts import read_cluster_summary
    from repgenr.core.manifest import GenomeRecord
    from repgenr.stages.cluster_summary import ClusterSummaryParams
    from repgenr.stages.cluster_summary import run as cluster_summary_run

    ctx = _setup_contract(workdir)

    def set_quality(completeness: float) -> None:
        ctx.manifest.replace_genomes(
            [
                GenomeRecord(
                    accession=f"GCA_00000{i}.1",
                    filename=name,
                    completeness=completeness,
                    contamination=1.0,
                )
                for i, name in enumerate(_GENOMES, start=1)
            ]
        )

    set_quality(97.5)
    cluster_summary_run(ctx, ClusterSummaryParams())
    derep_stock_run(ctx, DerepStockParams(action="pack", name="run1"))
    set_quality(50.0)
    derep_stock_run(ctx, DerepStockParams(action="unpack", name="run1"))
    rows = read_cluster_summary(ctx.derep_dir / "cluster_summary.tsv")
    assert {r.rep_completeness for r in rows} == {50.0}
    # The stored copy keeps the pack-time view.
    stored = read_cluster_summary(ctx.derep_dir / "stock" / "run1" / "cluster_summary.tsv")
    assert {r.rep_completeness for r in stored} == {97.5}


def test_deleting_an_unknown_run_exits_3_and_lists_the_stored_runs(workdir: Path) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    _setup_contract(workdir)
    runner = CliRunner()
    args = ["derep-stock", "-wd", str(workdir)]
    assert runner.invoke(app, [*args, "--action", "pack", "--name", "keep"]).exit_code == 0
    assert runner.invoke(app, [*args, "--action", "pack", "--name", "gone"]).exit_code == 0
    assert runner.invoke(app, [*args, "--action", "delete", "--name", "gone"]).exit_code == 0
    again = runner.invoke(app, [*args, "--action", "delete", "--name", "gone"])
    assert again.exit_code == 3, again.output
    assert "'gone'" in again.output and "keep" in again.output


def test_unpack_ignores_an_incomplete_dereplicate_record(workdir: Path) -> None:
    """An [interrupted] dereplicate record describes a run that did not finish;
    unpacking must not copy its tool or parameters onto the restored run."""
    ctx = _setup_contract(workdir)
    derep_stock_run(ctx, DerepStockParams(action="pack", name="run1"))
    ctx.config.record_stage("dereplicate", tool="vsearch", params={"tool": "vsearch"})
    ctx.save_config()

    derep_stock_run(ctx, DerepStockParams(action="unpack", name="run1"))
    record = ctx.config.stages["dereplicate"]
    assert record.completed
    assert record.tool is None
    assert record.params == {"stock": "run1"}


def test_an_interrupted_unpack_leaves_the_dereplicate_record_incomplete(
    workdir: Path, monkeypatch
) -> None:
    """Unpack replaces the dereplicate stage's outputs; if it stops half-way,
    `status` must not keep reporting the replaced dereplication as done."""
    from repgenr.core.config import Config
    from repgenr.stages import derep_stock

    ctx = _setup_contract(workdir)
    ctx.config.record_stage(
        "dereplicate", tool="skder", params={"tool": "skder"}, completed="t0", fingerprint="fp0"
    )
    ctx.save_config()
    derep_stock_run(ctx, DerepStockParams(action="pack", name="run1"))

    def killed(path):  # noqa: ANN001
        raise KeyboardInterrupt

    monkeypatch.setattr(derep_stock, "remove_tree", killed)
    with pytest.raises(KeyboardInterrupt):
        derep_stock_run(ctx, DerepStockParams(action="unpack", name="run1"))
    record = Config.load(workdir).stages["dereplicate"]
    assert record.completed is None and record.fingerprint is None

    monkeypatch.undo()
    derep_stock_run(ctx, DerepStockParams(action="unpack", name="run1"))
    record = Config.load(workdir).stages["dereplicate"]
    # The repeat keeps what the interrupted unpack was carrying over.
    assert record.completed and record.tool == "skder"
    assert record.params == {"tool": "skder", "stock": "run1"}


def test_unpack_links_the_representatives_like_dereplicate(workdir: Path) -> None:
    # dereplicate hardlinks representatives to genomes/ where it can; unpack
    # does the same instead of copying every genome again.
    ctx = _setup_contract(workdir)
    derep_stock_run(ctx, DerepStockParams(action="pack", name="run1"))
    derep_stock_run(ctx, DerepStockParams(action="unpack", name="run1"))
    for name in _REPS:
        restored = (ctx.representatives_dir / name).stat()
        assert restored.st_ino == (ctx.genomes_dir / name).stat().st_ino
