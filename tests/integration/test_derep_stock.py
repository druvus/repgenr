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
from repgenr.core.errors import UserInputError
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
    with pytest.raises(UserInputError):
        derep_stock_run(ctx, DerepStockParams(action="delete", name="missing"))
    with pytest.raises(UserInputError):
        derep_stock_run(ctx, DerepStockParams(action="bogus", name="run1"))


@pytest.mark.parametrize("bad_name", ["../escape", "/abs/path", "a/b", "..", ".hidden"])
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
