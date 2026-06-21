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
