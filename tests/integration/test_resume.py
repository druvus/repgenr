"""Resume/idempotency: a completed stage with unchanged params is skipped."""

from __future__ import annotations

import sys
import types
from dataclasses import dataclass
from pathlib import Path

from repgenr.cli import base as cli


@dataclass
class _P:
    a: int = 1


def _install_fake_stage(monkeypatch, calls: list[int]) -> None:
    fake = types.ModuleType("repgenr.stages.faketest")

    def run(ctx, params):  # noqa: ANN001
        calls.append(params.a)
        ctx.config.record_stage("faketest", tool="x", params={"a": params.a}, completed="t")
        ctx.save_config()

    fake.run = run  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "repgenr.stages.faketest", fake)


def test_resume_skips_then_force_and_param_change_rerun(tmp_path: Path, monkeypatch) -> None:
    calls: list[int] = []
    _install_fake_stage(monkeypatch, calls)

    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    cli._run("faketest", tmp_path, lambda: _P(), create=True)  # runs
    cli._run("faketest", tmp_path, lambda: _P(), create=True)  # skipped (same params)
    assert calls == [1]

    monkeypatch.setitem(cli._RUN_STATE, "force", True)  # --force re-runs
    cli._run("faketest", tmp_path, lambda: _P(), create=True)
    assert calls == [1, 1]

    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    cli._run("faketest", tmp_path, lambda: _P(a=2), create=True)  # changed params re-runs
    assert calls == [1, 1, 2]


def test_incomplete_stage_reruns(tmp_path: Path, monkeypatch) -> None:
    # A stage that never recorded `completed` (e.g. crashed) must re-run.
    calls: list[int] = []
    fake = types.ModuleType("repgenr.stages.faketest2")

    def run(ctx, params):  # noqa: ANN001
        calls.append(params.a)
        ctx.config.record_stage("faketest2", tool="x", params={})  # no completed stamp
        ctx.save_config()

    fake.run = run  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "repgenr.stages.faketest2", fake)
    monkeypatch.setitem(cli._RUN_STATE, "force", False)

    cli._run("faketest2", tmp_path, lambda: _P(), create=True)
    cli._run("faketest2", tmp_path, lambda: _P(), create=True)
    assert calls == [1, 1]  # not skipped, because completed was never set


def test_registered_input_change_reruns(tmp_path: Path, monkeypatch) -> None:
    """A stage with a STAGE_INPUTS entry reruns when its input dir changes."""
    calls: list[int] = []
    _install_fake_stage(monkeypatch, calls)
    input_dir = tmp_path / "genomes"
    input_dir.mkdir(parents=True)
    (input_dir / "g1.fasta").write_text(">g1\nACGT\n", encoding="utf-8")
    monkeypatch.setitem(cli.STAGE_INPUTS, "faketest", lambda ctx, p: [ctx.genomes_dir])
    monkeypatch.setitem(cli._RUN_STATE, "force", False)

    cli._run("faketest", tmp_path, lambda: _P(), create=True)  # runs
    cli._run("faketest", tmp_path, lambda: _P(), create=True)  # unchanged input -> skip
    assert calls == [1]

    (input_dir / "g2.fasta").write_text(">g2\nACGT\n", encoding="utf-8")
    cli._run("faketest", tmp_path, lambda: _P(), create=True)  # new input file -> rerun
    assert calls == [1, 1]

    cli._run("faketest", tmp_path, lambda: _P(), create=True)  # stable again -> skip
    assert calls == [1, 1]


def test_unregistered_stage_still_skips_on_params(tmp_path: Path, monkeypatch) -> None:
    """A stage without a STAGE_INPUTS entry fingerprints on params alone."""
    calls: list[int] = []
    _install_fake_stage(monkeypatch, calls)
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    assert "faketest" not in cli.STAGE_INPUTS

    cli._run("faketest", tmp_path, lambda: _P(), create=True)
    cli._run("faketest", tmp_path, lambda: _P(), create=True)
    assert calls == [1]


# -- deliverable checks: a completed stage whose outputs were deleted reruns --


def _flush_log(workdir: Path) -> str:
    import logging

    for h in logging.getLogger("repgenr").handlers:
        h.flush()
    return (workdir / "repgenr.log").read_text(encoding="utf-8")


def test_deliverable_directory_must_be_non_empty(tmp_path: Path) -> None:
    from repgenr.core.context import WorkdirContext

    ctx = WorkdirContext(tmp_path, create=True)
    try:
        ctx.genomes_dir.mkdir()
        (tmp_path / "selection.tsv").write_text("x\n", encoding="utf-8")
        (tmp_path / "manifest.sqlite").write_text("", encoding="utf-8")
        params = types.SimpleNamespace(genomes_dir=str(tmp_path / "src"))
        missing = cli.missing_deliverables(ctx, "ingest", params)
        assert missing == [ctx.genomes_dir]
        (ctx.genomes_dir / "a.fasta").write_text(">a\nA\n", encoding="utf-8")
        assert cli.missing_deliverables(ctx, "ingest", params) == []
    finally:
        ctx.close()


def _assert_rerun_on_missing_deliverable(
    workdir: Path, stage: str, build_params, deliverable: Path, run_marker: str
) -> None:
    import shutil

    cli._run(stage, workdir, build_params, create=True)  # cold run
    log = _flush_log(workdir)
    assert log.count(run_marker) == 1
    assert deliverable.exists()

    cli._run(stage, workdir, build_params, create=True)  # identical -> skip
    log = _flush_log(workdir)
    assert log.count("already completed") == 1
    assert log.count(run_marker) == 1

    if deliverable.is_dir():
        shutil.rmtree(deliverable)
    else:
        deliverable.unlink()
    cli._run(stage, workdir, build_params, create=True)  # deliverable gone -> rerun
    log = _flush_log(workdir)
    rel = deliverable.relative_to(workdir)
    assert f"Stage '{stage}': deliverable {rel} missing; re-running." in log
    assert log.count(run_marker) == 2
    assert log.count("already completed") == 1
    assert deliverable.exists()

    cli._run(stage, workdir, build_params, create=True)  # untouched -> skip again
    log = _flush_log(workdir)
    assert log.count("already completed") == 2
    assert log.count(run_marker) == 2


def _resume_env(monkeypatch) -> None:
    import logging

    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    monkeypatch.setitem(cli._RUN_STATE, "log_level", logging.INFO)


def test_ingest_reruns_when_genomes_dir_deleted(tmp_path: Path, monkeypatch) -> None:
    from repgenr.stages.ingest import IngestParams

    _resume_env(monkeypatch)
    src = tmp_path / "src"
    src.mkdir()
    for i in range(2):
        (src / f"Fam_Gen_sp_GCA_00000{i}.1.fasta").write_text(">s\nACGTACGT\n")
    wd = tmp_path / "wd"
    _assert_rerun_on_missing_deliverable(
        wd,
        "ingest",
        lambda: IngestParams(genomes_dir=str(src)),
        wd / "genomes",
        "Ingested 2",
    )


def test_dereplicate_reruns_when_clusters_tsv_deleted(tmp_path: Path, monkeypatch) -> None:
    from repgenr.core.context import WorkdirContext
    from repgenr.core.manifest import GenomeRecord
    from repgenr.core.plugins import ToolCapabilities
    from repgenr.dereplicators.base import (
        STATUS_REPRESENTATIVE,
        Dereplicator,
        DerepResult,
        registry,
    )
    from repgenr.stages.dereplicate import DereplicateParams

    class _NoRep(Dereplicator):
        capabilities = ToolCapabilities(name="deliverablenorep", supports_native_scaling=True)

        def preflight(self) -> dict[str, str]:
            return {"deliverablenorep": "1.0"}

        def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
            genomes = list(genomes)
            return DerepResult(
                representatives=list(genomes),
                clusters={g.name: [] for g in genomes},
                genome_status={g.name: STATUS_REPRESENTATIVE for g in genomes},
            )

    _resume_env(monkeypatch)
    registry._load()
    registry.register("deliverablenorep", _NoRep, replace=True)
    try:
        gdir = tmp_path / "genomes"
        gdir.mkdir()
        for i in range(2):
            (gdir / f"Fam_g_s_GCF_10000{i}.1.fasta").write_text(">x\nACGT\n", encoding="utf-8")
        seed = WorkdirContext(tmp_path, create=True)
        seed.manifest.upsert_many(
            [
                GenomeRecord(accession=f"GCF_10000{i}.1", filename=f"Fam_g_s_GCF_10000{i}.1.fasta")
                for i in range(2)
            ]
        )
        seed.close()
        _assert_rerun_on_missing_deliverable(
            tmp_path,
            "dereplicate",
            lambda: DereplicateParams(tool="deliverablenorep"),
            tmp_path / "derep" / "clusters.tsv",
            "Dereplicating 2 genomes",
        )
    finally:
        registry._classes.pop("deliverablenorep", None)


def test_tree2tax_reruns_when_genomes_map_deleted(tmp_path: Path, monkeypatch) -> None:
    from repgenr.stages.tree2tax import Tree2taxParams

    _resume_env(monkeypatch)
    tree_dir = tmp_path / "tree"
    tree_dir.mkdir()
    (tree_dir / "tree.nwk").write_text(
        "((Fam_gen_sp_GCA_000001:0.1,Fam_gen_sp_GCA_000002:0.1):0.2,Fam_gen_sp_GCA_000003:0.5);\n"
    )
    _assert_rerun_on_missing_deliverable(
        tmp_path,
        "tree2tax",
        lambda: Tree2taxParams(include_dereplicated=False),
        tmp_path / "genomes_map.tsv",
        "Wrote ",
    )


def test_deliverable_directory_with_only_dotfiles_is_missing(tmp_path: Path) -> None:
    # Finder leaves .DS_Store in a directory it showed; an emptied
    # representatives/ holding only that must not let dereplicate skip.
    from repgenr.core.context import WorkdirContext

    ctx = WorkdirContext(tmp_path, create=True)
    try:
        ctx.genomes_dir.mkdir()
        (ctx.genomes_dir / ".DS_Store").write_bytes(b"")
        (ctx.genomes_dir / "._a.fasta").write_bytes(b"")
        (tmp_path / "selection.tsv").write_text("x\n", encoding="utf-8")
        (tmp_path / "manifest.sqlite").write_text("", encoding="utf-8")
        params = types.SimpleNamespace(genomes_dir=str(tmp_path / "src"))
        assert cli.missing_deliverables(ctx, "ingest", params) == [ctx.genomes_dir]
    finally:
        ctx.close()


def test_dereplicate_deliverables_list_each_representative(tmp_path: Path) -> None:
    # A representative deleted by hand reruns dereplicate; it used to skip
    # while doctor asked for a rerun and phylo refused the directory.
    from repgenr.core.context import WorkdirContext
    from repgenr.core.contracts import write_clusters

    ctx = WorkdirContext(tmp_path, create=True)
    try:
        ctx.representatives_dir.mkdir(parents=True)
        write_clusters(ctx.derep_dir / "clusters.tsv", {"a.fasta": [], "b.fasta": ["c.fasta"]})
        for name in ("genome_status.tsv", "cluster_summary.tsv"):
            (ctx.derep_dir / name).write_text("x\n", encoding="utf-8")
        (ctx.representatives_dir / "a.fasta").write_text(">a\nA\n", encoding="utf-8")
        params = types.SimpleNamespace()
        assert cli.missing_deliverables(ctx, "dereplicate", params) == [
            ctx.representatives_dir / "b.fasta"
        ]
    finally:
        ctx.close()
