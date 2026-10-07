"""glance, derep-unpack and derep-stock record a stage (D-7): status lists
them, an identical repeat skips, and the stock listing stays a query."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from repgenr.cli import base as cli
from repgenr.cli.main import app
from repgenr.core.config import Config
from repgenr.core.contracts import CLUSTER_SUMMARY_TSV, CLUSTERS_TSV, read_cluster_summary

_runner = CliRunner()


def _derep_workdir(tmp_path: Path) -> Path:
    src = tmp_path / "src"
    src.mkdir()
    for name in ("Fam_Gen_sp1_GCA_000001.1.fasta", "Fam_Gen_sp2_GCA_000002.1.fasta"):
        (src / name).write_text(">s\nACGTACGT\n", encoding="utf-8")
    wd = tmp_path / "wd"
    assert _runner.invoke(app, ["ingest", "-wd", str(wd), "--genomes-dir", str(src)]).exit_code == 0
    derep = wd / "derep"
    (derep / "representatives").mkdir(parents=True)
    for name in ("Fam_Gen_sp1_GCA_000001.1.fasta",):
        (derep / "representatives" / name).write_text(">s\nACGTACGT\n", encoding="utf-8")
    (derep / CLUSTERS_TSV).write_text(
        "representative\tmember\n"
        "Fam_Gen_sp1_GCA_000001.1.fasta\tFam_Gen_sp1_GCA_000001.1.fasta\n"
        "Fam_Gen_sp1_GCA_000001.1.fasta\tFam_Gen_sp2_GCA_000002.1.fasta\n",
        encoding="utf-8",
    )
    (derep / "genome_status.tsv").write_text(
        "genome\tstatus\nFam_Gen_sp1_GCA_000001.1.fasta\trepresentative\n"
        "Fam_Gen_sp2_GCA_000002.1.fasta\tcontained\n",
        encoding="utf-8",
    )
    return wd


def test_derep_unpack_records_and_skips(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    wd = _derep_workdir(tmp_path)
    first = _runner.invoke(app, ["derep-unpack", "-wd", str(wd)])
    assert first.exit_code == 0, first.output
    rec = Config.load(wd).stages["derep_unpack"]
    assert rec.completed and rec.inputs, "recorded with digested inputs"
    stamp = rec.completed
    second = _runner.invoke(app, ["derep-unpack", "-wd", str(wd)])
    assert second.exit_code == 0 and Config.load(wd).stages["derep_unpack"].completed == stamp
    status = _runner.invoke(app, ["status", "-wd", str(wd)]).stdout
    assert "optional stages run" in status and "derep_unpack" in status


def test_cluster_summary_regenerates_from_clusters_tsv(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    wd = _derep_workdir(tmp_path)
    first = _runner.invoke(app, ["cluster-summary", "-wd", str(wd)])
    assert first.exit_code == 0, first.output
    (row,) = read_cluster_summary(wd / "derep" / CLUSTER_SUMMARY_TSV)
    assert row.representative == "Fam_Gen_sp1_GCA_000001.1.fasta"
    assert (row.n_members, row.n_species, row.species) == (1, 2, "sp1,sp2")
    assert row.rep_completeness is None and row.best_member == ""
    rec = Config.load(wd).stages["cluster_summary"]
    assert rec.completed and rec.inputs
    stamp = rec.completed
    second = _runner.invoke(app, ["cluster-summary", "-wd", str(wd)])
    assert second.exit_code == 0 and Config.load(wd).stages["cluster_summary"].completed == stamp


def test_derep_stock_list_is_a_query_but_pack_is_recorded(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    wd = _derep_workdir(tmp_path)
    listed = _runner.invoke(app, ["derep-stock", "-wd", str(wd), "--action", "list"])
    assert listed.exit_code == 0 and "derep_stock" not in Config.load(wd).stages
    packed = _runner.invoke(
        app, ["derep-stock", "-wd", str(wd), "--action", "pack", "--name", "r1"]
    )
    assert packed.exit_code == 0, packed.output
    rec = Config.load(wd).stages["derep_stock"]
    assert rec.completed and rec.params["action"] == "pack" and rec.params["name"] == "r1"
    assert (wd / "derep" / "stock" / "r1").is_dir()


def test_cluster_summary_on_missing_workdir_exits_3_without_creating_it(tmp_path: Path) -> None:
    # The manifest digest of the resume fingerprint must not create the
    # workdir (or its manifest) before the stage reports the missing input.
    wd = tmp_path / "absent"
    result = _runner.invoke(app, ["cluster-summary", "-wd", str(wd)])
    assert result.exit_code == 3, result.output
    assert not wd.exists()


def test_derep_stock_unpack_reruns_after_the_dereplication_changed(
    tmp_path: Path, monkeypatch
) -> None:
    # A repeat unpack of the same run must restore it when the live derep/
    # outputs changed in between (e.g. a new dereplicate), and skip otherwise.
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    wd = _derep_workdir(tmp_path)
    clusters = wd / "derep" / CLUSTERS_TSV
    stored = clusters.read_text(encoding="utf-8")

    def stock(action: str) -> None:
        result = _runner.invoke(
            app, ["derep-stock", "-wd", str(wd), "--action", action, "--name", "r1"]
        )
        assert result.exit_code == 0, result.output

    stock("pack")
    stock("unpack")
    stamp = Config.load(wd).stages["derep_stock"].completed
    stock("unpack")  # nothing changed: skipped
    assert Config.load(wd).stages["derep_stock"].completed == stamp

    # A different dereplication replaces the live outputs.
    clusters.write_text(
        "representative\tmember\n"
        "Fam_Gen_sp2_GCA_000002.1.fasta\tFam_Gen_sp2_GCA_000002.1.fasta\n"
        "Fam_Gen_sp2_GCA_000002.1.fasta\tFam_Gen_sp1_GCA_000001.1.fasta\n",
        encoding="utf-8",
    )
    stock("unpack")
    assert clusters.read_text(encoding="utf-8") == stored


def test_derep_stock_on_missing_workdir_exits_3(tmp_path: Path) -> None:
    # A mistyped -wd must not read as an empty store.
    wd = tmp_path / "absent"
    for args in (["--action", "list"], ["--action", "delete", "--name", "r1"]):
        result = _runner.invoke(app, ["derep-stock", "-wd", str(wd), *args])
        assert result.exit_code == 3, result.output
    assert not wd.exists()


def test_derep_stock_list_prints_the_runs_on_stdout_under_quiet(
    tmp_path: Path, monkeypatch
) -> None:
    # The list is the command's result, not a log message: it must reach
    # stdout, one name per line, also with --quiet. Delete says what it removed.
    import logging

    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    monkeypatch.setitem(cli._RUN_STATE, "log_level", logging.INFO)
    wd = _derep_workdir(tmp_path)
    args = ["derep-stock", "-wd", str(wd)]
    for name in ("r2", "r1", "r3"):
        assert _runner.invoke(app, [*args, "--action", "pack", "--name", name]).exit_code == 0
    deleted = _runner.invoke(app, [*args, "--action", "delete", "--name", "r3"])
    assert deleted.exit_code == 0 and "Deleted stored run 'r3'" in deleted.output
    listed = _runner.invoke(app, ["--quiet", *args, "--action", "list"])
    assert listed.exit_code == 0, listed.output
    assert listed.stdout.splitlines() == ["r1", "r2"]


def test_a_refused_derep_stock_call_leaves_the_last_record_complete(
    tmp_path: Path, monkeypatch
) -> None:
    # One derep_stock record serves every stored run, so a mistyped name or
    # an unknown run must be refused before the harness marks the record of
    # the last (finished) pack as interrupted; doctor then reports no failure.
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    wd = _derep_workdir(tmp_path)
    args = ["derep-stock", "-wd", str(wd)]
    assert _runner.invoke(app, [*args, "--action", "pack", "--name", "r1"]).exit_code == 0
    stamp = Config.load(wd).stages["derep_stock"].completed
    refusals = [
        (["--action", "pack", "--name", "a/b"], 2),
        (["--action", "unpack", "--name", "a/b"], 2),
        (["--action", "unpack", "--name", "missing"], 2),
    ]
    for extra, code in refusals:
        result = _runner.invoke(app, [*args, *extra])
        assert result.exit_code == code, (extra, result.output)
        assert Config.load(wd).stages["derep_stock"].completed == stamp, extra
    (wd / "derep" / CLUSTERS_TSV).unlink()
    result = _runner.invoke(app, [*args, "--action", "pack", "--name", "r2"])
    assert result.exit_code == 3, result.output
    assert Config.load(wd).stages["derep_stock"].completed == stamp


def test_derep_stock_pack_over_a_stored_name_warns(tmp_path: Path, monkeypatch) -> None:
    # Packing under a name already in the store replaces that run; say so.
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    wd = _derep_workdir(tmp_path)
    args = ["derep-stock", "-wd", str(wd), "--action", "pack", "--name", "r1"]
    first = _runner.invoke(app, args)
    assert first.exit_code == 0 and "Replacing" not in first.output
    (wd / "derep" / CLUSTERS_TSV).write_text(
        "representative\tmember\nFam_Gen_sp1_GCA_000001.1.fasta\tFam_Gen_sp1_GCA_000001.1.fasta\n",
        encoding="utf-8",
    )
    second = _runner.invoke(app, args)
    assert second.exit_code == 0, second.output
    assert "Replacing stored run 'r1'" in second.output
