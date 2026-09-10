"""glance, derep-unpack and derep-stock record a stage (D-7): status lists
them, an identical repeat skips, and the stock listing stays a query."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from repgenr.cli import base as cli
from repgenr.cli.main import app
from repgenr.core.config import Config
from repgenr.core.contracts import CLUSTERS_TSV

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
