"""The `repgenr status` command reports pipeline progress from repgenr.yaml."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core.config import Config

_runner = CliRunner()


def test_status_missing_workdir_exits_3(tmp_path: Path) -> None:
    missing = tmp_path / "missing"
    result = _runner.invoke(app, ["status", "-wd", str(missing)])
    assert result.exit_code == 3
    assert "Workdir not found" in result.output
    assert not missing.exists()


def test_status_empty_workdir_hints_at_entry_stages(tmp_path: Path) -> None:
    result = _runner.invoke(app, ["status", "-wd", str(tmp_path)])
    assert result.exit_code == 0
    assert "No RepGenR run found" in result.stdout


def test_doctor_missing_workdir_exits_3(tmp_path: Path) -> None:
    missing = tmp_path / "missing"
    result = _runner.invoke(app, ["doctor", "-wd", str(missing)])
    assert result.exit_code == 3
    assert "Workdir not found" in result.output
    assert not missing.exists()


def test_doctor_empty_workdir_warns(tmp_path: Path) -> None:
    result = _runner.invoke(app, ["doctor", "-wd", str(tmp_path)])
    assert result.exit_code == 0
    assert "No RepGenR run found" in result.stdout


def test_status_bacterial_progress(tmp_path: Path, write_deliverables) -> None:
    cfg = Config()
    cfg.record_stage("metadata", completed="2026-01-01T00:00:00")
    cfg.record_stage("genome", completed="2026-01-01T00:01:00")
    cfg.save(tmp_path)
    write_deliverables(tmp_path)

    result = _runner.invoke(app, ["status", "-wd", str(tmp_path)])
    assert result.exit_code == 0
    assert "Pipeline: bacterial" in result.stdout
    assert "[done]    metadata" in result.stdout
    assert "[done]    genome" in result.stdout
    # dereplicate is the first incomplete stage
    assert "[next] dereplicate" in result.stdout
    assert "Next: repgenr dereplicate" in result.stdout


def test_status_viral_detection_and_extras(tmp_path: Path) -> None:
    cfg = Config()
    cfg.record_stage("vmetadata", completed="2026-01-01T00:00:00")
    cfg.record_stage("vgenome", completed="2026-01-01T00:01:00")
    cfg.record_stage("snptype", tool="simple", completed="2026-01-01T00:02:00")
    cfg.save(tmp_path)

    result = _runner.invoke(app, ["status", "-wd", str(tmp_path)])
    assert result.exit_code == 0
    assert "Pipeline: viral" in result.stdout
    assert "optional stages run:" in result.stdout
    assert "snptype" in result.stdout


def test_status_all_complete(tmp_path: Path, write_deliverables) -> None:
    cfg = Config()
    for stage in ("metadata", "genome", "dereplicate", "phylo", "tree2tax"):
        cfg.record_stage(stage, completed="2026-01-01T00:00:00")
    cfg.save(tmp_path)
    write_deliverables(tmp_path)

    result = _runner.invoke(app, ["status", "-wd", str(tmp_path)])
    assert result.exit_code == 0
    assert "All stages complete" in result.stdout


def test_status_marks_interrupted_stage(tmp_path: Path) -> None:
    cfg = Config()
    cfg.record_stage("metadata", completed="2026-01-01T00:00:00")
    # dereplicate started a re-run and crashed: provenance present, no stamp
    cfg.record_stage("genome", tool="datasets", params={"total": 5})
    cfg.save(tmp_path)

    result = _runner.invoke(app, ["status", "-wd", str(tmp_path)])
    assert result.exit_code == 0
    assert "[interrupted] genome" in result.stdout
    assert "did not finish; outputs may be partial" in result.stdout


def test_status_marks_an_interrupted_optional_stage_like_a_chain_stage(tmp_path: Path) -> None:
    # A glance run killed mid-way leaves a record without a stamp; status
    # names it interrupted, as it does for the stages of the chain.
    cfg = Config()
    cfg.record_stage("ingest", completed="2026-01-01T00:00:00")
    cfg.record_stage("glance", tool="drep", params={"plot_min": 0.0})
    cfg.save(tmp_path)

    result = _runner.invoke(app, ["status", "-wd", str(tmp_path)])
    assert result.exit_code == 0
    assert "glance [drep]  [interrupted]" in result.stdout
    assert "did not finish; outputs may be partial" in result.stdout


def _derep_workdir(tmp_path: Path, reps: int) -> Path:
    """An ingested workdir with a finished dereplicate record over `reps` representatives."""
    from repgenr.core.contracts import write_clusters
    from repgenr.core.inputs import dir_stat_digest, manifest_digest_for_stage
    from repgenr.core.manifest import Manifest

    src = tmp_path / "src"
    src.mkdir()
    for i in range(4):
        (src / f"Fam_Gen_sp{i}_GCA_00000{i}.1.fasta").write_text(">s\nACGT\n", encoding="utf-8")
    wd = tmp_path / "wd"
    assert _runner.invoke(app, ["ingest", "-wd", str(wd), "--genomes-dir", str(src)]).exit_code == 0
    names = sorted(p.name for p in (wd / "genomes").iterdir())
    derep = wd / "derep"
    (derep / "representatives").mkdir(parents=True)
    write_clusters(derep / "clusters.tsv", {n: [] for n in names[:reps]})
    for n in names[:reps]:
        (derep / "representatives" / n).write_text(">s\nACGT\n", encoding="utf-8")
    for table in ("genome_status.tsv", "cluster_summary.tsv"):
        (derep / table).write_text("x\n", encoding="utf-8")
    manifest = Manifest.open_readonly(wd / "manifest.sqlite")
    try:
        manifest_digest = manifest_digest_for_stage("dereplicate", manifest)
    finally:
        manifest.close()
    cfg = Config.load(wd)
    cfg.record_stage(
        "dereplicate",
        tool="sourmash",
        completed="2026-01-01T00:00:00",
        inputs={"genomes": dir_stat_digest(wd / "genomes"), "manifest": manifest_digest},
    )
    cfg.save(wd)
    return wd


def test_status_marks_a_stage_whose_input_changed_as_stale(tmp_path: Path) -> None:
    # status used to show dereplicate as done after the genome set changed;
    # only doctor and the next run noticed.
    wd = _derep_workdir(tmp_path, reps=3)
    result = _runner.invoke(app, ["status", "-wd", str(wd)])
    assert "[done]    dereplicate [sourmash]" in result.stdout
    assert "Next: repgenr phylo" in result.stdout

    (tmp_path / "src" / "Fam_Gen_sp9_GCA_000009.1.fasta").write_text(">s\nA\n", encoding="utf-8")
    assert _runner.invoke(app, ["ingest", "-wd", str(wd), "--genomes-dir", str(tmp_path / "src")])
    result = _runner.invoke(app, ["status", "-wd", str(wd)])
    assert result.exit_code == 0
    assert "[stale]   dereplicate [sourmash]" in result.stdout
    assert "input changed: genomes" in result.stdout  # and the manifest
    assert "Next: repgenr dereplicate" in result.stdout
    doctor = _runner.invoke(app, ["doctor", "-wd", str(wd)])
    assert "[WARN] dereplicate: input(s) changed since completion (genomes" in doctor.stdout


def test_status_marks_a_stage_with_a_missing_output_as_stale(tmp_path: Path) -> None:
    wd = _derep_workdir(tmp_path, reps=3)
    (wd / "derep" / "genome_status.tsv").unlink()
    result = _runner.invoke(app, ["status", "-wd", str(wd)])
    assert "[stale]   dereplicate" in result.stdout
    assert "missing: derep/genome_status.tsv" in result.stdout


def test_status_warns_that_phylo_will_refuse_too_few_representatives(tmp_path: Path) -> None:
    wd = _derep_workdir(tmp_path, reps=2)
    result = _runner.invoke(app, ["status", "-wd", str(wd)])
    assert "Next: repgenr phylo" in result.stdout
    assert "derep/representatives holds 2 genome(s) and a tree needs at least 3" in result.stdout
    assert "--all-genomes" in result.stdout


def test_derep_stock_record_is_not_stale_after_a_new_dereplication(tmp_path: Path) -> None:
    wd = _derep_workdir(tmp_path, reps=3)
    cfg = Config.load(wd)
    cfg.record_stage(
        "derep_stock",
        params={"action": "pack", "name": "r1"},
        completed="2026-01-01T00:00:00",
        inputs={"derep/clusters.tsv": "an older run"},
    )
    cfg.save(wd)
    (wd / "derep" / "stock" / "r1").mkdir(parents=True)
    (wd / "derep" / "stock" / "r1" / "clusters.tsv").write_text("x\n", encoding="utf-8")
    result = _runner.invoke(app, ["status", "-wd", str(wd)])
    assert "[stale]" not in result.stdout
    doctor = _runner.invoke(app, ["doctor", "-wd", str(wd)])
    assert "derep_stock: input(s) changed" not in doctor.stdout
