"""`repgenr doctor`: read-only workdir health checks (outputs vs records)."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core.config import Config
from repgenr.core.contracts import SelectionRow, write_clusters, write_selection
from repgenr.core.doctor import diagnose
from repgenr.core.manifest import GenomeRecord, Manifest

_runner = CliRunner()


def _levels(findings, area: str | None = None) -> set[str]:
    return {f.level for f in findings if area is None or f.area == area}


def _messages(findings, level: str | None = None) -> str:
    return "\n".join(f.message for f in findings if level is None or f.level == level)


def _base_workdir(tmp_path: Path) -> Path:
    """A healthy two-genome workdir: selection + manifest + genomes + outgroup."""
    wd = tmp_path / "wd"
    wd.mkdir()
    rows = [
        SelectionRow("GCF_1.1", "Fam", "Gen", "sp1", False, "Fam_Gen_sp1_GCF_1.1.fasta"),
        SelectionRow("GCF_2.1", "Fam", "Gen", "sp2", False, "Fam_Gen_sp2_GCF_2.1.fasta"),
        SelectionRow("GCF_9.1", "Fam", "Out", "grp", True, "Fam_Out_grp_GCF_9.1.fasta"),
    ]
    write_selection(wd / "selection.tsv", rows)
    (wd / "genomes").mkdir()
    for r in rows[:2]:
        (wd / "genomes" / r.filename).write_text(">x\nACGT\n", encoding="utf-8")
    (wd / "outgroup").mkdir()
    (wd / "outgroup" / rows[2].filename).write_text(">og\nACGT\n", encoding="utf-8")
    (wd / "outgroup_accession.txt").write_text("GCF_9.1\n", encoding="utf-8")
    manifest = Manifest(wd / "manifest.sqlite")
    manifest.replace_genomes([
        GenomeRecord(accession="GCF_1.1", filename=rows[0].filename),
        GenomeRecord(accession="GCF_2.1", filename=rows[1].filename),
        GenomeRecord(accession="GCF_9.1", filename=rows[2].filename, is_outgroup=True),
    ])
    manifest.close()
    cfg = Config()
    cfg.record_stage("metadata", completed="2026-01-01T00:00:00")
    cfg.record_stage("genome", completed="2026-01-01T00:01:00")
    cfg.save(wd)
    return wd


def test_healthy_workdir_has_no_failures(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    findings = diagnose(wd)
    assert "fail" not in _levels(findings), _messages(findings, "fail")


def test_no_workdir_reports_cleanly(tmp_path: Path) -> None:
    findings = diagnose(tmp_path / "nope")
    assert findings and findings[0].level == "warn"


def test_interrupted_stage_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    cfg = Config.load(wd)
    cfg.record_stage("dereplicate", tool="skder", params={"x": 1})  # no completed
    cfg.save(wd)
    findings = diagnose(wd)
    assert any(f.level == "fail" and "dereplicate" in f.area for f in findings)
    assert "re-run" in _messages(findings, "fail")


def test_missing_genome_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "genomes" / "Fam_Gen_sp2_GCF_2.1.fasta").unlink()
    findings = diagnose(wd)
    assert any(
        f.level == "fail" and "Fam_Gen_sp2_GCF_2.1.fasta" in f.message for f in findings
    )


def test_excused_missing_accession_is_not_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "genomes" / "Fam_Gen_sp2_GCF_2.1.fasta").unlink()
    (wd / "missing_accessions.txt").write_text("GCF_2.1\n", encoding="utf-8")
    findings = diagnose(wd)
    assert "fail" not in _levels(findings), _messages(findings, "fail")


def test_non_fasta_genome_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "genomes" / "Fam_Gen_sp1_GCF_1.1.fasta").write_text(
        "<html>error</html>", encoding="utf-8"
    )
    findings = diagnose(wd)
    assert any(
        f.level == "fail" and "Fam_Gen_sp1_GCF_1.1.fasta" in f.message for f in findings
    )


def test_manifest_selection_drift_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    manifest = Manifest(wd / "manifest.sqlite")
    manifest.upsert_many([GenomeRecord(accession="GCF_STALE.1")])
    manifest.close()
    findings = diagnose(wd)
    assert any(f.level == "fail" and "GCF_STALE.1" in f.message for f in findings)


def test_representatives_clusters_mismatch_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    reps = wd / "derep" / "representatives"
    reps.mkdir(parents=True)
    (reps / "Fam_Gen_sp1_GCF_1.1.fasta").write_text(">x\nACGT\n", encoding="utf-8")
    write_clusters(wd / "derep" / "clusters.tsv", {
        "Fam_Gen_sp1_GCF_1.1.fasta": [],
        "Fam_Gen_sp2_GCF_2.1.fasta": [],  # listed but absent on disk
    })
    findings = diagnose(wd)
    assert any(
        f.level == "fail" and "Fam_Gen_sp2_GCF_2.1.fasta" in f.message for f in findings
    )


def test_truncated_tree_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "tree").mkdir()
    (wd / "tree" / "tree.nwk").write_text("(a,b", encoding="utf-8")  # no ';'
    findings = diagnose(wd)
    assert any(f.level == "fail" and "tree.nwk" in f.message for f in findings)


def test_mismatched_tree2tax_pair_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "tree2tax.tsv").write_text("child\tparent\n", encoding="utf-8")
    # genomes_map.tsv missing
    findings = diagnose(wd)
    assert any(f.level == "fail" and "genomes_map.tsv" in f.message for f in findings)


def test_unresolvable_outgroup_is_a_warning(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "outgroup" / "Fam_Out_grp_GCF_9.1.fasta").unlink()
    findings = diagnose(wd)
    assert any(f.level == "warn" and "GCF_9.1" in f.message for f in findings)


def test_leftover_temp_files_are_a_warning(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "tree").mkdir()
    (wd / "tree" / "tree.nwk.part").write_text("(a", encoding="utf-8")
    findings = diagnose(wd)
    assert any(f.level == "warn" and "tree.nwk.part" in f.message for f in findings)


def test_changed_inputs_are_a_warning(tmp_path: Path) -> None:
    """A completed stage whose recorded input digests no longer match reality
    is stale (it will re-run) -- doctor should say so."""
    wd = _base_workdir(tmp_path)
    cfg = Config.load(wd)
    cfg.record_stage(
        "dereplicate", tool="skder", params={},
        completed="2026-01-01T00:02:00", fingerprint="f",
        inputs={"genomes": "old-digest"},
    )
    cfg.save(wd)
    findings = diagnose(wd)
    assert any(
        f.level == "warn" and f.area == "dereplicate" and "changed" in f.message
        for f in findings
    )


def test_diagnose_is_read_only(tmp_path: Path) -> None:
    """Doctor must not create files in a workdir (notably no manifest.sqlite)."""
    wd = tmp_path / "wd"
    wd.mkdir()
    Config().save(wd)
    before = sorted(p.name for p in wd.iterdir())
    diagnose(wd)
    assert sorted(p.name for p in wd.iterdir()) == before


# --- CLI ----------------------------------------------------------------------


def test_cli_doctor_exit_codes(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    result = _runner.invoke(app, ["doctor", "-wd", str(wd)])
    assert result.exit_code == 0
    assert "OK" in result.stdout or "ok" in result.stdout

    (wd / "genomes" / "Fam_Gen_sp2_GCF_2.1.fasta").unlink()
    result = _runner.invoke(app, ["doctor", "-wd", str(wd)])
    assert result.exit_code == 1
    assert "FAIL" in result.stdout
