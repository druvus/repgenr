"""`repgenr doctor`: read-only workdir health checks (outputs vs records)."""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest
from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core.config import Config
from repgenr.core.contracts import SelectionRow, write_clusters, write_selection
from repgenr.core.doctor import diagnose
from repgenr.core.errors import WorkdirError
from repgenr.core.inputs import dir_stat_digest, manifest_digest_for_stage
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


def test_dereplicate_completion_with_derep_status_is_not_stale(tmp_path: Path) -> None:
    """doctor's stale-input check must use the same per-stage manifest digest
    dereplicate's own resume fingerprint stamps (include_derep=False, since
    dereplicate WRITES derep_status/representative itself) -- otherwise every
    dereplicated workdir falsely reports "input(s) changed" for dereplicate on
    every run, forever, since the two digests can never agree."""
    wd = _base_workdir(tmp_path)
    manifest = Manifest(wd / "manifest.sqlite")
    manifest.set_derep_status_many([
        ("GCF_1.1", "representative", None),
        ("GCF_2.1", "contained", "Fam_Gen_sp1_GCF_1.1.fasta"),
    ])
    manifest.close()

    genomes_digest = dir_stat_digest(wd / "genomes")
    ro = Manifest.open_readonly(wd / "manifest.sqlite")
    try:
        manifest_digest_value = manifest_digest_for_stage("dereplicate", ro)
    finally:
        ro.close()

    cfg = Config.load(wd)
    cfg.record_stage(
        "dereplicate", tool="skder", params={},
        completed="2026-01-01T00:02:00", fingerprint="f",
        inputs={"genomes": genomes_digest, "manifest": manifest_digest_value},
    )
    cfg.save(wd)

    findings = diagnose(wd)
    assert not any(f.area == "dereplicate" and f.level == "warn" for f in findings)

    # A real quality-only manifest edit must still be caught as stale.
    manifest2 = Manifest(wd / "manifest.sqlite")
    manifest2.upsert_many([
        GenomeRecord(
            accession="GCF_1.1", filename="Fam_Gen_sp1_GCF_1.1.fasta",
            completeness=98.0, contamination=0.5,
        )
    ])
    manifest2.close()

    findings2 = diagnose(wd)
    dereplicate_warnings = [
        f for f in findings2 if f.area == "dereplicate" and f.level == "warn"
    ]
    assert len(dereplicate_warnings) == 1
    assert "manifest" in dereplicate_warnings[0].message


def test_diagnose_is_read_only(tmp_path: Path) -> None:
    """Doctor must not create files in a workdir (notably no manifest.sqlite)."""
    wd = tmp_path / "wd"
    wd.mkdir()
    Config().save(wd)
    before = sorted(p.name for p in wd.iterdir())
    diagnose(wd)
    assert sorted(p.name for p in wd.iterdir()) == before


# --- CLI ----------------------------------------------------------------------


def _non_wal_listing(workdir: Path) -> list[str]:
    # A WAL-mode database opened read-only may cause sqlite to create a
    # "-shm" (and sometimes "-wal") sidecar purely to read the wal-index;
    # a read-only connection cannot check those back in on close. They are
    # not a write to the manifest itself, so they are excluded here -- the
    # mtime assertion below is what actually proves the manifest untouched.
    return sorted(
        p.name for p in workdir.iterdir()
        if not (p.name.endswith("-shm") or p.name.endswith("-wal"))
    )


def test_doctor_leaves_manifest_untouched(workdir: Path) -> None:
    workdir.mkdir(parents=True)
    with Manifest.open(workdir) as m:
        m.upsert(GenomeRecord(accession="GCA_1", filename="x.fasta"))
    (workdir / "selection.tsv").write_text("accession\tfilename\nGCA_1\tx.fasta\n")
    Config().save(workdir)  # so diagnose() reaches the manifest-drift check
    manifest = workdir / "manifest.sqlite"
    before = (manifest.stat().st_mtime_ns, _non_wal_listing(workdir))

    diagnose(workdir)

    after = (manifest.stat().st_mtime_ns, _non_wal_listing(workdir))
    assert before == after


def test_open_readonly_refuses_writes(workdir: Path) -> None:
    import sqlite3

    workdir.mkdir(parents=True)
    with Manifest.open(workdir) as m:
        m.upsert(GenomeRecord(accession="GCA_1"))
    ro = Manifest.open_readonly(workdir / "manifest.sqlite")
    try:
        assert [g.accession for g in ro.all_genomes()] == ["GCA_1"]
        with pytest.raises(sqlite3.OperationalError):
            ro.upsert(GenomeRecord(accession="GCA_2"))
    finally:
        ro.close()


def test_open_readonly_missing_file_raises_workdir_error(tmp_path: Path) -> None:
    with pytest.raises(WorkdirError):
        Manifest.open_readonly(tmp_path / "absent.sqlite")


def test_open_readonly_unwritable_dir_raises_workdir_error(workdir: Path) -> None:
    """A WAL manifest needs a writable directory for its -shm file even to be
    opened read-only; that failure must surface as WorkdirError, not a raw
    sqlite3.OperationalError, so doctor's callers get the documented contract."""
    if os.geteuid() == 0:
        pytest.skip("root bypasses directory permissions")
    workdir.mkdir(parents=True)
    with Manifest.open(workdir) as m:
        m.upsert(GenomeRecord(accession="GCA_1"))
    manifest_path = workdir / "manifest.sqlite"
    original_mode = stat.S_IMODE(workdir.stat().st_mode)
    workdir.chmod(0o555)
    try:
        with pytest.raises(WorkdirError):
            Manifest.open_readonly(manifest_path)
    finally:
        workdir.chmod(original_mode)


def test_open_readonly_unreadable_file_raises_workdir_error(workdir: Path) -> None:
    """Even connect() itself can fail (manifest file unreadable, e.g. chmod
    0o000); that must also surface as WorkdirError, not a raw sqlite error."""
    if os.geteuid() == 0:
        pytest.skip("root bypasses file permissions")
    workdir.mkdir(parents=True)
    with Manifest.open(workdir) as m:
        m.upsert(GenomeRecord(accession="GCA_1"))
    manifest_path = workdir / "manifest.sqlite"
    original_mode = stat.S_IMODE(manifest_path.stat().st_mode)
    manifest_path.chmod(0o000)
    try:
        with pytest.raises(WorkdirError):
            Manifest.open_readonly(manifest_path)
    finally:
        manifest_path.chmod(original_mode)


def test_cli_doctor_exit_codes(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    result = _runner.invoke(app, ["doctor", "-wd", str(wd)])
    assert result.exit_code == 0
    assert "OK" in result.stdout or "ok" in result.stdout

    (wd / "genomes" / "Fam_Gen_sp2_GCF_2.1.fasta").unlink()
    result = _runner.invoke(app, ["doctor", "-wd", str(wd)])
    assert result.exit_code == 1
    assert "FAIL" in result.stdout
