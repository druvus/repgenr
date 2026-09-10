"""Auxiliary commands and global behaviour on an ingested, dereplicated workdir:
derep-unpack, derep-stock, status, versions, doctor, --force, logging flags
and the REPGENR_* environment variables, all as real processes."""

from __future__ import annotations

from pathlib import Path

import pytest
from live_helpers import representatives

from repgenr.core.config import Config

pytestmark = [pytest.mark.live, pytest.mark.requires_binary("sourmash")]


@pytest.fixture
def derep_wd(synthetic_set, ingested_workdir, run_repgenr) -> Path:
    genomes = synthetic_set("clonal", n=8, length=50_000)
    wd = ingested_workdir(genomes, copy=True)
    run_repgenr("dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2")
    return wd


def test_derep_unpack_with_and_without_representant(run_repgenr, derep_wd: Path) -> None:
    run_repgenr("derep-unpack", "-wd", derep_wd)
    unpacked = derep_wd / "derep" / "unpacked"
    reps = representatives(derep_wd)
    assert {p.name for p in unpacked.iterdir()} == {Path(r).stem for r in reps}
    with_rep = {p.name for d in unpacked.iterdir() for p in d.iterdir()}
    assert reps <= with_rep

    run_repgenr("--force", "derep-unpack", "-wd", derep_wd, "--no-representant")
    without = {p.name for d in unpacked.iterdir() for p in d.iterdir()}
    assert not (reps & without), "representatives are left out of their cluster directories"


def test_derep_stock_round_trip(run_repgenr, derep_wd: Path) -> None:
    before = representatives(derep_wd)
    run_repgenr("derep-stock", "-wd", derep_wd, "--action", "pack", "--name", "run1")
    listed = run_repgenr("derep-stock", "-wd", derep_wd, "--action", "list")
    assert "run1" in listed.stdout + listed.stderr
    # A stricter run changes the representatives; unpack restores the stored set.
    run_repgenr("dereplicate", "-wd", derep_wd, "--tool", "sourmash", "-t", "2", "-sani", "0.999")
    assert representatives(derep_wd) != before
    run_repgenr("derep-stock", "-wd", derep_wd, "--action", "unpack", "--name", "run1")
    assert representatives(derep_wd) == before
    run_repgenr("derep-stock", "-wd", derep_wd, "--action", "delete", "--name", "run1")
    assert not (derep_wd / "derep" / "stock" / "run1").exists()
    missing = run_repgenr(
        "derep-stock", "-wd", derep_wd, "--action", "unpack", "--name", "run1", check=False
    )
    assert missing.returncode != 0


def test_status_and_versions(run_repgenr, derep_wd: Path, tmp_path: Path) -> None:
    status = run_repgenr("status", "-wd", derep_wd).stdout
    assert "Pipeline: local" in status
    assert "[done]    ingest" in status and "[done]    dereplicate [sourmash]" in status
    assert "[next] phylo" in status
    printed = run_repgenr("versions", "-wd", derep_wd).stdout
    assert printed.startswith("sourmash:") or "\nsourmash:" in printed
    out = tmp_path / "versions.yml"
    run_repgenr("versions", "-wd", derep_wd, "--versions-out", out)
    assert "    sourmash:" in out.read_text(encoding="utf-8")


def test_doctor_passes_then_fails_on_a_corrupt_genome(run_repgenr, derep_wd: Path) -> None:
    ok = run_repgenr("doctor", "-wd", derep_wd)
    assert "0 failure(s)" in ok.stdout
    victim = next(iter(sorted((derep_wd / "genomes").glob("*.fasta"))))
    victim.write_text("this is not a fasta file\n", encoding="utf-8")
    bad = run_repgenr("doctor", "-wd", derep_wd, check=False)
    assert bad.returncode == 1 and "genomes" in bad.stdout


def test_second_run_skips_and_force_reruns(run_repgenr, derep_wd: Path) -> None:
    stamp = Config.load(derep_wd).stages["dereplicate"].completed
    again = run_repgenr("dereplicate", "-wd", derep_wd, "--tool", "sourmash", "-t", "2")
    assert "skipping (use --force to re-run)" in again.stderr + again.stdout
    assert Config.load(derep_wd).stages["dereplicate"].completed == stamp

    forced = run_repgenr("--force", "dereplicate", "-wd", derep_wd, "--tool", "sourmash", "-t", "2")
    assert "skipping" not in forced.stderr + forced.stdout
    assert Config.load(derep_wd).stages["dereplicate"].completed != stamp

    env_stamp = Config.load(derep_wd).stages["dereplicate"].completed
    run_repgenr(
        "dereplicate", "-wd", derep_wd, "--tool", "sourmash", "-t", "2", env={"REPGENR_FORCE": "1"}
    )
    assert Config.load(derep_wd).stages["dereplicate"].completed != env_stamp


def test_logging_flags_and_env(run_repgenr, derep_wd: Path) -> None:
    quiet = run_repgenr("--quiet", "--force", "dereplicate", "-wd", derep_wd, "--tool", "sourmash")
    assert " INFO " not in quiet.stderr + quiet.stdout, "quiet suppresses INFO on the console"
    log = derep_wd / "repgenr.log"
    size_before = log.stat().st_size
    verbose = run_repgenr(
        "--verbose", "--force", "dereplicate", "-wd", derep_wd, "--tool", "sourmash"
    )
    assert "DEBUG" in verbose.stderr + verbose.stdout
    assert "DEBUG" in log.read_text(encoding="utf-8")[size_before:]
    size_before = log.stat().st_size
    run_repgenr(
        "--force",
        "dereplicate",
        "-wd",
        derep_wd,
        "--tool",
        "sourmash",
        env={"REPGENR_LOG_LEVEL": "DEBUG"},
    )
    assert "DEBUG" in log.read_text(encoding="utf-8")[size_before:]
