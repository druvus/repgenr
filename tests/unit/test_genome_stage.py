"""End-to-end tests of the genome stage with a faked datasets CLI (WP2-4).

``_run_cmd`` is replaced by a fake that writes what NCBI datasets would:
a dehydrated zip, rehydrated .fna files, and the outgroup zip. The stage's
organization logic (naming, pruning, manifest write-back, FASTA sanity check)
runs for real.
"""

from __future__ import annotations

import zipfile
from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.errors import WorkdirError
from repgenr.core.manifest import GenomeRecord
from repgenr.stages import genome
from repgenr.stages.genome import GenomeParams

_SELECTED = [
    GenomeRecord(accession="GCF_000001.1", source="gtdb",
                 family="Francisellaceae", genus="Francisella", species="tularensis"),
    GenomeRecord(accession="GCF_000002.1", source="gtdb",
                 family="Francisellaceae", genus="Francisella", species="tularensis"),
]
_OUTGROUP = GenomeRecord(accession="GCF_000010.1", source="gtdb", is_outgroup=True,
                         family="Francisellaceae", genus="Francisella",
                         species="philomiragia")


@pytest.fixture()
def ctx(tmp_path, monkeypatch) -> WorkdirContext:
    monkeypatch.setattr(genome, "preflight", lambda caps: {"datasets": "16.0"})
    monkeypatch.setattr(genome, "_check_disk", lambda *a, **k: None)
    c = WorkdirContext(tmp_path / "wd", create=True)
    c.manifest.upsert_many([*_SELECTED, _OUTGROUP])
    return c


def _fake_run_cmd(monkeypatch, fasta: bytes = b">seq\nACGT\n") -> list[list[str]]:
    """datasets download -> zip with a marker; rehydrate -> per-accession .fna."""
    calls: list[list[str]] = []

    def fake(cmd, *, logger, **kw):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        if cmd[:3] == ["datasets", "download", "genome"] and "--dehydrated" in cmd:
            acc_file = Path(cmd[cmd.index("--inputfile") + 1])
            zip_path = Path(cmd[cmd.index("--filename") + 1])
            with zipfile.ZipFile(zip_path, "w") as zf:
                zf.writestr("README.md", acc_file.read_text(encoding="utf-8"))
        elif cmd[:2] == ["datasets", "rehydrate"]:
            extract = Path(cmd[cmd.index("--directory") + 1])
            readme = extract / "README.md"
            for acc in readme.read_text(encoding="utf-8").splitlines():
                if acc:
                    fna = extract / "data" / acc / f"{acc}_genomic.fna"
                    fna.parent.mkdir(parents=True, exist_ok=True)
                    fna.write_bytes(fasta)
        elif cmd[:3] == ["datasets", "download", "genome"]:  # outgroup, hydrated
            acc = cmd[4]
            zip_path = Path(cmd[cmd.index("--filename") + 1])
            with zipfile.ZipFile(zip_path, "w") as zf:
                zf.writestr(f"data/{acc}/{acc}_genomic.fna", fasta)
        else:
            raise AssertionError(f"unexpected datasets command: {cmd}")
        return 0

    monkeypatch.setattr(genome, "_run_cmd", fake)
    return calls


def test_run_downloads_and_organizes(ctx, monkeypatch) -> None:
    _fake_run_cmd(monkeypatch)
    count = genome.run(ctx, GenomeParams())
    assert count == 2

    fastas = sorted(p.name for p in ctx.genomes_dir.iterdir())
    assert fastas == [
        "Francisellaceae_Francisella_tularensis_GCF_000001.1.fasta",
        "Francisellaceae_Francisella_tularensis_GCF_000002.1.fasta",
    ]
    out = list(ctx.outgroup_dir.iterdir())
    assert [p.name for p in out] == [
        "Francisellaceae_Francisella_philomiragia_GCF_000010.1.fasta"
    ]
    # filenames recorded back into the manifest
    by_acc = {g.accession: g for g in ctx.manifest.all_genomes(include_outgroup=False)}
    assert by_acc["GCF_000001.1"].filename == fastas[0]
    assert "genome" in ctx.config.stages


def test_run_skips_present_and_prunes_stale(ctx, monkeypatch) -> None:
    calls = _fake_run_cmd(monkeypatch)
    ctx.genomes_dir.mkdir(parents=True, exist_ok=True)
    present = ctx.genomes_dir / "Francisellaceae_Francisella_tularensis_GCF_000001.1.fasta"
    present.write_text(">seq\nACGT\n", encoding="utf-8")
    stale = ctx.genomes_dir / "Old_Stale_genome_GCF_999999.9.fasta"
    stale.write_text(">seq\nACGT\n", encoding="utf-8")
    keepme = ctx.genomes_dir / "notes.txt"
    keepme.write_text("user file", encoding="utf-8")

    genome.run(ctx, GenomeParams())
    assert not stale.exists(), "stale genome FASTA is pruned"
    assert keepme.exists(), "non-FASTA files are never pruned"
    downloaded = [c for c in calls if "--dehydrated" in c]
    acc_list = (ctx.workdir / "ncbi_acc_download_list.txt").read_text(encoding="utf-8")
    assert "GCF_000001.1" not in acc_list, "already-present genome is not re-downloaded"
    assert downloaded, "the missing genome still downloads"


def test_run_accession_list_only_stops_before_download(ctx, monkeypatch) -> None:
    calls = _fake_run_cmd(monkeypatch)
    assert genome.run(ctx, GenomeParams(accession_list_only=True)) == 0
    assert (ctx.workdir / "ncbi_acc_download_list.txt").exists()
    assert calls == []


def test_run_empty_manifest_raises(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(genome, "preflight", lambda caps: {})
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(WorkdirError, match="metadata stage"):
        genome.run(ctx, GenomeParams())


def test_non_fasta_download_raises(ctx, monkeypatch) -> None:
    _fake_run_cmd(monkeypatch, fasta=b"<html>Service unavailable</html>")
    with pytest.raises(WorkdirError, match="not FASTA"):
        genome.run(ctx, GenomeParams())


def test_missing_accessions_logged_not_fatal(ctx, monkeypatch, caplog) -> None:
    """NCBI returning no genome for an accession warns and continues."""
    import logging

    def fake(cmd, *, logger, **kw):
        cmd = [str(c) for c in cmd]
        if "--dehydrated" in cmd:
            zip_path = Path(cmd[cmd.index("--filename") + 1])
            with zipfile.ZipFile(zip_path, "w") as zf:
                zf.writestr("README.md", "")
        elif cmd[:2] == ["datasets", "rehydrate"]:
            extract = Path(cmd[cmd.index("--directory") + 1])
            # only the first accession comes back
            fna = extract / "data" / "GCF_000001.1" / "x.fna"
            fna.parent.mkdir(parents=True, exist_ok=True)
            fna.write_bytes(b">seq\nACGT\n")
        elif cmd[:3] == ["datasets", "download", "genome"]:
            zip_path = Path(cmd[cmd.index("--filename") + 1])
            with zipfile.ZipFile(zip_path, "w") as zf:
                zf.writestr("data/GCF_000010.1/x.fna", b">seq\nACGT\n")
        return 0

    monkeypatch.setattr(genome, "_run_cmd", fake)
    with caplog.at_level(logging.WARNING):
        count = genome.run(ctx, GenomeParams())
    assert count == 2
    assert any("no genome" in r.getMessage() for r in caplog.records)
