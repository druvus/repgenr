"""Genome sketches with a real sourmash: written by ingest, completed by
``repgenr sketch``, and copied by ``ingest --from-workdir``."""

from __future__ import annotations

import csv
import io
import zipfile
from pathlib import Path

import pytest

from repgenr.core.contracts import list_fasta, record_name

pytestmark = [pytest.mark.live, pytest.mark.requires_binary("sourmash")]


def _signatures(path: Path) -> list[tuple[str, int, str]]:
    """(name, ksize, scaled) of each signature in a sourmash zip collection."""
    with zipfile.ZipFile(path) as zf:
        text = zf.read("SOURMASH-MANIFEST.csv").decode("utf-8")
    rows = csv.DictReader(io.StringIO("\n".join(text.splitlines()[1:])))
    return sorted((r["name"], int(r["ksize"]), r["scaled"]) for r in rows)


def test_ingest_sketches_and_sketch_command(synthetic_set, ingested_workdir, run_repgenr) -> None:
    genomes = synthetic_set("balanced", n=4, length=30_000)
    wd = ingested_workdir(genomes)
    sketches = wd / "sketches"
    names = {record_name(g) for g in list_fasta(wd / "genomes")}
    assert {p.name for p in sketches.glob("*.sig.zip")} == {f"{n}.sig.zip" for n in names}
    for name in names:
        assert _signatures(sketches / f"{name}.sig.zip") == [
            (name, 21, "1000"),
            (name, 31, "1000"),
            (name, 51, "1000"),
        ]

    # Every sketch is current: the command writes none.
    out = run_repgenr("sketch", "-wd", wd, "-t", "2")
    assert "4 present, 0 written, 0 stale replaced, 0 removed" in out.stdout + out.stderr

    # A sketch removed by hand is written again, and only that one.
    victim = sorted(sketches.glob("*.sig.zip"))[0]
    victim.unlink()
    out = run_repgenr("sketch", "-wd", wd, "-t", "2")
    assert "3 present, 1 written" in out.stdout + out.stderr
    assert victim.is_file()
    assert "sketches: 4/4" in run_repgenr("status", "-wd", wd).stdout


def test_no_sketch_and_from_workdir_copy(synthetic_set, ingested_workdir, run_repgenr, tmp_path):
    genomes = synthetic_set("balanced", n=3, length=30_000)
    first = ingested_workdir(genomes, name="first")
    assert len(list((first / "sketches").glob("*.sig.zip"))) == 3

    bare = tmp_path / "bare"
    run_repgenr("ingest", "-wd", bare, "--genomes-dir", genomes, "--no-sketch")
    assert not (bare / "sketches").exists()

    second = tmp_path / "second"
    out = run_repgenr("ingest", "-wd", second, "--from-workdir", first)
    assert "3 copied" in out.stdout + out.stderr
    assert len(list((second / "sketches").glob("*.sig.zip"))) == 3
