"""Manifest quality as an adapter input: galah's input order, dRep's --genomeInfo.

galah prefers earlier genomes as representatives when it has no quality, and
dRep scores genomes with CheckM unless it is given their quality. Both adapters
read ``DerepParams.quality`` (filename -> completeness, contamination), which
the stage fills from the manifest and the chunk/merge steps from selection.tsv.
"""

from __future__ import annotations

import csv
import gzip
import logging
import shutil
from pathlib import Path

import pytest

import repgenr.dereplicators.drep as drep_mod
import repgenr.dereplicators.galah as galah_mod
from repgenr.dereplicators.base import DerepParams

_LOG = logging.getLogger("test")


def _write(path: Path, length: int) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">{path.stem}\n{'ACGT' * (length // 4)}\n", encoding="utf-8")
    return path


# --- galah -------------------------------------------------------------------


def _galah_first_listed_is_rep(seen: list[list[str]]):
    """Fake galah that, like galah without quality, keeps the first listed genome."""

    def run_tool(caps, command, **kwargs):
        cmd = [str(c) for c in command]
        listed = (
            Path(cmd[cmd.index("--genome-fasta-list") + 1]).read_text(encoding="utf-8").splitlines()
        )
        seen.append([Path(p).name for p in listed])
        rep = listed[0]
        rows = [f"{rep}\t{p}" for p in listed]
        Path(cmd[cmd.index("--output-cluster-definition") + 1]).write_text(
            "\n".join(rows) + "\n", encoding="utf-8"
        )
        return 0

    return run_tool


def test_galah_without_quality_lists_the_largest_genome_first(tmp_path, monkeypatch) -> None:
    # A 40 percent fragment whose name sorts first must not be listed first.
    fragment = _write(tmp_path / "g" / "a_fragment.fasta", 400)
    complete = _write(tmp_path / "g" / "b_complete.fasta", 1000)
    middle = _write(tmp_path / "g" / "c_middle.fasta", 800)
    seen: list[list[str]] = []
    monkeypatch.setattr(galah_mod, "run_tool", _galah_first_listed_is_rep(seen))

    result = galah_mod.GalahDereplicator().dereplicate(
        [fragment, complete, middle], tmp_path / "out", DerepParams(), _LOG
    )
    assert seen == [["b_complete.fasta", "c_middle.fasta", "a_fragment.fasta"]]
    assert [p.name for p in result.representatives] == ["b_complete.fasta"]


def test_galah_with_quality_for_every_genome_keeps_the_given_order(tmp_path, monkeypatch) -> None:
    # With quality, --keeper quality re-picks each cluster's representative.
    fragment = _write(tmp_path / "g" / "a_fragment.fasta", 400)
    complete = _write(tmp_path / "g" / "b_complete.fasta", 1000)
    seen: list[list[str]] = []
    monkeypatch.setattr(galah_mod, "run_tool", _galah_first_listed_is_rep(seen))
    quality = {"a_fragment.fasta": (40.0, 0.0), "b_complete.fasta": (99.0, 0.1)}

    galah_mod.GalahDereplicator().dereplicate(
        [fragment, complete], tmp_path / "out", DerepParams(quality=quality), _LOG
    )
    assert seen == [["a_fragment.fasta", "b_complete.fasta"]]


def test_galah_with_partial_quality_orders_by_size(tmp_path, monkeypatch) -> None:
    fragment = _write(tmp_path / "g" / "a_fragment.fasta", 400)
    complete = _write(tmp_path / "g" / "b_complete.fasta", 1000)
    seen: list[list[str]] = []
    monkeypatch.setattr(galah_mod, "run_tool", _galah_first_listed_is_rep(seen))

    galah_mod.GalahDereplicator().dereplicate(
        [fragment, complete],
        tmp_path / "out",
        DerepParams(quality={"b_complete.fasta": (99.0, 0.1)}),
        _LOG,
    )
    assert seen == [["b_complete.fasta", "a_fragment.fasta"]]


# --- dRep --------------------------------------------------------------------


def _fake_drep(seen: list[dict]):
    """Record argv and the --genomeInfo table; write a one-cluster result."""

    def run_tool(caps, command, **kwargs):
        cmd = [str(c) for c in command]
        info = None
        if "--genomeInfo" in cmd:
            with open(cmd[cmd.index("--genomeInfo") + 1], encoding="utf-8", newline="") as fo:
                info = list(csv.reader(fo))
        seen.append({"cmd": cmd, "info": info})
        wd = Path(cmd[2])
        staged = [
            Path(p) for p in Path(cmd[cmd.index("-g") + 1]).read_text(encoding="utf-8").splitlines()
        ]
        (wd / "dereplicated_genomes").mkdir(parents=True)
        shutil.copy2(staged[0], wd / "dereplicated_genomes" / staged[0].name)
        (wd / "data_tables").mkdir()
        rows = ["genome,secondary_cluster", *(f"{p.name},1_1" for p in staged)]
        (wd / "data_tables" / "Cdb.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")
        return 0

    return run_tool


@pytest.fixture()
def drep_genomes(tmp_path) -> list[Path]:
    out = [_write(tmp_path / "g" / f"g{i}.fasta", 100) for i in (1, 2, 3)]
    gz = tmp_path / "g" / "g4.fasta.gz"
    gz.write_bytes(gzip.compress(b">g4\nACGTACGT\n"))
    return [*out, gz]


_QUALITY = {
    "g1.fasta": (99.5, 0.2),
    "g2.fasta": (97.0, 1.5),
    "g3.fasta": (88.25, 3.0),
    "g4.fasta.gz": (91.0, 0.0),
}


def test_drep_receives_manifest_quality_as_genome_info(drep_genomes, tmp_path, monkeypatch) -> None:
    seen: list[dict] = []
    monkeypatch.setattr(drep_mod, "run_tool", _fake_drep(seen))
    drep_mod.DrepDereplicator().dereplicate(
        drep_genomes, tmp_path / "out", DerepParams(quality=_QUALITY), _LOG
    )
    (call,) = seen
    assert "--ignoreGenomeQuality" not in call["cmd"]
    # dRep matches rows by the basename of the file it reads: the gzipped
    # genome is staged decompressed, so its row names the staged copy.
    assert call["info"] == [
        ["genome", "completeness", "contamination"],
        ["g1.fasta", "99.5", "0.2"],
        ["g2.fasta", "97", "1.5"],
        ["g3.fasta", "88.25", "3"],
        ["g4.fasta", "91", "0"],
    ]


def test_drep_without_quality_for_every_genome_gets_no_genome_info(
    drep_genomes, tmp_path, monkeypatch, caplog
) -> None:
    seen: list[dict] = []
    monkeypatch.setattr(drep_mod, "run_tool", _fake_drep(seen))
    partial = {k: v for k, v in _QUALITY.items() if k != "g3.fasta"}
    with caplog.at_level(logging.INFO, logger="test"):
        drep_mod.DrepDereplicator().dereplicate(
            drep_genomes, tmp_path / "out", DerepParams(quality=partial), _LOG
        )
    (call,) = seen
    assert "--genomeInfo" not in call["cmd"]
    assert any("1 of 4" in r.getMessage() and "g3.fasta" in r.getMessage() for r in caplog.records)


def test_drep_without_any_quality_gets_no_genome_info(drep_genomes, tmp_path, monkeypatch) -> None:
    seen: list[dict] = []
    monkeypatch.setattr(drep_mod, "run_tool", _fake_drep(seen))
    drep_mod.DrepDereplicator().dereplicate(drep_genomes, tmp_path / "out", DerepParams(), _LOG)
    assert "--genomeInfo" not in seen[0]["cmd"]


def test_drep_virus_mode_ignores_quality(drep_genomes, tmp_path, monkeypatch) -> None:
    seen: list[dict] = []
    monkeypatch.setattr(drep_mod, "run_tool", _fake_drep(seen))
    params = DerepParams(quality=_QUALITY, extra={"virus": True})
    drep_mod.DrepDereplicator().dereplicate(drep_genomes, tmp_path / "out", params, _LOG)
    (call,) = seen
    assert "--ignoreGenomeQuality" in call["cmd"]
    assert "--genomeInfo" not in call["cmd"]
