"""dRep genome staging: hardlink plain FASTAs, decompress only .gz."""

from __future__ import annotations

import gzip
import os
from pathlib import Path

from repgenr.dereplicators.drep import _stage_genome


def test_stage_plain_fasta_is_hardlinked(tmp_path: Path) -> None:
    src = tmp_path / "g.fasta"
    src.write_text(">x\nACGT\n")
    dest = tmp_path / "staged"
    dest.mkdir()
    out = _stage_genome(src, dest)
    assert out.read_text() == ">x\nACGT\n"
    # hardlink shares the inode (no extra disk), not a separate copy
    assert os.path.samefile(out, src) or out.stat().st_ino == src.stat().st_ino


def test_stage_gz_is_decompressed(tmp_path: Path) -> None:
    src = tmp_path / "g.fasta.gz"
    with gzip.open(src, "wb") as fo:
        fo.write(b">x\nACGT\n")
    dest = tmp_path / "staged"
    dest.mkdir()
    out = _stage_genome(src, dest)
    assert out.name == "g.fasta"  # .gz stripped
    assert out.read_text() == ">x\nACGT\n"


def test_dereplicate_passes_genome_fofn(tmp_path: Path, monkeypatch) -> None:
    # dRep -g accepts a text file of genome paths; passing thousands of paths
    # on argv would hit the OS ARG_MAX limit at scale.
    import logging

    from repgenr.dereplicators import drep as drep_mod
    from repgenr.dereplicators.base import DerepParams
    from repgenr.dereplicators.drep import DrepDereplicator

    calls: list[list[str]] = []

    def fake_run_tool(caps, cmd, **kwargs):
        parts = [str(c) for c in cmd]
        calls.append(parts)
        wd = Path(parts[2])
        (wd / "dereplicated_genomes").mkdir(parents=True)
        (wd / "dereplicated_genomes" / "g1.fasta").write_text(">x\nA\n", encoding="utf-8")
        dt = wd / "data_tables"
        dt.mkdir()
        (dt / "Cdb.csv").write_text(
            "genome,secondary_cluster\ng1.fasta,1_1\ng2.fasta,1_1\n", encoding="utf-8"
        )
        return 0

    monkeypatch.setattr(drep_mod, "run_tool", fake_run_tool)
    genomes = []
    for name in ("g1.fasta", "g2.fasta"):
        p = tmp_path / name
        p.write_text(">x\nACGT\n", encoding="utf-8")
        genomes.append(p)

    params = DerepParams(primary_ani=0.9, secondary_ani=0.99, threads=1)
    result = DrepDereplicator().dereplicate(
        genomes, tmp_path / "out", params, logging.getLogger("test")
    )

    parts = calls[0]
    gidx = parts.index("-g")
    fofn = Path(parts[gidx + 1])
    assert parts[gidx + 2] == "--processors"  # exactly one value after -g
    assert fofn.name.endswith(".fofn")
    listed = fofn.read_text(encoding="utf-8").splitlines()
    assert len(listed) == 2
    assert [p.name for p in result.representatives] == ["g1.fasta"]
