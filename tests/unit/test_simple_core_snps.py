"""Tests for the simple SNP typer's core-SNP reduction."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.snptypers.simple import _write_core_snps


def test_core_snp_reduction(tmp_path: Path) -> None:
    consensuses = {
        "ref": "ACGTACGT",
        "s1": "ACGAACGT",  # differs at col 3
        "s2": "ACGTACGA",  # differs at col 7
    }
    core = tmp_path / "core.fasta"
    matrix = tmp_path / "dist.tsv"
    n = _write_core_snps(consensuses, core, matrix)
    assert n == 2  # columns 3 and 7 are variable

    records = _read_fasta(core)
    assert records["ref"] == "TT"  # ref bases at the two variable columns
    assert records["s1"] == "AT"
    assert records["s2"] == "TA"

    # distance matrix: ref vs s1 = 1, ref vs s2 = 1, s1 vs s2 = 2
    lines = matrix.read_text().splitlines()
    assert lines[0].split("\t")[1:] == ["ref", "s1", "s2"]


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:]
            records[name] = ""
        elif name is not None:
            records[name] += line.strip()
    return records


def test_thread_split_prefers_concurrent_genomes() -> None:
    """The budget buys workers first; threads only once genomes run out."""
    from repgenr.snptypers.simple import _split_threads

    assert _split_threads(8, 68) == (8, 1)
    assert _split_threads(8, 2) == (2, 4)
    assert _split_threads(1, 68) == (1, 1)
    assert _split_threads(8, 0) == (1, 1)
    assert _split_threads(0, 4) == (1, 1)


def test_per_genome_chain_is_threaded_and_clears_its_scratch(tmp_path: Path, monkeypatch) -> None:
    """Every tool gets the thread count, and the intermediates go once read."""
    from repgenr.snptypers import simple as mod

    calls: list[list[str]] = []

    def fake_run_tool(caps, cmd, **kw):  # noqa: ANN001
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        tool = cmd[0] if cmd[0] != "bcftools" else f"bcftools {cmd[1]}"
        out = kw.get("stdout_path")
        if out is None and "-o" in cmd:
            out = cmd[cmd.index("-o") + 1]
        if tool == "bcftools consensus":
            Path(out).write_text(">ref\nACGT\n", encoding="utf-8")
        elif tool == "bcftools index":
            Path(cmd[-1] + ".csi").write_text("", encoding="utf-8")
        elif tool == "samtools index":
            Path(cmd[-1] + ".bai").write_text("", encoding="utf-8")
        elif out is not None:
            Path(out).write_text("", encoding="utf-8")

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    work = tmp_path / "per_genome"
    work.mkdir()
    ref = tmp_path / "reference.fasta"
    ref.write_text(">ref\nACGT\n", encoding="utf-8")
    genome = tmp_path / "g1.fasta"
    genome.write_text(">g1\nACGA\n", encoding="utf-8")

    from repgenr.snptypers.base import SnpParams

    seq = mod._call_one(genome, ref, work, 4, SnpParams(), logging.getLogger("t"))
    assert seq == "ACGT"
    for cmd in calls:
        if cmd[0] == "minimap2":
            assert cmd[cmd.index("-t") + 1] == "4"
        if cmd[0] == "samtools":
            assert cmd[cmd.index("-@") + 1] == "4"
    pileup = [c for c in calls if c[:2] == ["bcftools", "mpileup"]][0]
    assert "-Ob" in pileup, "the pileup is written compressed"
    assert list(work.iterdir()) == [], "nothing is left behind once the consensus is read"
