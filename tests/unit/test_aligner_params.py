"""Aligner parameter passthrough: divergence-tuning flags reach the tool argv."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.aligners import progressivemauve as pm
from repgenr.aligners import sibeliaz as sz
from repgenr.aligners.base import AlignParams

_LOG = logging.getLogger("test")


def _genomes(tmp_path: Path, n: int) -> list[Path]:
    out = []
    for i in range(n):
        p = tmp_path / f"g{i}.fasta"
        p.write_text(">x\nACGTACGT\n")
        out.append(p)
    return out


def test_sibeliaz_forwards_tuning_flags(monkeypatch, tmp_path: Path) -> None:
    captured: dict[str, list[str]] = {}

    def fake_run(caps, cmd, **kw):
        captured["cmd"] = [str(c) for c in cmd]
        raise RuntimeError("stop after capture")

    monkeypatch.setattr(sz, "run_tool", fake_run)
    monkeypatch.setattr(sz, "_sibeliaz_invocation", lambda out_dir, logger: ["sibeliaz"])

    aln = sz.SibeliazAligner()
    params = AlignParams(threads=2, extra={"kmer": 15, "bubble": 200})
    with pytest.raises(RuntimeError):
        aln.align(_genomes(tmp_path, 3), None, tmp_path / "out", params, _LOG)

    cmd = captured["cmd"]
    assert "-k" in cmd and cmd[cmd.index("-k") + 1] == "15"
    assert "-b" in cmd and cmd[cmd.index("-b") + 1] == "200"
    # default: no -k when unset
    captured.clear()
    with pytest.raises(RuntimeError):
        aln.align(_genomes(tmp_path, 3), None, tmp_path / "out2", AlignParams(), _LOG)
    assert "-k" not in captured["cmd"]


def test_progressivemauve_forwards_seed_weight(monkeypatch, tmp_path: Path) -> None:
    captured: dict[str, list[str]] = {}

    def fake_run(caps, cmd, **kw):
        captured["cmd"] = [str(c) for c in cmd]
        raise RuntimeError("stop after capture")

    monkeypatch.setattr(pm, "run_tool", fake_run)

    aln = pm.ProgressiveMauveAligner()
    # 3 genomes (mauve needs >=3); threads=1 -> sequential -> exception propagates
    params = AlignParams(threads=1, extra={"seed_weight": 11})
    with pytest.raises(RuntimeError):
        aln.align(_genomes(tmp_path, 3), None, tmp_path / "out", params, _LOG)

    cmd = captured["cmd"]
    assert "--seed-weight" in cmd and cmd[cmd.index("--seed-weight") + 1] == "11"


# The upstream wrapper lines the macOS patch rewrites, verbatim.
_UPSTREAM_WRAPPER = """#!/bin/bash
outdir=$1
find $outdir -name "*.tmp" -printf "%p\\n" | sort > "$outdir/tmp_list.txt"
find $outdir -name "*.msa"  -print0 | sort -z | xargs -0 cat >> "$outdir/out.maf"
fasta=$(mktemp --suffix=.fa $outdir/block.XXXXX)
echo ">s" > $fasta
rm $fasta
"""


def _run_bsd_wrapper(monkeypatch, tmp_path: Path, work: Path) -> None:
    import subprocess
    import types

    upstream = tmp_path / "sibeliaz"
    upstream.write_text(_UPSTREAM_WRAPPER)
    monkeypatch.setattr(sz.sys, "platform", "darwin")
    monkeypatch.setattr(sz, "get_config", lambda: types.SimpleNamespace(active=False))
    monkeypatch.setattr(sz.shutil, "which", lambda name: str(upstream))
    argv = sz._sibeliaz_invocation(tmp_path, _LOG)
    assert argv[0] == "bash"
    subprocess.run([*argv, str(work)], check=True)


def test_bsd_wrapper_skips_appledouble_files(monkeypatch, tmp_path: Path) -> None:
    """On a non-HFS volume macOS writes a binary ``._`` sibling for every file;
    the wrapper's finds must not feed those into the block list or the MAF."""
    work = tmp_path / "work"
    work.mkdir()
    (work / "a.tmp").write_text("x")
    (work / "._a.tmp").write_bytes(b"\x00\x05\x16\x07Mac OS X\xb0")
    (work / "a.msa").write_text("a\ns x 0 4 + 4 ACGT\n")
    (work / "._a.msa").write_bytes(b"\x00\x05\x16\x07Mac OS X\xb0")
    _run_bsd_wrapper(monkeypatch, tmp_path, work)
    assert (work / "out.maf").read_bytes() == b"a\ns x 0 4 + 4 ACGT\n"
    assert (work / "tmp_list.txt").read_text().split() == [str(work / "a.tmp")]


def test_bsd_wrapper_leaves_no_block_temp_files(monkeypatch, tmp_path: Path) -> None:
    """BSD mktemp has no --suffix; the substitute must not leave the bare
    ``block.XXXXX`` file behind for every block it aligns."""
    work = tmp_path / "work"
    work.mkdir()
    _run_bsd_wrapper(monkeypatch, tmp_path, work)
    assert not list(work.glob("block.*"))
