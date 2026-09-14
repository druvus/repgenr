"""Assembler adapters: argument vectors, output contract, platform mapping, selection."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.assemblers.base import (
    AssembleParams,
    ReadSet,
    registry,
    select_assembler,
)
from repgenr.core.errors import UserInputError

_LOG = logging.getLogger("test")
_CONTIGS = ">c1\nACGTACGTACGT\n>c2\nGGCC\n"


def _flag(cmd: list[str], flag: str) -> str:
    return cmd[cmd.index(flag) + 1]


mounts: list[list[str]] = []


def _fake_run_tool(recorded: list[list[str]]):
    """Write each assembler's contigs where the real tool would."""

    def run_tool(caps, command, *, logger, cwd=None, **kwargs):
        cmd = [str(c) for c in command]
        recorded.append(cmd)
        mounts.append([str(m) for m in kwargs.get("extra_mounts", ())])
        if cmd[0] == "skesa":
            Path(_flag(cmd, "--contigs_out")).write_text(_CONTIGS, encoding="utf-8")
        elif cmd[0] == "shovill":
            out = Path(_flag(cmd, "--outdir"))
            out.mkdir(parents=True, exist_ok=True)
            (out / "contigs.fa").write_text(_CONTIGS, encoding="utf-8")
        elif cmd[0] == "flye":
            out = Path(_flag(cmd, "-o"))
            out.mkdir(parents=True, exist_ok=True)
            (out / "assembly.fasta").write_text(_CONTIGS, encoding="utf-8")
        return 0

    return run_tool


@pytest.fixture()
def recorded(monkeypatch, tmp_path) -> list[list[str]]:
    import sys

    calls: list[list[str]] = []
    fake = _fake_run_tool(calls)
    for name in registry.names():
        module = sys.modules[registry.get(name).__module__]
        monkeypatch.setattr(module, "run_tool", fake, raising=False)
    return calls


def _reads(tmp_path: Path, platform="ILLUMINA", layout="PAIRED", model="Illumina MiSeq") -> ReadSet:
    files = []
    for i in (1, 2) if layout == "PAIRED" else (1,):
        f = tmp_path / f"run_{i}.fastq.gz"
        f.write_bytes(b"x")
        files.append(f)
    return ReadSet("SRR1", platform, model, layout, tuple(files), bases=1_000_000)


_TOKENS = {
    "skesa": ["--cores", "5", "--memory", "7"],
    "shovill": ["--cpus", "5", "--ram", "8", "--assembler", "spades"],  # 7 GB raised to 8
    "flye": ["--threads", "5", "--nano-hq"],
}


@pytest.mark.parametrize("tool", sorted(_TOKENS))
def test_assembler_contract(tool, recorded, tmp_path) -> None:
    adapter = registry.create(tool)
    reads = (
        _reads(tmp_path, "OXFORD_NANOPORE", "SINGLE", "GridION")
        if tool == "flye"
        else _reads(tmp_path)
    )
    result = adapter.assemble(reads, tmp_path / "out", AssembleParams(threads=5, memory_gb=7), _LOG)
    assert result.contigs.exists() and result.contigs.read_text(encoding="utf-8") == _CONTIGS
    flat = [tok for cmd in recorded for tok in cmd]
    for token in _TOKENS[tool]:
        assert token in flat, f"{token!r} missing from argv for {tool}"
    # the reads' directory is declared for the container backend's mounts
    assert str(tmp_path.resolve()) in mounts[-1]


def test_every_builtin_assembler_has_contract_coverage() -> None:
    builtin = {"skesa", "shovill", "flye"}
    assert builtin <= set(registry.names())
    assert builtin & set(registry.names()) <= set(_TOKENS)


def test_skesa_single_end_passes_one_file(recorded, tmp_path) -> None:
    registry.create("skesa").assemble(
        _reads(tmp_path, layout="SINGLE"), tmp_path / "out", AssembleParams(), _LOG
    )
    cmd = recorded[0]
    assert _flag(cmd, "--reads").endswith("run_1.fastq.gz") and "," not in _flag(cmd, "--reads")


def test_shovill_refuses_single_end(recorded, tmp_path) -> None:
    with pytest.raises(UserInputError, match="paired"):
        registry.create("shovill").assemble(
            _reads(tmp_path, layout="SINGLE"), tmp_path / "out", AssembleParams(), _LOG
        )


@pytest.mark.parametrize(
    ("platform", "model", "extra", "expected"),
    [
        ("OXFORD_NANOPORE", "GridION", {}, "--nano-hq"),
        ("OXFORD_NANOPORE", "MinION", {"mode": "nano-raw"}, "--nano-raw"),
        ("PACBIO_SMRT", "Sequel II", {}, "--pacbio-hifi"),
        ("PACBIO_SMRT", "Revio", {}, "--pacbio-hifi"),
        ("PACBIO_SMRT", "PacBio RS II", {}, "--pacbio-raw"),
    ],
)
def test_flye_read_mode_follows_platform_and_instrument(
    recorded, tmp_path, platform, model, extra, expected
) -> None:
    registry.create("flye").assemble(
        _reads(tmp_path, platform, "SINGLE", model),
        tmp_path / "out",
        AssembleParams(extra=extra),
        _LOG,
    )
    assert expected in recorded[0]


def test_select_assembler_by_platform_and_layout(monkeypatch) -> None:
    from repgenr.core import plugins

    monkeypatch.setattr(plugins, "_tool_available", lambda caps: True)
    assert select_assembler(registry, ReadSet("r", "ILLUMINA", "MiSeq", "PAIRED", ())) == "skesa"
    assert select_assembler(registry, ReadSet("r", "ILLUMINA", "MiSeq", "SINGLE", ())) == "skesa"
    assert (
        select_assembler(registry, ReadSet("r", "OXFORD_NANOPORE", "GridION", "SINGLE", ()))
        == "flye"
    )
    assert select_assembler(registry, ReadSet("r", "PACBIO_SMRT", "Revio", "SINGLE", ())) == "flye"
    assert select_assembler(registry, ReadSet("r", "ION_TORRENT", "S5", "SINGLE", ())) is None


def test_adapters_declare_read_types_and_layouts() -> None:
    assert registry.get("skesa").read_types == frozenset({"ILLUMINA"})
    assert registry.get("shovill").layouts == frozenset({"PAIRED"})
    assert registry.get("flye").read_types == frozenset({"OXFORD_NANOPORE", "PACBIO_SMRT"})


def test_shovill_raises_the_ram_cap_to_its_minimum(recorded, tmp_path, caplog) -> None:
    import logging as _logging

    with caplog.at_level(_logging.WARNING):
        registry.create("shovill").assemble(
            _reads(tmp_path), tmp_path / "out", AssembleParams(memory_gb=6), _LOG
        )
    assert _flag(recorded[0], "--ram") == "8"
    assert any("refuses a RAM cap below 8" in r.message for r in caplog.records)
