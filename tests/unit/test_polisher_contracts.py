"""Polisher adapters: argument vectors, output contract, model detection, selection."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import pytest

from repgenr.assemblers.base import ReadSet
from repgenr.core.errors import UserInputError
from repgenr.polishers.base import PolishParams, registry, select_polisher

_LOG = logging.getLogger("test")
_POLISHED = ">c1\nACGTACGTACGA\n"


def _flag(cmd: list[str], flag: str) -> str:
    return cmd[cmd.index(flag) + 1]


@pytest.fixture()
def recorded(monkeypatch):
    """Record every tool call; write each polisher's output where the tool would."""
    import sys

    calls: list[dict] = []

    def run_tool(caps, command, *, logger, cwd=None, stdout_path=None, **kwargs):
        cmd = [str(c) for c in command]
        calls.append(
            {"cmd": cmd, "caps": caps.name, "mounts": list(kwargs.get("extra_mounts", ()))}
        )
        if cmd[0] == "medaka_consensus":
            out = Path(_flag(cmd, "-o"))
            out.mkdir(parents=True, exist_ok=True)
            (out / "consensus.fasta").write_text(_POLISHED, encoding="utf-8")
        elif cmd[0] == "minimap2":
            Path(_flag(cmd, "-o")).write_text("paf\n", encoding="utf-8")
        elif cmd[0] == "racon":
            assert stdout_path is not None, "racon writes the consensus to stdout"
            Path(stdout_path).write_text(_POLISHED, encoding="utf-8")
        return 0

    for name in registry.names():
        module = sys.modules[registry.get(name).__module__]
        monkeypatch.setattr(module, "run_tool", run_tool, raising=False)
    return calls


def _ont(tmp_path: Path, header: str = "@read1 ch=1") -> tuple[ReadSet, Path]:
    fq = tmp_path / "SRR1_1.fastq.gz"
    with gzip.open(fq, "wt", encoding="utf-8") as fo:
        fo.write(f"{header}\nACGT\n+\nIIII\n")
    draft = tmp_path / "draft.fa"
    draft.write_text(">c1\nACGTACGTACGT\n", encoding="utf-8")
    return ReadSet("SRR1", "OXFORD_NANOPORE", "GridION", "SINGLE", (fq,), 4), draft


def test_medaka_runs_consensus_with_the_bacterial_model_and_auto_detection(recorded, tmp_path):
    reads, draft = _ont(
        tmp_path, "@read1 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
    )
    result = registry.create("medaka").polish(
        reads, draft, tmp_path / "pol", PolishParams(threads=6), _LOG
    )
    assert result.contigs.read_text(encoding="utf-8") == _POLISHED and result.rounds == 1
    [call] = recorded
    cmd = call["cmd"]
    assert cmd[0] == "medaka_consensus" and _flag(cmd, "-t") == "6" and "--bacteria" in cmd
    assert _flag(cmd, "-i") == str(reads.files[0]) and _flag(cmd, "-d") == str(draft)
    assert "-m" not in cmd  # the model comes from the read headers
    assert result.tool_stats["model_source"] == "read headers"
    assert str(tmp_path.resolve()) in call["mounts"]


def test_medaka_assumes_the_bacterial_model_for_sra_reads(recorded, tmp_path, caplog):
    """SRA rewrites FASTQ headers, so mirrored reads never name their basecaller."""
    from repgenr.polishers.medaka import BACTERIAL_MODEL

    reads, draft = _ont(tmp_path, "@SRR1.1 32921a8b-d4e2-417a-bf31-034a6749932f/1")
    with caplog.at_level(logging.WARNING):
        result = registry.create("medaka").polish(
            reads, draft, tmp_path / "pol", PolishParams(threads=2), _LOG
        )
    cmd = recorded[-1]["cmd"]
    assert _flag(cmd, "-m") == BACTERIAL_MODEL and "--bacteria" not in cmd
    assert result.tool_stats["model_source"].startswith("assumed")
    assert any("basecaller" in r.getMessage() for r in caplog.records)
    result = registry.create("medaka").polish(
        reads, draft, tmp_path / "pol2", PolishParams(threads=2, extra={"bacteria": "false"}), _LOG
    )
    assert "-m" not in recorded[-1]["cmd"] and result.tool_stats["model_source"] == "medaka default"


def test_medaka_takes_an_explicit_model(recorded, tmp_path):
    reads, draft = _ont(tmp_path)
    result = registry.create("medaka").polish(
        reads,
        draft,
        tmp_path / "pol2",
        PolishParams(threads=2, extra={"model": "r941_min_hac_g507"}),
        _LOG,
    )
    assert _flag(recorded[-1]["cmd"], "-m") == "r941_min_hac_g507"
    assert result.tool_stats["model_source"] == "--tool-arg"


def _two_files(tmp_path: Path) -> tuple[ReadSet, Path]:
    """An ONT run that ENA lists as two FASTQ files (labelled PAIRED)."""
    reads, draft = _ont(tmp_path)
    second = tmp_path / "SRR1_2.fastq.gz"
    with gzip.open(second, "wt", encoding="utf-8") as fo:
        fo.write("@read2\nTTTT\n+\nIIII\n")
    return ReadSet("SRR1", "OXFORD_NANOPORE", "GridION", "PAIRED", (*reads.files, second), 8), draft


@pytest.mark.parametrize("name", ["medaka", "racon"])
def test_polishers_use_every_read_file_of_a_run(recorded, tmp_path, name):
    """A long-read run listed as several files is polished with all of them, in one file."""
    reads, draft = _two_files(tmp_path)
    registry.create(name).polish(reads, draft, tmp_path / "pol", PolishParams(), _LOG)
    cmd = recorded[0]["cmd"]
    query = Path(_flag(cmd, "-i")) if name == "medaka" else Path(cmd[-1])
    assert query not in reads.files
    with gzip.open(query, "rt", encoding="utf-8") as fi:
        assert [line for line in fi.read().splitlines() if line.startswith("@")] == [
            "@read1 ch=1",
            "@read2",
        ]


def test_read_files_of_mixed_compression_are_refused(recorded, tmp_path):
    reads, draft = _ont(tmp_path)
    plain = tmp_path / "SRR1_2.fastq"
    plain.write_text("@r\nA\n+\nI\n", encoding="utf-8")
    mixed = ReadSet("SRR1", "OXFORD_NANOPORE", "GridION", "SINGLE", (*reads.files, plain), 4)
    with pytest.raises(UserInputError, match="compressed"):
        registry.create("medaka").polish(mixed, draft, tmp_path / "pol", PolishParams(), _LOG)


def test_racon_alternates_minimap2_and_racon_for_each_round(recorded, tmp_path):
    reads, draft = _ont(tmp_path)
    result = registry.create("racon").polish(
        reads, draft, tmp_path / "pol", PolishParams(threads=4, rounds=2), _LOG
    )
    assert result.rounds == 2 and result.contigs.read_text(encoding="utf-8") == _POLISHED
    tools = [c["cmd"][0] for c in recorded]
    assert tools == ["minimap2", "racon", "minimap2", "racon"]
    first_map, first_racon = recorded[0]["cmd"], recorded[1]["cmd"]
    assert _flag(first_map, "-x") == "map-ont" and _flag(first_map, "-t") == "4"
    assert first_map[-2:] == [str(draft), str(reads.files[0])]  # target then query
    assert first_racon[-3:] == [str(reads.files[0]), _flag(first_map, "-o"), str(draft)]
    # round two polishes round one's consensus, and each tool runs in its own image
    assert recorded[2]["cmd"][-2] != str(draft)
    assert {c["caps"] for c in recorded} == {"minimap2", "racon"}


def test_racon_uses_the_pacbio_preset_for_clr_reads(recorded, tmp_path):
    reads, draft = _ont(tmp_path)
    clr = ReadSet("SRR1", "PACBIO_SMRT", "PacBio RS II", "SINGLE", reads.files, 4)
    registry.create("racon").polish(clr, draft, tmp_path / "pol", PolishParams(), _LOG)
    assert _flag(recorded[0]["cmd"], "-x") == "map-pb"


@pytest.mark.parametrize(
    ("platform", "model", "expected"),
    [
        ("OXFORD_NANOPORE", "GridION", "medaka"),
        ("PACBIO_SMRT", "PacBio RS II", "racon"),
        ("PACBIO_SMRT", "Sequel II", None),  # HiFi needs no polishing
        ("PACBIO_SMRT", "Revio", None),
        ("ILLUMINA", "Illumina MiSeq", None),
    ],
)
def test_select_polisher_by_platform_and_instrument(monkeypatch, platform, model, expected):
    monkeypatch.setattr("repgenr.polishers.base.tool_available", lambda caps: True)
    reads = ReadSet("SRR1", platform, model, "SINGLE", (Path("r.fq.gz"),), 4)
    assert select_polisher(registry, reads) == expected


def test_polishers_are_a_registered_family_and_listed() -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    assert {"medaka", "racon"} <= set(registry.names())
    result = CliRunner().invoke(app, ["list-tools"])
    assert result.exit_code == 0 and "polishers: medaka, racon" in result.output


def test_medaka_retries_without_bacteria_when_medaka_refuses(monkeypatch, tmp_path, caplog):
    import repgenr.polishers.medaka as md
    from repgenr.core.errors import ToolExecutionError

    calls: list[list[str]] = []

    def run_tool(caps, command, *, logger, **kwargs):
        cmd = [str(c) for c in command]
        calls.append(cmd)
        if "--bacteria" in cmd:
            raise ToolExecutionError(
                cmd, 1, output="ERROR: --bacteria was specified but input data was not compatible."
            )
        out = Path(_flag(cmd, "-o"))
        out.mkdir(parents=True, exist_ok=True)
        (out / "consensus.fasta").write_text(_POLISHED, encoding="utf-8")
        return 0

    monkeypatch.setattr(md, "run_tool", run_tool)
    reads, draft = _ont(tmp_path, "@r1 basecall_model_version_id=dna_r9.4.1_e8_sup@v3.3")
    with caplog.at_level(logging.WARNING):
        result = registry.create("medaka").polish(
            reads, draft, tmp_path / "pol", PolishParams(threads=2), _LOG
        )
    assert [("--bacteria" in c) for c in calls] == [True, False]
    assert result.tool_stats["bacteria"] is False and result.contigs.exists()
    assert any("no bacterial model" in r.getMessage() for r in caplog.records)


def test_racon_records_the_minimap2_image_under_a_container_backend(monkeypatch) -> None:
    # racon runs minimap2 from its own pinned image; only the racon image
    # reached the stage record and `versions`.
    from repgenr.core import containers
    from repgenr.polishers.racon import _MINIMAP2, RaconPolisher

    monkeypatch.setattr(containers, "check_engine_ready", lambda config=None: None)
    monkeypatch.setattr(
        "repgenr.core.plugins.check_binaries", lambda specs: {s.name: "x" for s in specs}
    )
    try:
        containers.configure_container("docker")
        versions = RaconPolisher().preflight()
    finally:
        containers.configure_container("none")
    assert versions == {
        "racon": RaconPolisher.capabilities.container,
        "minimap2": _MINIMAP2.container,
    }
