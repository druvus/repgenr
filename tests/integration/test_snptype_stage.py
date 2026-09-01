"""SNP typing stage test using an in-process fake SNP typer."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import CORE_SNP_FASTA
from repgenr.core.errors import UserInputError
from repgenr.core.plugins import ToolCapabilities
from repgenr.maskers.base import Masker
from repgenr.maskers.base import registry as masker_registry
from repgenr.snptypers.base import SnpResult, SnpTyper, registry
from repgenr.stages.snptype import SnptypeParams, run


class _FakeTyper(SnpTyper):
    capabilities = ToolCapabilities(name="faketyper")
    requires_reference = False

    def preflight(self) -> dict[str, str]:
        return {"faketyper": "1.0"}

    def call(self, genomes, reference, out_dir, params, logger) -> SnpResult:  # noqa: ANN001
        core = out_dir / "core.fasta"
        core.write_text("".join(f">{g.stem}\nACGT\n" for g in genomes))
        return SnpResult(core_snp_fasta=core)


@pytest.fixture
def fake_typer():
    registry._load()
    registry.register("faketyper", _FakeTyper, replace=True)
    yield
    registry._classes.pop("faketyper", None)


def test_snptype_writes_core_snp(workdir: Path, genome_files, fake_typer) -> None:
    ctx = WorkdirContext(workdir, create=True)
    result = run(ctx, SnptypeParams(tool="faketyper", all_genomes=True))

    assert result.core_snp_fasta == ctx.snp_dir / CORE_SNP_FASTA
    assert result.core_snp_fasta.exists()
    # one record per input genome
    assert result.core_snp_fasta.read_text().count(">") == len(genome_files)
    assert ctx.config.stages["snptype"].tool == "faketyper"


class _FullTyper(_FakeTyper):
    def call(self, genomes, reference, out_dir, params, logger) -> SnpResult:  # noqa: ANN001
        core = out_dir / "core.fasta"
        full = out_dir / "full.fasta"
        core.write_text("".join(f">{g.stem}\nACGT\n" for g in genomes))
        full.write_text("".join(f">{g.stem}\nACGTACGTACGT\n" for g in genomes))
        return SnpResult(core_snp_fasta=core, full_alignment=full)


class _RecordingMasker(Masker):
    capabilities = ToolCapabilities(name="fakemask")
    seen: dict = {}

    def preflight(self) -> dict[str, str]:
        return {"fakemask": "1"}

    def mask(self, full_alignment, out_dir, params, logger):  # noqa: ANN001
        _RecordingMasker.seen = {"input": full_alignment, "threads": params.threads}
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / "masked.fasta"
        out.write_text(">a\nAC\n")
        return out


@pytest.fixture
def fake_masker(register_tool):
    register_tool(masker_registry, "fakemask", _RecordingMasker)
    yield


def test_masker_receives_full_alignment(workdir, genome_files, register_tool, fake_masker):
    register_tool(registry, "fulltyper", _FullTyper)
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, SnptypeParams(tool="fulltyper", all_genomes=True, mask="fakemask", threads=3))
    assert _RecordingMasker.seen["input"].name == "full.fasta"
    assert _RecordingMasker.seen["threads"] == 3
    assert (ctx.snp_dir / CORE_SNP_FASTA).read_text() == ">a\nAC\n"


def test_mask_refused_without_full_alignment(workdir, genome_files, fake_typer, fake_masker):
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(UserInputError, match="whole-genome alignment"):
        run(ctx, SnptypeParams(tool="faketyper", all_genomes=True, mask="fakemask"))
