"""SNP typing stage test using an in-process fake SNP typer."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import CORE_SNP_FASTA
from repgenr.core.plugins import ToolCapabilities
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
    registry._classes["faketyper"] = _FakeTyper
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
