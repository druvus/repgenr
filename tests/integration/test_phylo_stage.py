"""Phylo stage composition tests using in-process fakes (no external tools)."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.aligners.base import Aligner, AlignResult
from repgenr.aligners.base import registry as aligner_registry
from repgenr.core.context import WorkdirContext
from repgenr.core.plugins import ToolCapabilities
from repgenr.stages.phylo import PhyloParams, run
from repgenr.treebuilders.base import InputKind, TreeBuilder
from repgenr.treebuilders.base import registry as tb_registry


class _GenomesTreeBuilder(TreeBuilder):
    capabilities = ToolCapabilities(name="faketree_genomes")
    input_kind = InputKind.GENOMES

    def preflight(self):
        return {"faketree": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "tree.nwk"
        leaves = [Path(g).stem for g in msa_or_genomes]
        tree.write_text("(" + ",".join(leaves) + ");\n")
        return tree


class _MsaTreeBuilder(TreeBuilder):
    capabilities = ToolCapabilities(name="faketree_msa")
    input_kind = InputKind.MSA_FASTA

    def preflight(self):
        return {"faketree": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "tree.nwk"
        tree.write_text("(from_msa);\n")
        return tree


class _FakeAligner(Aligner):
    capabilities = ToolCapabilities(name="fakealigner")

    def preflight(self):
        return {"fakealigner": "1.0"}

    def align(self, genomes, reference, out_dir, params, logger) -> AlignResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        msa = out_dir / "msa.fasta"
        msa.write_text("".join(f">{Path(g).stem}\nACGT\n" for g in genomes))
        return AlignResult(msa_fasta=msa)


@pytest.fixture
def fake_phylo_tools():
    tb_registry._load()
    aligner_registry._load()
    tb_registry._classes["faketree_genomes"] = _GenomesTreeBuilder
    tb_registry._classes["faketree_msa"] = _MsaTreeBuilder
    aligner_registry._classes["fakealigner"] = _FakeAligner
    yield
    for n in ("faketree_genomes", "faketree_msa"):
        tb_registry._classes.pop(n, None)
    aligner_registry._classes.pop("fakealigner", None)


def _make_reps(workdir: Path) -> None:
    reps = workdir / "derep" / "representatives"
    reps.mkdir(parents=True)
    for i in range(1, 4):
        (reps / f"Fam_gen_sp_GCA_00000{i}.fasta").write_text(f">s{i}\nACGTACGT\n")


def test_alignment_free_path(workdir: Path, fake_phylo_tools) -> None:
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    tree = run(ctx, PhyloParams(treebuilder="faketree_genomes", no_outgroup=True))
    assert tree.exists()
    assert tree.read_text().startswith("(")
    assert ctx.config.stages["phylo"].tool == "faketree_genomes"


def test_aligner_msa_path(workdir: Path, fake_phylo_tools) -> None:
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    tree = run(
        ctx,
        PhyloParams(
            treebuilder="faketree_msa",
            msa_source="aligner",
            aligner="fakealigner",
            no_outgroup=True,
        ),
    )
    assert tree.read_text().strip() == "(from_msa);"
    assert (ctx.align_dir / "msa.fasta").exists()
