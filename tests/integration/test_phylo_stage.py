"""Phylo stage composition tests using in-process fakes (no external tools)."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.aligners.base import Aligner, AlignResult
from repgenr.aligners.base import registry as aligner_registry
from repgenr.core.context import WorkdirContext
from repgenr.core.plugins import ToolCapabilities
from repgenr.snptypers.base import SnpResult, SnpTyper
from repgenr.snptypers.base import registry as snp_registry
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
    seen_extra: dict | None = None

    def preflight(self):
        return {"faketree": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        type(self).seen_extra = dict(params.extra)
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "tree.nwk"
        tree.write_text("(from_msa);\n")
        return tree


class _FakeAligner(Aligner):
    capabilities = ToolCapabilities(name="fakealigner", accepted_extras=frozenset({"kmer"}))
    seen_extra: dict | None = None

    def preflight(self):
        return {"fakealigner": "1.0"}

    def align(self, genomes, reference, out_dir, params, logger) -> AlignResult:
        type(self).seen_extra = dict(params.extra)
        out_dir.mkdir(parents=True, exist_ok=True)
        msa = out_dir / "msa.fasta"
        msa.write_text("".join(f">{Path(g).stem}\nACGT\n" for g in genomes))
        return AlignResult(msa_fasta=msa)


class _FakeSnpTyper(SnpTyper):
    capabilities = ToolCapabilities(name="fakesnptyper", accepted_extras=frozenset({"kmer"}))
    requires_reference = False
    seen_extra: dict | None = None

    def preflight(self):
        return {"fakesnptyper": "1.0"}

    def call(self, genomes, reference, out_dir, params, logger) -> SnpResult:
        type(self).seen_extra = dict(params.extra)
        out_dir.mkdir(parents=True, exist_ok=True)
        core = out_dir / "core.fasta"
        core.write_text("".join(f">{Path(g).stem}\nACGT\n" for g in genomes))
        return SnpResult(core_snp_fasta=core)


@pytest.fixture
def fake_phylo_tools():
    tb_registry._load()
    aligner_registry._load()
    tb_registry.register("faketree_genomes", _GenomesTreeBuilder, replace=True)
    tb_registry.register("faketree_msa", _MsaTreeBuilder, replace=True)
    aligner_registry.register("fakealigner", _FakeAligner, replace=True)
    yield
    for n in ("faketree_genomes", "faketree_msa"):
        tb_registry._classes.pop(n, None)
    aligner_registry._classes.pop("fakealigner", None)


@pytest.fixture
def fake_snptyper():
    snp_registry._load()
    snp_registry.register("fakesnptyper", _FakeSnpTyper, replace=True)
    yield
    snp_registry._classes.pop("fakesnptyper", None)


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


def test_aligner_receives_extra(workdir: Path, fake_phylo_tools) -> None:
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    _FakeAligner.seen_extra = None
    run(
        ctx,
        PhyloParams(
            treebuilder="faketree_msa",
            msa_source="aligner",
            aligner="fakealigner",
            no_outgroup=True,
            extra={"kmer": "15"},
        ),
    )
    assert _FakeAligner.seen_extra == {"kmer": "15"}


def test_mask_key_stripped_before_aligner_and_treebuilder(
    workdir: Path, fake_phylo_tools
) -> None:
    """'mask' is a phylo-stage-owned key: the aligner and tree builder must
    never see it, even though it rides along in params.extra."""
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    _FakeAligner.seen_extra = None
    _MsaTreeBuilder.seen_extra = None
    run(
        ctx,
        PhyloParams(
            treebuilder="faketree_msa",
            msa_source="aligner",
            aligner="fakealigner",
            no_outgroup=True,
            extra={"seed_weight": 11, "mask": "gubbins"},
        ),
    )
    assert _FakeAligner.seen_extra == {"seed_weight": 11}
    assert _MsaTreeBuilder.seen_extra == {"seed_weight": 11}


def _extras_warnings(caplog, key: str) -> list[str]:
    return [
        r.getMessage()
        for r in caplog.records
        if r.levelno == logging.WARNING and key in r.getMessage()
    ]


def test_aligner_key_does_not_warn_from_the_tree_builder(
    workdir: Path, fake_phylo_tools, caplog
) -> None:
    """An extra the aligner reads is consumed by the stage, so the tree builder
    (which declares no extras of its own) must not report it as ignored."""
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=logging.getLogger("test-phylo"))
    with caplog.at_level(logging.WARNING, logger="test-phylo"):
        run(
            ctx,
            PhyloParams(
                treebuilder="faketree_msa",
                msa_source="aligner",
                aligner="fakealigner",
                no_outgroup=True,
                extra={"kmer": "15"},
            ),
        )
    assert _extras_warnings(caplog, "kmer") == []


def test_extra_read_by_no_tool_warns_once(workdir: Path, fake_phylo_tools, caplog) -> None:
    """A key neither the aligner nor the tree builder declares is reported once."""
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=logging.getLogger("test-phylo"))
    with caplog.at_level(logging.WARNING, logger="test-phylo"):
        run(
            ctx,
            PhyloParams(
                treebuilder="faketree_msa",
                msa_source="aligner",
                aligner="fakealigner",
                no_outgroup=True,
                extra={"nosuchkey": "1"},
            ),
        )
    assert len(_extras_warnings(caplog, "nosuchkey")) == 1


def test_snptyper_key_does_not_warn_from_the_tree_builder(
    workdir: Path, fake_phylo_tools, fake_snptyper, caplog
) -> None:
    """Same for the snptype MSA source: the typer's key is consumed, not ignored."""
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=logging.getLogger("test-phylo"))
    with caplog.at_level(logging.WARNING, logger="test-phylo"):
        run(
            ctx,
            PhyloParams(
                treebuilder="faketree_msa",
                msa_source="snptype",
                snptyper="fakesnptyper",
                no_outgroup=True,
                extra={"kmer": "15"},
            ),
        )
    assert _extras_warnings(caplog, "kmer") == []


def test_snptype_receives_extra_without_mask(
    workdir: Path, fake_phylo_tools, fake_snptyper
) -> None:
    """A non-mask extra key (e.g. from --aligner-arg) must reach the SNP
    typer's params.extra when msa_source=snptype; 'mask' must not (it is
    consumed by the phylo stage itself, via SnptypeParams.mask)."""
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    _FakeSnpTyper.seen_extra = None
    run(
        ctx,
        PhyloParams(
            treebuilder="faketree_msa",
            msa_source="snptype",
            snptyper="fakesnptyper",
            no_outgroup=True,
            extra={"kmer": "15", "mask": "none"},
        ),
    )
    assert _FakeSnpTyper.seen_extra == {"kmer": "15"}
