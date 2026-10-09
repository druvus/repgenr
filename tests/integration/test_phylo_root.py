"""phylo --root: rooting of tree/tree.nwk, and how tree2tax reads the result.

A fake genome-input builder writes a fixed Newick with branch lengths and
support values over four ingroup genomes and the outgroup, so each rooting
mode can be checked without an external tool.
"""

from __future__ import annotations

import logging
from pathlib import Path

import dendropy
import pytest
from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import TREE_NWK, TREE_UNROOTED_NWK, record_name
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.core.plugins import ToolCapabilities
from repgenr.stages.phylo import PhyloBuildParams, PhyloParams, phylo_build
from repgenr.stages.phylo import run as phylo_run
from repgenr.stages.tree2tax import Tree2taxParams, Tree2taxStepParams, tree2tax_relations
from repgenr.stages.tree2tax import run as tree2tax_run
from repgenr.tree.rooting import (
    is_rooted_on,
    midpoint_root,
    read_newick,
    root_on_outgroup,
    write_newick,
)
from repgenr.treebuilders.base import InputKind, TreeBuilder
from repgenr.treebuilders.base import registry as tb_registry

_LOG = logging.getLogger("test")
_OUTGROUP = "Fam_Out_grp_GCF_000009.1"
_INGROUP = [f"Fam_gen_sp_GCA_00000{i}.1" for i in range(1, 5)]

# Unrooted (trifurcating root), with supports 95 and 88 on the two internal
# branches and lengths on every branch: {A,B}|{C,D,OG} has 95, {A,B,C}|{D,OG} 88.
_UNROOTED = "({a}:0.1,{b}:0.2,({c}:0.3,({d}:0.4,{og}:1.0)88:0.05)95:0.06);\n"
# As IQ-TREE writes with -o: the outgroup is a child of a trifurcating root.
_IQTREE_STYLE = "({og}:0.5,{a}:0.1,(({b}:0.2,{c}:0.3)70:0.04,{d}:0.4)80:0.06);\n"
# Already rooted on the outgroup: bifurcating root, outgroup on one side.
_ROOTED = "({og}:0.5,(({a}:0.1,{b}:0.2)95:0.06,({c}:0.3,{d}:0.4)88:0.05):0.5);\n"


def _fill(template: str, names: list[str]) -> str:
    a, b, c, d, og = names
    return template.format(a=a, b=b, c=c, d=d, og=og)


class _FixedTreeBuilder(TreeBuilder):
    capabilities = ToolCapabilities(name="faketree_fixed")
    input_kind = InputKind.GENOMES
    template = _UNROOTED
    calls = 0

    def preflight(self):
        return {"faketree": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        type(self).calls += 1
        names = [record_name(g) for g in msa_or_genomes]
        if len(names) == 4:  # built without an outgroup
            a, b, c, d = names
            text = f"({a}:0.1,{b}:0.2,({c}:0.3,{d}:0.4)95:0.06);\n"
        else:
            text = _fill(type(self).template, names)
        out = out_dir / "builder.nwk"
        out.write_text(text, encoding="utf-8")
        return out


@pytest.fixture
def fixed_builder(register_tool):
    tb_registry._load()
    register_tool(tb_registry, "faketree_fixed", _FixedTreeBuilder)
    _FixedTreeBuilder.template = _UNROOTED
    _FixedTreeBuilder.calls = 0
    yield _FixedTreeBuilder
    _FixedTreeBuilder.template = _UNROOTED


def _workdir(workdir: Path, *, outgroup: bool = True) -> WorkdirContext:
    """A workdir with four representatives and, unless ``outgroup`` is False,
    a staged outgroup; the stage logs go to the 'test-root' logger."""
    reps = workdir / "derep" / "representatives"
    reps.mkdir(parents=True)
    for name in _INGROUP:
        (reps / f"{name}.fasta").write_text(">s\nACGTACGT\n")
    if outgroup:
        og = workdir / "outgroup"
        og.mkdir()
        (og / f"{_OUTGROUP}.fasta").write_text(">o\nACGTACGT\n")
        (workdir / "outgroup_accession.txt").write_text("GCF_000009.1\n")
    return WorkdirContext(workdir, create=True, logger=logging.getLogger("test-root"))


def _splits(text: str) -> dict[frozenset[str], str]:
    """Support label by split, each split given by the side without the first leaf."""
    tree = read_newick(text)
    leaves = sorted(lf.taxon.label for lf in tree.leaf_node_iter())
    first = leaves[0]
    out: dict[frozenset[str], str] = {}
    for node in tree.internal_nodes():
        if node is tree.seed_node or node.label is None:
            continue
        side = frozenset(lf.taxon.label for lf in node.leaf_iter())
        if first in side:
            side = frozenset(leaves) - side
        out[side] = node.label
    return out


def _total_length(text: str) -> float:
    tree = read_newick(text)
    return sum(e.length or 0.0 for e in tree.postorder_edge_iter() if e.tail_node is not None)


def _root_children(text: str) -> list[dendropy.Node]:
    return read_newick(text).seed_node.child_nodes()


# --- the shared rooting functions -------------------------------------------------

_NAMES = [*_INGROUP, _OUTGROUP]


def test_outgroup_rooting_keeps_supports_on_their_splits_and_the_lengths() -> None:
    text = _fill(_UNROOTED, _NAMES)
    tree = read_newick(text)
    assert root_on_outgroup(tree, _OUTGROUP, source="t") is True
    rooted = write_newick(tree)
    assert is_rooted_on(read_newick(rooted), _OUTGROUP)
    assert _splits(rooted) == _splits(text)
    assert _total_length(rooted) == pytest.approx(_total_length(text))
    # The outgroup's branch is split in two halves.
    lengths = {
        (n.taxon.label if n.taxon else "ingroup"): n.edge.length for n in _root_children(rooted)
    }
    assert lengths == {_OUTGROUP: pytest.approx(0.5), "ingroup": pytest.approx(0.5)}


def test_outgroup_rooting_is_a_no_op_on_a_tree_rooted_on_the_outgroup() -> None:
    text = _fill(_ROOTED, _NAMES)
    tree = read_newick(text)
    assert root_on_outgroup(tree, _OUTGROUP, source="t") is False
    assert write_newick(tree) == write_newick(read_newick(text))


def test_outgroup_that_is_not_a_leaf_is_an_error() -> None:
    tree = read_newick(_fill(_UNROOTED, _NAMES))
    with pytest.raises(WorkdirError, match="not a leaf"):
        root_on_outgroup(tree, "absent", source="t", hint=" Hint.")


def test_midpoint_rooting_keeps_supports_and_lengths() -> None:
    text = _fill(_UNROOTED, _NAMES)
    tree = read_newick(text)
    midpoint_root(tree, source="t")
    rooted = write_newick(tree)
    assert len(_root_children(rooted)) == 2
    assert _splits(rooted) == _splits(text)
    assert _total_length(rooted) == pytest.approx(_total_length(text))
    # The longest path runs from the outgroup (1.0 + 0.05) to B (0.2 + 0.06);
    # its midpoint lies on the outgroup's branch.
    assert is_rooted_on(read_newick(rooted), _OUTGROUP)


def test_midpoint_rooting_needs_branch_lengths() -> None:
    tree = read_newick("(a,b,(c,d));")
    with pytest.raises(WorkdirError, match="branch lengths"):
        midpoint_root(tree, source="t")


def test_written_tree_has_no_rooting_tag_and_unquoted_underscores() -> None:
    tree = read_newick(_fill(_UNROOTED, _NAMES))
    root_on_outgroup(tree, _OUTGROUP, source="t")
    text = write_newick(tree)
    assert "[&R]" not in text and "'" not in text
    assert _OUTGROUP in text


# --- phylo ---------------------------------------------------------------------


def test_default_roots_on_the_outgroup_and_keeps_the_builder_output(
    workdir: Path, fixed_builder
) -> None:
    ctx = _workdir(workdir)
    tree = phylo_run(ctx, PhyloParams(treebuilder="faketree_fixed"))
    rooted = tree.read_text(encoding="utf-8")
    assert is_rooted_on(read_newick(rooted), _OUTGROUP)
    unrooted = (ctx.tree_dir / TREE_UNROOTED_NWK).read_text(encoding="utf-8")
    assert unrooted == _fill(_UNROOTED, _NAMES)
    assert _splits(rooted) == _splits(unrooted)
    record = ctx.config.stages["phylo"]
    assert record.params["root"] == "outgroup"
    assert record.params["outgroup"] == _OUTGROUP


def test_root_midpoint(workdir: Path, fixed_builder) -> None:
    ctx = _workdir(workdir)
    tree = phylo_run(ctx, PhyloParams(treebuilder="faketree_fixed", root="midpoint"))
    rooted = tree.read_text(encoding="utf-8")
    assert len(_root_children(rooted)) == 2
    assert _splits(rooted) == _splits(_fill(_UNROOTED, _NAMES))
    assert (ctx.tree_dir / TREE_UNROOTED_NWK).is_file()
    assert ctx.config.stages["phylo"].params["root"] == "midpoint"


def test_root_none_leaves_the_builder_output(workdir: Path, fixed_builder) -> None:
    ctx = _workdir(workdir)
    tree = phylo_run(ctx, PhyloParams(treebuilder="faketree_fixed", root="none"))
    assert tree.read_text(encoding="utf-8") == _fill(_UNROOTED, _NAMES)
    assert not (ctx.tree_dir / TREE_UNROOTED_NWK).exists()
    assert ctx.config.stages["phylo"].params["root"] == "none"


def test_no_outgroup_defaults_to_no_rooting(workdir: Path, fixed_builder) -> None:
    ctx = _workdir(workdir)
    phylo_run(ctx, PhyloParams(treebuilder="faketree_fixed", no_outgroup=True))
    assert not (ctx.tree_dir / TREE_UNROOTED_NWK).exists()
    assert ctx.config.stages["phylo"].params["root"] == "none"


def test_no_outgroup_with_midpoint_is_allowed(workdir: Path, fixed_builder) -> None:
    ctx = _workdir(workdir)
    tree = phylo_run(
        ctx, PhyloParams(treebuilder="faketree_fixed", no_outgroup=True, root="midpoint")
    )
    assert len(_root_children(tree.read_text(encoding="utf-8"))) == 2
    assert ctx.config.stages["phylo"].params["root"] == "midpoint"


def test_tree_already_rooted_on_the_outgroup_is_left_alone(
    workdir: Path, fixed_builder, caplog
) -> None:
    fixed_builder.template = _ROOTED
    ctx = _workdir(workdir)
    with caplog.at_level(logging.INFO, logger="test-root"):
        tree = phylo_run(ctx, PhyloParams(treebuilder="faketree_fixed"))
    assert tree.read_text(encoding="utf-8") == _fill(_ROOTED, _NAMES)
    assert not (ctx.tree_dir / TREE_UNROOTED_NWK).exists()
    assert "already rooted on the outgroup" in caplog.text
    assert ctx.config.stages["phylo"].params["root"] == "outgroup"


def test_outgroup_at_a_trifurcating_root_is_rerooted_with_one_info_line(
    workdir: Path, fixed_builder, caplog
) -> None:
    # IQ-TREE -o and RAxML-NG --outgroup put the outgroup at a trifurcating
    # root; the ingroup is not one clade there, so the tree is rerooted.
    fixed_builder.template = _IQTREE_STYLE
    ctx = _workdir(workdir)
    with caplog.at_level(logging.INFO, logger="test-root"):
        tree = phylo_run(ctx, PhyloParams(treebuilder="faketree_fixed"))
    rooted = tree.read_text(encoding="utf-8")
    assert is_rooted_on(read_newick(rooted), _OUTGROUP)
    assert _splits(rooted) == _splits(_fill(_IQTREE_STYLE, _NAMES))
    assert len([r for r in caplog.records if r.getMessage().startswith("Rooted the ")]) == 1


def test_root_outgroup_without_a_staged_outgroup_is_a_usage_error(
    workdir: Path, fixed_builder
) -> None:
    ctx = _workdir(workdir, outgroup=False)
    with pytest.raises(UserInputError, match="needs an outgroup"):
        phylo_run(ctx, PhyloParams(treebuilder="faketree_fixed", root="outgroup"))
    assert fixed_builder.calls == 0
    assert not (ctx.tree_dir / TREE_NWK).exists()


def test_root_outgroup_with_no_outgroup_is_a_usage_error(workdir: Path, fixed_builder) -> None:
    ctx = _workdir(workdir)
    with pytest.raises(UserInputError, match="--no-outgroup"):
        phylo_run(ctx, PhyloParams(treebuilder="faketree_fixed", no_outgroup=True, root="outgroup"))
    assert fixed_builder.calls == 0


def test_rebuild_with_root_none_removes_the_previous_unrooted_copy(
    workdir: Path, fixed_builder
) -> None:
    ctx = _workdir(workdir)
    phylo_run(ctx, PhyloParams(treebuilder="faketree_fixed"))
    assert (ctx.tree_dir / TREE_UNROOTED_NWK).exists()
    phylo_run(ctx, PhyloParams(treebuilder="faketree_fixed", root="none"))
    assert not (ctx.tree_dir / TREE_UNROOTED_NWK).exists()


def test_phylo_build_step_roots_on_the_outgroup(tmp_path: Path, fixed_builder) -> None:
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    for name in _INGROUP:
        (genomes / f"{name}.fasta").write_text(">s\nACGT\n")
    og = tmp_path / "outgroup"
    og.mkdir()
    (og / f"{_OUTGROUP}.fasta").write_text(">o\nACGT\n")
    acc = tmp_path / "outgroup_accession.txt"
    acc.write_text("GCF_000009.1\n")
    out = tmp_path / "out"
    tree = phylo_build(
        PhyloBuildParams(
            genomes_dir=genomes,
            out_dir=out,
            outgroup_dir=og,
            outgroup_accession=acc,
            phylo=PhyloParams(treebuilder="faketree_fixed"),
        ),
        _LOG,
    )
    assert is_rooted_on(read_newick(tree.read_text(encoding="utf-8")), _OUTGROUP)
    assert (out / "tree" / TREE_UNROOTED_NWK).is_file()


# --- tree2tax --------------------------------------------------------------------


def _relations(workdir: Path) -> tuple[str, str]:
    return (
        (workdir / "tree2tax.tsv").read_text(encoding="utf-8"),
        (workdir / "genomes_map.tsv").read_text(encoding="utf-8"),
    )


def test_tree2tax_relations_unchanged_by_rooting_in_phylo(
    tmp_path: Path, fixed_builder, caplog
) -> None:
    """tree2tax on phylo's outgroup-rooted tree equals tree2tax rooting the
    builder's unrooted output itself (the behaviour before --root)."""
    rooted_wd = tmp_path / "rooted"
    ctx = _workdir(rooted_wd)
    phylo_run(ctx, PhyloParams(treebuilder="faketree_fixed"))
    with caplog.at_level(logging.INFO, logger="test-root"):
        tree2tax_run(ctx, Tree2taxParams())
    assert "already rooted on the outgroup" in caplog.text

    plain_wd = tmp_path / "plain"
    ctx2 = _workdir(plain_wd)
    phylo_run(ctx2, PhyloParams(treebuilder="faketree_fixed", root="none"))
    tree2tax_run(ctx2, Tree2taxParams())

    t2t_a, map_a = _relations(rooted_wd)
    t2t_b, map_b = _relations(plain_wd)
    assert sorted(t2t_a.splitlines()) == sorted(t2t_b.splitlines())
    assert sorted(map_a.splitlines()) == sorted(map_b.splitlines())


def _root_children_from_relations(t2t: str, root_name: str) -> set[str]:
    return {
        line.split("\t")[0] for line in t2t.splitlines()[1:] if line.split("\t")[1] == root_name
    }


def test_tree2tax_keeps_a_midpoint_root(workdir: Path, fixed_builder) -> None:
    ctx = _workdir(workdir)
    tree = phylo_run(ctx, PhyloParams(treebuilder="faketree_fixed", root="midpoint"))
    tree2tax_run(ctx, Tree2taxParams(root_name="top"))
    t2t, _ = _relations(workdir)
    parsed = read_newick(tree.read_text(encoding="utf-8"))
    expected = {c.taxon.label for c in parsed.seed_node.child_nodes() if c.is_leaf() and c.taxon}
    children = _root_children_from_relations(t2t, "top")
    assert len(children) == 2
    assert expected <= children


def test_tree2tax_midpoint_without_an_outgroup(workdir: Path, fixed_builder) -> None:
    ctx = _workdir(workdir)
    phylo_run(ctx, PhyloParams(treebuilder="faketree_fixed", no_outgroup=True, root="midpoint"))
    tree2tax_run(ctx, Tree2taxParams())
    t2t, gmap = _relations(workdir)
    assert len(_root_children_from_relations(t2t, "root")) == 2
    assert {line.split("\t")[1] for line in gmap.splitlines()} >= set(_INGROUP)


def test_tree2tax_midpoint_root_with_remove_outgroup(workdir: Path, fixed_builder) -> None:
    ctx = _workdir(workdir)
    phylo_run(ctx, PhyloParams(treebuilder="faketree_fixed", root="midpoint"))
    tree2tax_run(ctx, Tree2taxParams(remove_outgroup=True))
    t2t, _ = _relations(workdir)
    assert _OUTGROUP not in {line.split("\t")[0] for line in t2t.splitlines()}


def test_tree2tax_relations_step_keep_root(tmp_path: Path) -> None:
    # A midpoint-rooted tree whose root does not sit on the outgroup's branch.
    names = [*_INGROUP, _OUTGROUP]
    text = "(({a}:0.1,{b}:0.2)90:0.3,({c}:0.3,({d}:0.4,{og}:0.2)80:0.1)70:0.3);\n"
    tree = tmp_path / "tree.nwk"
    tree.write_text(_fill(text, names))
    og = tmp_path / "outgroup"
    og.mkdir()
    (og / f"{_OUTGROUP}.fasta").write_text(">o\nACGT\n")
    acc = tmp_path / "outgroup_accession.txt"
    acc.write_text("GCF_000009.1\n")

    def children(keep_root: bool) -> set[str]:
        out = tmp_path / ("kept" if keep_root else "rerooted")
        t2t, _ = tree2tax_relations(
            Tree2taxStepParams(
                tree=tree,
                out_dir=out,
                outgroup_dir=og,
                outgroup_accession=acc,
                keep_root=keep_root,
            ),
            _LOG,
        )
        return _root_children_from_relations(t2t.read_text(encoding="utf-8"), "root")

    assert _OUTGROUP not in children(True)
    assert _OUTGROUP in children(False)


# --- CLI -------------------------------------------------------------------------

_runner = CliRunner()


def test_phylo_build_cli_refuses_root_outgroup_with_no_outgroup(tmp_path: Path) -> None:
    result = _runner.invoke(
        app,
        [
            "phylo-build",
            "--genomes-dir",
            str(tmp_path),
            "-o",
            str(tmp_path / "out"),
            "--no-outgroup",
            "--root",
            "outgroup",
        ],
    )
    assert result.exit_code == 2
    assert "--no-outgroup" in result.output + str(result.exception or "")


def test_phylo_cli_refuses_an_unknown_root(tmp_path: Path) -> None:
    (tmp_path / "wd").mkdir()
    result = _runner.invoke(app, ["phylo", "-wd", str(tmp_path / "wd"), "--root", "tip"])
    assert result.exit_code != 0
    assert "--root" in result.output + str(result.exception or "")
