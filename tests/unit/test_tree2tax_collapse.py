"""Collapsing weak splits before tree2tax names nodes.

A node merges into its parent when its branch is shorter than the length
threshold or its support is below the support threshold; the root, the root's
children and leaves never collapse.
"""

from __future__ import annotations

import logging

import dendropy
import pytest

from repgenr.stages import tree2tax as mod

_LOG = logging.getLogger("test.collapse")


def _tree(newick: str) -> dendropy.Tree:
    return dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)


def _children_of_root(tree: dendropy.Tree) -> list[str]:
    out = []
    for child in tree.seed_node.child_nodes():
        out.append(child.taxon.label if child.taxon else f"clade{len(child.leaf_nodes())}")
    return sorted(out)


def _leaf_sets(tree: dendropy.Tree) -> set[frozenset[str]]:
    return {
        frozenset(lf.taxon.label for lf in n.leaf_nodes())
        for n in tree.internal_nodes()
        if n is not tree.seed_node
    }


def test_short_branch_collapses_into_a_multifurcation() -> None:
    # (c,d) hang on a 0.0001 branch under the clade ((c,d),e); the clade itself is kept.
    tree = _tree("(out:0.5,(a:0.1,((c:0.1,d:0.1):0.0001,e:0.1):0.2):0.5);")
    stats = mod._collapse_weak_nodes(tree, support=None, length=0.001, logger=_LOG)
    assert stats.collapsed == 1 and stats.by_length == 1 and stats.by_support == 0
    assert _leaf_sets(tree) == {frozenset("cde"), frozenset("acde")}
    # The collapsed node's branch length is carried down to its children.
    c = tree.find_node_with_taxon_label("c")
    assert c.edge.length == pytest.approx(0.1001)


def test_low_support_collapses_on_a_fraction_scale_tree() -> None:
    tree = _tree("(out:0.5,(a:0.1,((c:0.1,d:0.1)0.42:0.05,e:0.1)0.98:0.2)1.0:0.5);")
    stats = mod._collapse_weak_nodes(tree, support=0.7, length=None, logger=_LOG)
    assert stats.collapsed == 1 and stats.by_support == 1
    assert _leaf_sets(tree) == {frozenset("cde"), frozenset("acde")}


def test_percent_scale_supports_are_normalised(caplog) -> None:
    tree = _tree("(out:0.5,(a:0.1,((c:0.1,d:0.1)42:0.05,e:0.1)98:0.2)100:0.5);")
    with caplog.at_level(logging.INFO, logger=_LOG.name):
        stats = mod._collapse_weak_nodes(tree, support=0.7, length=None, logger=_LOG)
    assert stats.collapsed == 1
    assert any("percent" in r.message for r in caplog.records)


def test_root_children_and_leaves_never_collapse() -> None:
    # Every internal branch is tiny, including the ingroup clade's edge at the root.
    tree = _tree("(out:0.00001,((a:0.1,b:0.1):0.00001,(c:0.1,d:0.1):0.00001):0.00001);")
    before_root = _children_of_root(tree)
    stats = mod._collapse_weak_nodes(tree, support=None, length=0.001, logger=_LOG)
    assert _children_of_root(tree) == before_root == ["clade4", "out"]
    assert stats.collapsed == 2  # the two cherries, not the ingroup clade
    assert len(tree.leaf_nodes()) == 5


def test_either_criterion_suffices() -> None:
    tree = _tree(
        "(out:0.5,(a:0.1,((c:0.1,d:0.1)0.99:0.0001,(e:0.1,f:0.1)0.3:0.2)0.99:0.2)1.0:0.5);"
    )
    stats = mod._collapse_weak_nodes(tree, support=0.7, length=0.001, logger=_LOG)
    assert stats.collapsed == 2 and stats.by_length == 1 and stats.by_support == 1


def test_support_threshold_on_a_tree_without_supports_warns_and_keeps_topology(caplog) -> None:
    tree = _tree("(out:0.5,(a:0.1,((c:0.1,d:0.1):0.05,e:0.1):0.2):0.5);")
    before = _leaf_sets(tree)
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        stats = mod._collapse_weak_nodes(tree, support=0.7, length=None, logger=_LOG)
    assert stats.collapsed == 0
    assert _leaf_sets(tree) == before
    assert any("no support values" in r.message for r in caplog.records)


def test_no_thresholds_is_a_no_op() -> None:
    tree = _tree("(out:0.5,(a:0.1,((c:0.1,d:0.1):0.0001,e:0.1):0.2):0.5);")
    before = _leaf_sets(tree)
    stats = mod._collapse_weak_nodes(tree, support=None, length=None, logger=_LOG)
    assert stats.collapsed == 0 and _leaf_sets(tree) == before
