"""Rooting of Newick trees, shared by the phylo and tree2tax stages.

Two methods are offered: on the branch leading to an outgroup leaf, and at the
midpoint of the longest leaf-to-leaf path. Both are carried out with dendropy.

Internal node labels (bootstrap or other support values) describe the split
of the branch above the node, not the node. A reroot reverses the branches on
the path between the old and the new root, so dendropy, which keeps a label on
its node, would attach those labels to the wrong branches. The labels are
therefore recorded by the split they describe before the reroot and written
back by split afterwards. Branch lengths are moved by dendropy itself.
"""

from __future__ import annotations

import dendropy

from ..core.errors import WorkdirError

ROOT_OUTGROUP = "outgroup"
ROOT_MIDPOINT = "midpoint"
ROOT_NONE = "none"
ROOT_METHODS = (ROOT_OUTGROUP, ROOT_MIDPOINT, ROOT_NONE)


def read_newick(text: str) -> dendropy.Tree:
    """Parse one Newick tree, keeping '_' in labels (genome names hold them)."""
    return dendropy.Tree.get(data=text, schema="newick", preserve_underscores=True)


def write_newick(tree: dendropy.Tree) -> str:
    """Newick text of ``tree``: no [&R] tag, underscores left unquoted.

    The rooting is expressed by a bifurcating root, which is how the tree
    builders' output and most readers express it as well. Lengths are written
    with 12 significant digits, so a halved branch does not show float noise.
    """
    return tree.as_string(
        schema="newick",
        suppress_rooting=True,
        unquoted_underscores=True,
        preserve_spaces=True,
        real_value_format_specifier=".12g",
    )


def is_rooted_on(tree: dendropy.Tree, leaf_label: str) -> bool:
    """True when the root is bifurcating and one of its two children is the leaf."""
    children = tree.seed_node.child_nodes()
    if len(children) != 2:
        return False
    return any(
        child.is_leaf() and child.taxon is not None and child.taxon.label == leaf_label
        for child in children
    )


def root_on_outgroup(
    tree: dendropy.Tree, leaf_label: str, *, source: object, hint: str = ""
) -> bool:
    """Root ``tree`` on the branch leading to the leaf ``leaf_label``.

    The branch is split in two halves, so the root is bifurcating with the
    outgroup on one side and the whole ingroup on the other. Rooting at the
    outgroup's parent node instead would leave, on an unrooted (trifurcating)
    tree, a root with three children, and the ingroup would not appear as one
    clade. Returns False, leaving the tree unchanged, when it is already
    rooted that way; raises WorkdirError when the label is not a leaf.
    """
    node = tree.find_node_with_taxon_label(leaf_label)
    if node is None or not node.is_leaf():
        raise WorkdirError(f"Outgroup {leaf_label} is not a leaf of the tree {source}.{hint}")
    if is_rooted_on(tree, leaf_label):
        return False
    labels = _labels_by_split(tree)
    length = node.edge.length
    half = None if length is None else length / 2
    tree.reroot_at_edge(node.edge, length1=half, length2=half, update_bipartitions=False)
    _restore_labels(tree, labels)
    return True


def midpoint_root(tree: dendropy.Tree, *, source: object) -> None:
    """Root ``tree`` at the midpoint of its longest leaf-to-leaf path.

    Needs branch lengths; a tree without them has no midpoint.
    """
    if any(
        node.edge.length is None for node in tree.preorder_node_iter() if node is not tree.seed_node
    ):
        raise WorkdirError(
            f"The tree {source} lacks branch lengths, so it cannot be rooted at the "
            "midpoint. Use --root outgroup or --root none."
        )
    labels = _labels_by_split(tree)
    tree.reroot_at_midpoint(update_bipartitions=False)
    _restore_labels(tree, labels)


def _leaf_bits(tree: dendropy.Tree) -> dict[dendropy.Node, int]:
    """Leaf-set bitmask below each node, with leaves numbered by label."""
    leaves = sorted(tree.leaf_node_iter(), key=lambda n: n.taxon.label if n.taxon else "")
    index = {leaf: i for i, leaf in enumerate(leaves)}
    bits: dict[dendropy.Node, int] = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            bits[node] = 1 << index[node]
        else:
            mask = 0
            for child in node.child_nodes():
                mask |= bits[child]
            bits[node] = mask
    return bits


def _split_key(mask: int, full: int) -> int:
    """The side of a split that does not hold leaf 0, so both sides compare equal."""
    return full ^ mask if mask & 1 else mask


def _labels_by_split(tree: dendropy.Tree) -> dict[int, str]:
    # A bifurcating input rooted away from the outgroup has two root branches
    # for one split, and so possibly two labels for it; the last node visited
    # wins. No current tree builder writes such a tree with differing labels.
    bits = _leaf_bits(tree)
    full = bits[tree.seed_node]
    out: dict[int, str] = {}
    for node in tree.postorder_internal_node_iter():
        if node is tree.seed_node or node.label is None:
            continue
        out[_split_key(bits[node], full)] = node.label
    return out


def _restore_labels(tree: dendropy.Tree, labels: dict[int, str]) -> None:
    bits = _leaf_bits(tree)
    full = bits[tree.seed_node]
    for node in tree.postorder_internal_node_iter():
        if node is tree.seed_node:
            node.label = None
            continue
        node.label = labels.get(_split_key(bits[node], full))
