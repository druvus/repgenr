"""tree2tax stage: turn a rooted tree into FlexTaxD-compatible relations.

Turns a rooted tree into a FlexTaxD child->parent table (via dendropy): read
``tree/tree.nwk``, root by the outgroup, name internal nodes (a user basename or
a leaf-derived hash), then emit ``tree2tax.tsv`` (child -> parent) and
``genomes_map.tsv`` (accession -> leaf). Optionally lists redundant
(dereplicated) genomes under their representative leaf, read from the derep
``clusters.tsv`` contract.
"""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import dendropy

from ..core.context import WorkdirContext
from ..core.contracts import (
    CLUSTERS_TSV,
    GENOMES_MAP_TSV,
    TREE2TAX_TSV,
    TREE_NWK,
    accession_from_filename,
    read_clusters,
    write_genomes_map,
    write_tree2tax,
)
from ..core.errors import WorkdirError


@dataclass
class Tree2taxParams:
    node_basename: str | None = None
    root_name: str = "root"
    remove_outgroup: bool = False
    all_genomes: bool = False
    include_dereplicated: bool = False


@dataclass
class Tree2taxStepParams:
    """Inputs for the stateless tree2tax step (explicit paths, no workdir)."""

    tree: Path
    out_dir: Path
    clusters: Path | None = None
    outgroup_dir: Path | None = None
    outgroup_accession: Path | None = None
    node_basename: str | None = None
    root_name: str = "root"
    remove_outgroup: bool = False
    include_dereplicated: bool = False


def _emit_relations(
    tree_text: str,
    outgroup_leaf: str | None,
    redundant: dict[str, list[str]],
    *,
    node_basename: str | None,
    root_name: str,
    remove_outgroup: bool,
    out_tree2tax: Path,
    out_map: Path,
    logger: logging.Logger,
) -> tuple[Path, Path]:
    """Build FlexTaxD relations from a tree and write the two output tables.

    Stateless core shared by :func:`run` (workdir-bound) and
    :func:`tree2tax_relations` (data-channel step): it roots, names nodes and
    emits ``tree2tax.tsv`` + ``genomes_map.tsv`` to the given paths.
    """
    # preserve_underscores: genome leaf names contain '_' (Family_Genus_species_Acc)
    # and newick otherwise turns underscores into spaces.
    tree = dendropy.Tree.get(data=tree_text, schema="newick", preserve_underscores=True)

    if outgroup_leaf is not None:
        _set_outgroup(tree, outgroup_leaf, logger)

    leaves_nodes = _name_nodes(tree, node_basename)
    _build_paths(tree, leaves_nodes)

    if remove_outgroup and outgroup_leaf in leaves_nodes:
        del leaves_nodes[outgroup_leaf]
        for nodes in leaves_nodes.values():
            if nodes:
                nodes[-1] = root_name

    write_tree2tax(out_tree2tax, _edges(leaves_nodes))
    write_genomes_map(out_map, _genome_map(leaves_nodes, redundant))
    logger.info("Wrote %s and %s", out_tree2tax.name, out_map.name)
    return out_tree2tax, out_map


def tree2tax_relations(
    params: Tree2taxStepParams, logger: logging.Logger
) -> tuple[Path, Path]:
    """Emit FlexTaxD relations from explicit inputs (stateless; no config)."""
    if not params.tree.exists():
        raise WorkdirError(f"Tree not found: {params.tree}. Run the phylo step first.")
    outgroup_leaf = None
    if params.outgroup_dir is not None and params.outgroup_accession is not None:
        outgroup_leaf = _resolve_outgroup_leaf_from(
            params.outgroup_dir, params.outgroup_accession, logger
        )
    redundant = (
        _load_redundant_from(params.clusters)
        if params.include_dereplicated and params.clusters is not None
        else {}
    )
    params.out_dir.mkdir(parents=True, exist_ok=True)
    return _emit_relations(
        params.tree.read_text().strip(),
        outgroup_leaf,
        redundant,
        node_basename=params.node_basename,
        root_name=params.root_name,
        remove_outgroup=params.remove_outgroup,
        out_tree2tax=params.out_dir / TREE2TAX_TSV,
        out_map=params.out_dir / GENOMES_MAP_TSV,
        logger=logger,
    )


def run(ctx: WorkdirContext, params: Tree2taxParams) -> tuple[Path, Path]:
    logger = ctx.logger
    tree_file = ctx.tree_dir / TREE_NWK
    if not tree_file.exists():
        raise WorkdirError(f"Tree not found: {tree_file}. Run the phylo stage first.")

    outgroup_leaf = _resolve_outgroup_leaf(ctx, logger)
    redundant = _load_redundant(ctx) if params.include_dereplicated else {}

    out_tree2tax, out_map = _emit_relations(
        tree_file.read_text().strip(),
        outgroup_leaf,
        redundant,
        node_basename=params.node_basename,
        root_name=params.root_name,
        remove_outgroup=params.remove_outgroup,
        out_tree2tax=ctx.workdir / TREE2TAX_TSV,
        out_map=ctx.workdir / GENOMES_MAP_TSV,
        logger=logger,
    )

    ctx.config.record_stage(
        "tree2tax",
        params={
            "remove_outgroup": params.remove_outgroup,
            "include_dereplicated": params.include_dereplicated,
            "root_name": params.root_name,
        },
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return out_tree2tax, out_map


def _resolve_outgroup_leaf(ctx: WorkdirContext, logger) -> str | None:
    acc_file = ctx.workdir / "outgroup_accession.txt"
    if not acc_file.exists() or not ctx.outgroup_dir.exists():
        logger.warning("No outgroup available; tree is left unrooted")
        return None
    return _resolve_outgroup_leaf_from(ctx.outgroup_dir, acc_file, logger)


def _resolve_outgroup_leaf_from(
    outgroup_dir: Path, accession_file: Path, logger: logging.Logger
) -> str | None:
    if not accession_file.exists() or not outgroup_dir.exists():
        logger.warning("No outgroup available; tree is left unrooted")
        return None
    accession = accession_file.read_text().strip()
    if not accession:
        logger.warning("No outgroup accession recorded; tree is left unrooted")
        return None
    for f in sorted(outgroup_dir.iterdir()):
        if accession in f.name:
            return f.stem
    logger.warning("Outgroup accession %s not present among tree leaves", accession)
    return None


def _leaf_label(node) -> str:
    return node.taxon.label if node.taxon is not None else ""


def _set_outgroup(tree: dendropy.Tree, leaf_label: str, logger) -> None:
    node = tree.find_node_with_taxon_label(leaf_label)
    if node is None:
        logger.warning("Outgroup leaf %s not in tree; leaving unrooted", leaf_label)
        return
    tree.to_outgroup_position(node, update_bipartitions=False)


def _name_nodes(tree: dendropy.Tree, node_basename: str | None) -> dict[str, list[str]]:
    """Record leaf names; give every internal (non-root) node a stable name.

    Internal nodes are named by a hash of their (sorted) descendant leaf labels,
    or ``<basename><n>`` if a basename is given -- so the tree topology becomes a
    set of named clade nodes. Unlike the old ete3 code this keys on ``is_leaf``,
    so an internal support-value label can no longer be mistaken for a taxon.
    """
    leaves_nodes: dict[str, list[str]] = {}
    counter = 0
    for node in tree.preorder_node_iter():
        if node is tree.seed_node:
            continue
        if node.is_leaf():
            name = _leaf_label(node)
            if name:
                leaves_nodes[name] = []
            continue
        if node_basename:
            node.label = f"{node_basename}{counter}"
            counter += 1
        else:
            leaf_labels = sorted(_leaf_label(lf) for lf in node.leaf_iter())
            node.label = hashlib.md5(" ".join(leaf_labels).encode("utf-8")).hexdigest()
    return leaves_nodes


def _build_paths(tree: dendropy.Tree, leaves_nodes: dict[str, list[str]]) -> None:
    for node in tree.leaf_node_iter():
        name = _leaf_label(node)
        if name not in leaves_nodes:
            continue
        for ancestor in node.ancestor_iter(inclusive=False):
            if ancestor is tree.seed_node:
                ancestor.label = "root"
                leaves_nodes[name].append("root")
            else:
                leaves_nodes[name].append(ancestor.label or "")


def _edges(leaves_nodes: dict[str, list[str]]) -> list[tuple[str, str]]:
    edges: list[tuple[str, str]] = []
    for leaf, nodes in leaves_nodes.items():
        path = [leaf, *nodes]
        for i in range(len(path) - 1):
            edges.append((path[i], path[i + 1]))
    return edges


def _load_redundant(ctx: WorkdirContext) -> dict[str, list[str]]:
    return _load_redundant_from(ctx.derep_dir / CLUSTERS_TSV)


def _load_redundant_from(clusters_file: Path) -> dict[str, list[str]]:
    if not clusters_file.exists():
        return {}
    clusters = read_clusters(clusters_file)
    # map representative leaf-stem -> redundant leaf-stems
    out: dict[str, list[str]] = {}
    for rep, members in clusters.items():
        rep_stem = _stem(rep)
        redundant = [_stem(m) for m in members]
        if redundant:
            out[rep_stem] = redundant
    return out


def _genome_map(
    leaves_nodes: dict[str, list[str]], redundant: dict[str, list[str]]
) -> list[tuple[str, str]]:
    mapping: list[tuple[str, str]] = []
    for leaf in leaves_nodes:
        mapping.append((_accession(leaf), leaf))
        for red in redundant.get(leaf, []):
            mapping.append((_accession(red), leaf))
    return mapping


def _stem(filename: str) -> str:
    for suffix in (".fasta.gz", ".fasta", ".fa", ".fna", ".fas"):
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    return filename


def _accession(leaf: str) -> str:
    return accession_from_filename(leaf) or leaf
