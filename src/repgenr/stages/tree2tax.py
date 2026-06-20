"""tree2tax stage: turn a rooted tree into FlexTaxD-compatible relations.

Ports the ete3 logic of the old ``tree2tax.py`` to the new contracts: read
``tree/tree.nwk``, root by the outgroup, name internal nodes (a user basename or
a leaf-derived hash), then emit ``tree2tax.tsv`` (child -> parent) and
``genomes_map.tsv`` (accession -> leaf). Optionally lists redundant
(dereplicated) genomes under their representative leaf, read from the derep
``clusters.tsv`` contract.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import ete3

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


def run(ctx: WorkdirContext, params: Tree2taxParams) -> tuple[Path, Path]:
    logger = ctx.logger
    tree_file = ctx.tree_dir / TREE_NWK
    if not tree_file.exists():
        raise WorkdirError(f"Tree not found: {tree_file}. Run the phylo stage first.")
    tree = ete3.Tree(tree_file.read_text().strip())

    outgroup_leaf = _resolve_outgroup_leaf(ctx, logger)
    if outgroup_leaf is not None:
        tree.set_outgroup(outgroup_leaf)

    leaves_nodes = _name_nodes(tree, params.node_basename)
    _build_paths(tree, leaves_nodes)

    if params.remove_outgroup and outgroup_leaf in leaves_nodes:
        del leaves_nodes[outgroup_leaf]
        for nodes in leaves_nodes.values():
            if nodes:
                nodes[-1] = params.root_name

    redundant = _load_redundant(ctx) if params.include_dereplicated else {}

    edges = _edges(leaves_nodes)
    out_tree2tax = ctx.workdir / TREE2TAX_TSV
    write_tree2tax(out_tree2tax, edges)

    mapping = _genome_map(leaves_nodes, redundant)
    out_map = ctx.workdir / GENOMES_MAP_TSV
    write_genomes_map(out_map, mapping)

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
    logger.info("Wrote %s and %s", out_tree2tax.name, out_map.name)
    return out_tree2tax, out_map


def _resolve_outgroup_leaf(ctx: WorkdirContext, logger) -> str | None:
    acc_file = ctx.workdir / "outgroup_accession.txt"
    if not acc_file.exists() or not ctx.outgroup_dir.exists():
        logger.warning("No outgroup available; tree is left unrooted")
        return None
    accession = acc_file.read_text().strip()
    for f in ctx.outgroup_dir.iterdir():
        if accession in f.name:
            return f.stem
    logger.warning("Outgroup accession %s not present among tree leaves", accession)
    return None


def _name_nodes(tree: ete3.Tree, node_basename: str | None) -> dict[str, list[str]]:
    leaves_nodes: dict[str, list[str]] = {}
    counter = 0
    for node in tree.iter_descendants():
        if not node.name:
            if node_basename:
                node.name = f"{node_basename}{counter}"
                counter += 1
            else:
                leaf_names = [leaf.name for leaf in node.get_leaves()]
                node.name = hashlib.md5(" ".join(leaf_names).encode("utf-8")).hexdigest()
        else:
            leaves_nodes[node.name] = []
    return leaves_nodes


def _build_paths(tree: ete3.Tree, leaves_nodes: dict[str, list[str]]) -> None:
    for node in tree.iter_descendants():
        if node.name in leaves_nodes:
            for ancestor in node.get_ancestors():
                if ancestor.is_root():
                    ancestor.name = "root"
                leaves_nodes[node.name].append(ancestor.name)


def _edges(leaves_nodes: dict[str, list[str]]) -> list[tuple[str, str]]:
    edges: list[tuple[str, str]] = []
    for leaf, nodes in leaves_nodes.items():
        path = [leaf, *nodes]
        for i in range(len(path) - 1):
            edges.append((path[i], path[i + 1]))
    return edges


def _load_redundant(ctx: WorkdirContext) -> dict[str, list[str]]:
    clusters_file = ctx.derep_dir / CLUSTERS_TSV
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
