"""sourmash dereplication adapter (k-mer / MinHash ANI clustering).

Sketches each genome and greedily picks representatives: walk genomes ordered by
connectivity, and absorb any genome whose Jaccard similarity to a chosen
representative is above the threshold. Fast and low-memory, useful as a scalable
alternative to alignment-based ANI.

Two back-ends compute the pairwise similarities, picked automatically:

* **Sparse (preferred at scale)**: when the ``sourmash_plugin_branchwater``
  plugin is installed, ``sourmash scripts manysketch`` + ``sourmash scripts
  pairwise`` produce only the above-threshold edges (an edge list), never the
  dense N x N matrix. This keeps memory roughly linear in the number of close
  pairs rather than quadratic in the number of genomes, which matters at 10k+.
* **Dense (fallback)**: plain ``sourmash sketch`` + ``sourmash compare`` build
  the full N x N similarity matrix. Used when the plugin is absent (e.g. inside a
  stock BioContainer). For the same threshold both back-ends pick the same
  representatives on well-separated inputs.
"""

from __future__ import annotations

import csv
import logging
import os
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import numpy as np
import numpy.typing as npt

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import MissingBinaryError, ToolExecutionError, WorkdirError
from ..core.plugins import ToolCapabilities
from ..core.process import write_fofn
from .base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepParams,
    DerepResult,
)

# Above this, the dense N x N float64 matrix is too large to hold in memory
# (~0.2 GB at 5k, ~20 GB at 50k); require the sparse branchwater path instead.
_DENSE_MAX_GENOMES = 5000


class SourmashDereplicator(Dereplicator):
    capabilities = ToolCapabilities(
        name="sourmash",
        conda=("bioconda::sourmash",),
        required_binaries=(
            BinarySpec("sourmash", version_args=("--version",), min_version="4.0"),
        ),
        default_params={"ksize": 31, "scaled": 1000},
        recommended_max_genomes=None,
        supports_native_scaling=True,
        threads_param=None,
    )

    def dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: DerepParams,
        logger: logging.Logger,
    ) -> DerepResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        ksize = int(params.extra.get("ksize", self.capabilities.default_params["ksize"]))
        scaled = int(params.extra.get("scaled", self.capabilities.default_params["scaled"]))
        sani = params.secondary_ani
        threshold = sani if sani <= 1.0 else sani / 100

        # Signatures depend only on the genome set + ksize/scaled, not the ANI
        # threshold. The --target-reps search passes a shared sketch_cache dir so
        # the (expensive) sketching runs once and every threshold iteration
        # reuses it instead of re-sketching the same genomes.
        cache = params.extra.get("sketch_cache")
        sketch_cache = Path(cache) if cache else None

        if _branchwater_available(self.capabilities, logger):
            try:
                clusters, status = self._sparse_dereplicate(
                    genomes, out_dir, ksize, scaled, threshold, params, logger, sketch_cache
                )
            except (ToolExecutionError, WorkdirError) as exc:
                logger.warning(
                    "sourmash branchwater sparse path failed (%s); "
                    "falling back to dense compare",
                    exc,
                )
                clusters, status = self._dense_dereplicate(
                    genomes, out_dir, ksize, scaled, threshold, logger, sketch_cache
                )
        else:
            clusters, status = self._dense_dereplicate(
                genomes, out_dir, ksize, scaled, threshold, logger, sketch_cache
            )

        rep_paths = [p for p in genomes if p.name in clusters]
        return DerepResult(
            representatives=sorted(rep_paths),
            clusters=clusters,
            genome_status=status,
        )

    def _dense_dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        ksize: int,
        scaled: int,
        threshold: float,
        logger: logging.Logger,
        sketch_cache: Path | None = None,
    ) -> tuple[dict[str, list[str]], dict[str, str]]:
        """Stock sourmash sketch + N x N compare (no plugin needed)."""
        if len(genomes) > _DENSE_MAX_GENOMES:
            raise WorkdirError(
                f"sourmash dense compare needs an N x N matrix for {len(genomes)} genomes "
                f"(~{len(genomes) ** 2 * 8 / 1e9:.0f} GB). Install the branchwater plugin "
                "(pip install sourmash_plugin_branchwater) for the sparse path, or use a "
                "tool that scales better at this size (e.g. --tool skder)."
            )
        sig_dir = sketch_cache if sketch_cache is not None else (out_dir / "signatures")
        sig_dir.mkdir(parents=True, exist_ok=True)

        def _sigs() -> list[Path]:
            # Skip macOS AppleDouble companions ("._*") on exFAT/NTFS volumes.
            return [
                p for p in (sorted(sig_dir.glob("*.sig")) + sorted(sig_dir.glob("*.sig.gz")))
                if not p.name.startswith("._")
            ]

        cached = _sigs()
        if len(cached) >= len(genomes):
            logger.info("Reusing %d cached sourmash signatures", len(cached))
        else:
            fofn = write_fofn(genomes, out_dir / "genomes.fofn")
            # The genome paths live inside the fofn, not in argv, so the container
            # backend cannot infer their mounts; declare their directories (using
            # un-resolved abspaths to match write_fofn and the backend's bind logic).
            genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in genomes})
            # one signature file per genome, named by genome basename
            run_tool(self.capabilities,
                [
                    "sourmash", "sketch", "dna",
                    "-p", f"k={ksize},scaled={scaled}",
                    "--from-file", fofn,
                    "--outdir", sig_dir,
                ],
                logger=logger,
                log_prefix="sourmash",
                extra_mounts=genome_dirs,
            )

        matrix_csv = out_dir / "compare.csv"
        sig_files = _sigs()
        if not sig_files:
            raise WorkdirError(f"sourmash produced no signatures under {sig_dir}")
        # Pass signatures via --from-file, never on argv (ARG_MAX at scale).
        compare_fofn = write_fofn(sig_files, out_dir / "signatures.fofn")
        run_tool(self.capabilities,
            [
                "sourmash", "compare",
                "-k", str(ksize),
                "--csv", matrix_csv,
                "--from-file", compare_fofn,
            ],
            logger=logger,
            log_prefix="sourmash",
        )

        labels, sim = _read_compare_csv(matrix_csv)
        name_by_label = _match_labels_to_genomes(labels, genomes)
        return _greedy_cluster(labels, sim, name_by_label, threshold)

    def _sparse_dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        ksize: int,
        scaled: int,
        threshold: float,
        params: DerepParams,
        logger: logging.Logger,
        sketch_cache: Path | None = None,
    ) -> tuple[dict[str, list[str]], dict[str, str]]:
        """Branchwater manysketch + pairwise: emit only above-threshold edges.

        ``pairwise -t`` is a *containment* threshold for which pairs to report.
        Containment is always >= Jaccard, so reporting at the Jaccard ``threshold``
        yields a superset of the edges we want; we then keep only pairs whose
        Jaccard column is >= ``threshold``, matching the dense ``compare`` graph.
        """
        threads = str(params.threads)
        sigs_zip = (
            sketch_cache / "signatures.zip" if sketch_cache is not None
            else out_dir / "signatures.zip"
        )
        if sketch_cache is not None:
            sketch_cache.mkdir(parents=True, exist_ok=True)

        if sigs_zip.exists():
            logger.info("Reusing cached sourmash signatures.zip")
        else:
            # manysketch reads a CSV of (name, genome_filename, protein_filename).
            # The name becomes the signature name, which is what pairwise reports --
            # plain sketch leaves it empty, so the edge list would be unlabelled.
            sketch_csv = out_dir / "manysketch.csv"
            lines = ["name,genome_filename,protein_filename"]
            for g in genomes:
                lines.append(f"{g.stem},{os.path.abspath(g)},")
            sketch_csv.write_text("\n".join(lines) + "\n")

            genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in genomes})
            run_tool(self.capabilities,
                [
                    "sourmash", "scripts", "manysketch", sketch_csv,
                    "-o", sigs_zip,
                    "-p", f"dna,k={ksize},scaled={scaled}",
                    "-c", threads,
                ],
                logger=logger,
                log_prefix="sourmash",
                extra_mounts=[*genome_dirs, str(sketch_csv)],
            )
            if not sigs_zip.exists():
                raise WorkdirError(f"sourmash manysketch produced no signatures at {sigs_zip}")

        pairwise_csv = out_dir / "pairwise.csv"
        run_tool(self.capabilities,
            [
                "sourmash", "scripts", "pairwise", sigs_zip,
                "-o", pairwise_csv,
                "-t", f"{threshold:g}",
                "-k", str(ksize),
                "-c", threads,
            ],
            logger=logger,
            log_prefix="sourmash",
        )

        name_by_label = {g.stem: g.name for g in genomes}
        # Iterate in genome-basename order so the greedy tie-break (which member of
        # a mutually-similar group becomes the representative) matches the dense
        # ``compare`` path, whose label order is the sorted signature-file glob.
        labels = [g.stem for g in sorted(genomes, key=lambda g: g.name)]
        neighbors = _parse_pairwise_csv(pairwise_csv, threshold, set(labels))
        return _sparse_greedy_cluster(labels, neighbors, name_by_label)


_BRANCHWATER_CACHE: dict[tuple[str, int], bool] = {}


def _branchwater_available(caps: ToolCapabilities, logger: logging.Logger) -> bool:
    """True when the branchwater plugin (``sourmash scripts pairwise``) is usable.

    Probed through ``run_tool`` so the answer reflects where the tool actually
    runs: native (plugin pip-installed) vs a container image that may not bundle
    it. Any failure (missing subcommand, missing image, missing binary) means the
    dense fallback is used instead.

    The probe is a real subprocess/container start, so memoize it per (tool,
    container config) for the process lifetime -- chunked dereplication calls this
    once per chunk otherwise, paying a container cold-start each time.
    """
    from ..core.containers import get_config

    key = (caps.name, id(get_config()))
    cached = _BRANCHWATER_CACHE.get(key)
    if cached is not None:
        return cached
    try:
        rc = run_tool(
            caps,
            ["sourmash", "scripts", "pairwise", "--help"],
            logger=logger,
            check=False,
            stdout_path=os.devnull,
            log_prefix="sourmash",
        )
        result = rc == 0
    except (ToolExecutionError, MissingBinaryError, FileNotFoundError, OSError):
        result = False
    _BRANCHWATER_CACHE[key] = result
    return result


def _parse_pairwise_csv(
    path: Path, threshold: float, known: set[str]
) -> dict[str, set[str]]:
    """Parse a branchwater ``pairwise`` edge list into a symmetric neighbour map.

    Keeps pairs whose Jaccard is >= ``threshold`` (mirroring the dense graph) and
    drops self-edges. Only labels in ``known`` are kept, so a stray name cannot
    introduce a phantom node.
    """
    neighbors: dict[str, set[str]] = {}
    with open(path, newline="") as fo:
        reader = csv.DictReader(fo)
        for row in reader:
            q = row.get("query_name", "")
            m = row.get("match_name", "")
            if q == m or q not in known or m not in known:
                continue
            try:
                jaccard = float(row.get("jaccard", "") or "nan")
            except ValueError:
                continue
            if jaccard < threshold:
                continue
            neighbors.setdefault(q, set()).add(m)
            neighbors.setdefault(m, set()).add(q)
    return neighbors


def _sparse_greedy_cluster(
    labels: Sequence[str],
    neighbors: Mapping[str, Iterable[str]],
    name_by_label: Mapping[str, str],
) -> tuple[dict[str, list[str]], dict[str, str]]:
    """Greedy representative pick over a sparse adjacency map.

    Mirrors :func:`_greedy_cluster`: prefer the most-connected genome (ties broken
    by input order) as the representative, then claim its still-free neighbours.
    Genomes with no above-threshold neighbour become their own representatives.
    """
    order = sorted(
        range(len(labels)),
        key=lambda i: (-len(set(neighbors.get(labels[i], ()))), i),
    )
    assigned: dict[str, str] = {}  # label -> representative label
    reps: list[str] = []
    for i in order:
        lab = labels[i]
        if lab in assigned:
            continue
        reps.append(lab)
        assigned[lab] = lab
        for nb in neighbors.get(lab, ()):
            if nb not in assigned:
                assigned[nb] = lab

    clusters: dict[str, list[str]] = {}
    status: dict[str, str] = {}
    for rep in reps:
        clusters[name_by_label[rep]] = []
        status[name_by_label[rep]] = STATUS_REPRESENTATIVE
    for lab in labels:
        rep = assigned[lab]
        if lab == rep:
            continue
        clusters[name_by_label[rep]].append(name_by_label[lab])
        status[name_by_label[lab]] = STATUS_CONTAINED
    return clusters, status


def _read_compare_csv(path: Path) -> tuple[list[str], npt.NDArray[np.float64]]:
    with open(path, newline="") as fo:
        labels = next(csv.reader(fo))
        # Parse the N x N body straight into a contiguous float array. At 1000s of
        # genomes this is far smaller and faster than a Python list-of-lists
        # (a 10k x 10k matrix is ~0.8 GB as float64 vs several GB of Python floats).
        matrix = np.loadtxt(fo, delimiter=",", ndmin=2)
    return labels, matrix


def _match_labels_to_genomes(labels: Sequence[str], genomes: Sequence[Path]) -> dict[str, str]:
    """Map a sourmash column label to a genome basename.

    sourmash labels are signature names (often the file path or basename). Match
    by checking which genome basename the label ends with / contains.
    """
    out: dict[str, str] = {}
    by_name = {g.name: g.name for g in genomes}
    stems = {g.stem: g.name for g in genomes}
    for label in labels:
        base = Path(label).name
        if base in by_name:
            out[label] = base
        elif Path(label).stem in stems:
            out[label] = stems[Path(label).stem]
        else:
            out[label] = base
    return out


def _greedy_cluster(
    labels: list[str],
    sim: npt.NDArray[np.float64] | list[list[float]],
    name_by_label: dict[str, str],
    threshold: float,
) -> tuple[dict[str, list[str]], dict[str, str]]:
    n = len(labels)
    # Boolean adjacency (>= threshold); the per-row vectorized ops below replace
    # the previous pure-Python O(n^2) double loops (much faster at 1000s genomes).
    adj = np.asarray(sim, dtype=float) >= threshold
    np.fill_diagonal(adj, False)
    connectivity = adj.sum(axis=1)
    # prefer well-connected genomes as representatives (stable: -connectivity, idx)
    order = np.lexsort((np.arange(n), -connectivity))

    assigned = np.full(n, -1, dtype=np.int64)  # member idx -> representative idx (-1 = free)
    reps: list[int] = []
    for i in order:
        if assigned[i] != -1:
            continue
        reps.append(int(i))
        assigned[i] = i
        # claim every still-free neighbour of i in one vectorized step
        claim = adj[i] & (assigned == -1)
        assigned[claim] = i

    clusters: dict[str, list[str]] = {}
    status: dict[str, str] = {}
    for rep_idx in reps:
        rep_name = name_by_label[labels[rep_idx]]
        clusters[rep_name] = []
        status[rep_name] = STATUS_REPRESENTATIVE
    for member_idx in range(n):
        rep_idx = int(assigned[member_idx])
        if member_idx == rep_idx:
            continue
        rep_name = name_by_label[labels[rep_idx]]
        member_name = name_by_label[labels[member_idx]]
        clusters[rep_name].append(member_name)
        status[member_name] = STATUS_CONTAINED
    return clusters, status
