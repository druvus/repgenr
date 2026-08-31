"""sourmash dereplication adapter (k-mer / MinHash ANI clustering).

Sketches each genome and greedily picks representatives: walk genomes ordered by
connectivity, and absorb any genome whose *estimated ANI* to a chosen
representative is above the threshold. Fast and low-memory, useful as a scalable
alternative to alignment-based ANI.

The threshold is on sourmash's ANI estimate, the same scale skder and galah
use — not on raw sketch similarity. At k=31 a genome pair at 99.6 percent ANI
has Jaccard only ~0.8, so thresholding the raw similarity at an ANI-style 0.99
would leave nearly everything a singleton (verified against skani on synthetic
sets; see docs/scaling-audit.md).

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
import hashlib
import logging
import os
import uuid
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import numpy as np
import numpy.typing as npt

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import MissingBinaryError, ToolExecutionError, WorkdirError
from ..core.plugins import ToolCapabilities, parse_extra_int
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
        accepted_extras=frozenset({"ksize", "scaled", "sketch_cache", "dense_fallback"}),
        required_binaries=(
            BinarySpec("sourmash", version_args=("--version",), min_version="4.0"),
        ),
        default_params={"ksize": 31, "scaled": 1000},
        recommended_max_genomes=None,
        supports_native_scaling=True,
    )

    def dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: DerepParams,
        logger: logging.Logger,
    ) -> DerepResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        ksize = parse_extra_int(params.extra, "ksize", self.capabilities.default_params["ksize"])
        scaled = parse_extra_int(params.extra, "scaled", self.capabilities.default_params["scaled"])
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
                # Availability was pre-probed, so this is a genuine tool failure
                # (bad parameters, OOM, broken image) -- surface it. Silently
                # switching to the dense back-end can change which
                # representatives are picked, so the fallback is opt-in.
                if not _dense_fallback_requested(params):
                    raise
                logger.warning(
                    "sourmash branchwater sparse path failed (%s); "
                    "dense fallback requested, retrying with dense compare",
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
                f"(~{len(genomes) ** 2 * 8 / 1e9:.1f} GB). Install the branchwater plugin "
                "(pip install sourmash_plugin_branchwater) for the sparse path, or use a "
                "tool that scales better at this size (e.g. --tool skder)."
            )
        sig_dir = sketch_cache if sketch_cache is not None else (out_dir / "signatures")
        sig_dir.mkdir(parents=True, exist_ok=True)

        # Reuse is decided per genome, never by counting files: the cache dir may
        # be shared with other chunks (disjoint genome sets), so only this call's
        # genomes may be sketched, matched, and compared.
        sig_by_genome = _find_signatures(sig_dir, genomes)
        missing = [g for g in genomes if g not in sig_by_genome]
        if missing:
            if len(missing) < len(genomes):
                logger.info(
                    "Reusing %d cached sourmash signatures, sketching %d",
                    len(genomes) - len(missing), len(missing),
                )
            fofn = write_fofn(missing, out_dir / "genomes.fofn")
            # The genome paths live inside the fofn, not in argv, so the container
            # backend cannot infer their mounts; declare their directories (using
            # un-resolved abspaths to match write_fofn and the backend's bind logic).
            genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in missing})
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
            sig_by_genome = _find_signatures(sig_dir, genomes)
        else:
            logger.info("Reusing %d cached sourmash signatures", len(genomes))

        unsketched = [g.name for g in genomes if g not in sig_by_genome]
        if unsketched:
            raise WorkdirError(
                f"sourmash produced no signatures under {sig_dir} for: "
                + ", ".join(unsketched[:5])
                + ("..." if len(unsketched) > 5 else "")
            )
        matrix_csv = out_dir / "compare.csv"
        sig_files = sorted(sig_by_genome[g] for g in genomes)
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

        # Convert the Jaccard matrix to ANI estimates so the threshold applies
        # on the same scale as skder/galah (raw Jaccard at k=31 is ~0.8 for a
        # 99.6 pct ANI pair). Converting here, rather than via ``compare
        # --ani``, keeps identical-content pairs at 1.0 even for tiny sketches,
        # where sourmash's estimator refuses and reports 0.
        labels, sim = _read_compare_csv(matrix_csv)
        ani = _jaccard_to_ani_matrix(np.asarray(sim, dtype=float), ksize)
        name_by_label = _match_labels_to_genomes(labels, genomes)
        return _greedy_cluster(labels, ani, name_by_label, threshold)

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
        An ANI of ``threshold`` corresponds to containment ``threshold ** ksize``
        (ANI estimate = containment^(1/k)), and average containment never exceeds
        max containment, so prefiltering at ``threshold ** ksize`` yields a
        superset of the wanted edges; the parser then keeps pairs whose
        containment-derived ANI is >= ``threshold``, matching the dense
        ``compare --ani`` graph.
        """
        threads = str(params.threads)
        if sketch_cache is not None:
            # The cache dir may be shared with other chunks (disjoint genome
            # sets), so the zip is keyed by the genome set + sketch params: a
            # different set can never be mistaken for this one, and concurrent
            # writers target distinct temp files promoted atomically.
            sketch_cache.mkdir(parents=True, exist_ok=True)
            digest = _genome_set_digest(genomes, ksize, scaled)
            sigs_zip = sketch_cache / f"signatures-{digest}.zip"
        else:
            sigs_zip = out_dir / "signatures.zip"

        if sigs_zip.exists():
            logger.info("Reusing cached sourmash %s", sigs_zip.name)
        else:
            # manysketch reads a CSV of (name, genome_filename, protein_filename).
            # The name becomes the signature name, which is what pairwise reports --
            # plain sketch leaves it empty, so the edge list would be unlabelled.
            sketch_csv = out_dir / "manysketch.csv"
            lines = ["name,genome_filename,protein_filename"]
            for g in genomes:
                lines.append(f"{g.stem},{os.path.abspath(g)},")
            sketch_csv.write_text("\n".join(lines) + "\n")

            # Unique temp name beside the final zip (same dir, so the replace is
            # atomic and concurrent writers of the same set cannot collide).
            tmp_zip = sigs_zip.parent / f".{sigs_zip.stem}.{uuid.uuid4().hex}.partial.zip"
            genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in genomes})
            run_tool(self.capabilities,
                [
                    "sourmash", "scripts", "manysketch", sketch_csv,
                    "-o", tmp_zip,
                    "-p", f"dna,k={ksize},scaled={scaled}",
                    "-c", threads,
                ],
                logger=logger,
                log_prefix="sourmash",
                extra_mounts=[*genome_dirs, str(sketch_csv), str(tmp_zip.parent)],
            )
            if not tmp_zip.exists():
                raise WorkdirError(f"sourmash manysketch produced no signatures at {tmp_zip}")
            os.replace(tmp_zip, sigs_zip)

        pairwise_csv = out_dir / "pairwise.csv"
        run_tool(self.capabilities,
            [
                "sourmash", "scripts", "pairwise", sigs_zip,
                "-o", pairwise_csv,
                "-t", f"{threshold ** ksize:g}",
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
        neighbors = _parse_pairwise_csv(pairwise_csv, threshold, set(labels), ksize)
        return _sparse_greedy_cluster(labels, neighbors, name_by_label)


def _find_signatures(sig_dir: Path, genomes: Sequence[Path]) -> dict[Path, Path]:
    """Map each genome to its signature file in ``sig_dir``, by exact name.

    ``sourmash sketch --outdir`` names signatures after the input file
    (``<name>.sig``); some producers use the stem instead, so both are accepted.
    Only exact per-genome matches count -- a shared cache dir may hold
    signatures for other genome sets, which must never satisfy a lookup.
    """
    out: dict[Path, Path] = {}
    for g in genomes:
        for cand in (f"{g.name}.sig", f"{g.name}.sig.gz", f"{g.stem}.sig", f"{g.stem}.sig.gz"):
            p = sig_dir / cand
            if p.exists():
                out[g] = p
                break
    return out


def _genome_set_digest(genomes: Sequence[Path], ksize: int, scaled: int) -> str:
    """Short digest identifying a genome set + sketch params for cache keying."""
    blob = "\n".join(sorted(g.name for g in genomes)) + f"|k={ksize}|scaled={scaled}"
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:12]


def _dense_fallback_requested(params: DerepParams) -> bool:
    """Opt-in to retrying a failed sparse run with the dense back-end."""
    if params.extra.get("dense_fallback"):
        return True
    return os.environ.get("REPGENR_SOURMASH_DENSE_FALLBACK", "") not in ("", "0")


_BRANCHWATER_CACHE: dict[tuple, bool] = {}


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

    key = (caps.name, get_config().cache_key())
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


def _jaccard_to_ani_matrix(sim: npt.NDArray[np.float64], ksize: int) -> npt.NDArray[np.float64]:
    """Elementwise ANI estimate from a Jaccard similarity matrix.

    Jaccard's containment equivalent is 2j/(1+j); the ANI estimate is its
    k-th root (the quantity ``sourmash compare --ani`` reports, within 1e-4
    on real genomes). Zero stays zero.
    """
    containment = 2 * sim / (1 + sim)
    return np.where(sim > 0, containment ** (1.0 / ksize), 0.0)


def _parse_pairwise_csv(
    path: Path, threshold: float, known: set[str], ksize: int
) -> dict[str, set[str]]:
    """Parse a branchwater ``pairwise`` edge list into a symmetric neighbour map.

    Keeps pairs whose containment-derived ANI estimate (containment^(1/k), the
    quantity ``sourmash compare --ani`` reports) is >= ``threshold``, mirroring
    the dense graph. Rows without an ``average_containment`` column fall back to
    the containment equivalent of Jaccard, 2j/(1+j). Self-edges are dropped, and
    only labels in ``known`` are kept, so a stray name cannot introduce a
    phantom node.
    """
    min_containment = threshold ** ksize
    neighbors: dict[str, set[str]] = {}
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.DictReader(fo)
        for row in reader:
            q = row.get("query_name", "")
            m = row.get("match_name", "")
            if q == m or q not in known or m not in known:
                continue
            try:
                raw = row.get("average_containment")
                if raw:
                    containment = float(raw)
                else:
                    jaccard = float(row.get("jaccard", "") or "nan")
                    containment = 2 * jaccard / (1 + jaccard)
            except (ValueError, ZeroDivisionError):
                continue
            if not containment >= min_containment:  # NaN also fails here
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
    with open(path, encoding="utf-8", newline="") as fo:
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
