"""Gubbins recombination masking of a whole-genome alignment."""

from __future__ import annotations

import logging
import shlex
import shutil
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool, runs_on_host
from ..core.errors import ToolExecutionError, WorkdirError
from ..core.plugins import ToolCapabilities
from .base import Masker, MaskParams

_ALLOWED = set("ACGTNacgtn-")


def sanitise_alignment(source: Path, dest: Path, logger: logging.Logger) -> Path:
    """Copy ``source`` with every character outside ACGTN- replaced by N.

    Gubbins accepts only those symbols; SNP typers can leave IUPAC ambiguity
    codes (e.g. a heterozygous call rendered by ``bcftools consensus``) in a
    whole-genome alignment. Returns ``dest`` (or ``source`` when unchanged).
    """
    replaced = 0
    with open(source, encoding="utf-8") as fi, open(dest, "w", encoding="utf-8") as fo:
        for line in fi:
            if line.startswith(">"):
                fo.write(line)
                continue
            body = line.rstrip("\n")
            fixed = "".join(c if c in _ALLOWED else "N" for c in body)
            replaced += sum(1 for a, b in zip(body, fixed, strict=True) if a != b)
            fo.write(fixed + "\n")
    if replaced:
        logger.warning(
            "Replaced %d ambiguous base(s) with N before Gubbins (only ACGTN- are accepted).",
            replaced,
        )
        return dest
    dest.unlink(missing_ok=True)
    return source


def read_fasta(path: Path) -> dict[str, str]:
    """Record id -> sequence (single-line or wrapped), in file order."""
    records: dict[str, str] = {}
    name: str | None = None
    parts: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        records[name] = "".join(parts)
    return records


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n{seq}\n")


def read_recombination_gff(path: Path) -> dict[str, list[tuple[int, int]]]:
    """Recombinant regions per taxon from Gubbins' GFF3 (1-based, inclusive).

    Each feature line carries ``taxa="  a  b"`` in its attributes; a region
    applies to every taxon listed there.
    """
    regions: dict[str, list[tuple[int, int]]] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            start, end = int(cols[3]), int(cols[4])
            attrs = cols[8]
            key = 'taxa="'
            i = attrs.find(key)
            if i < 0:
                continue
            taxa = attrs[i + len(key) : attrs.find('"', i + len(key))].split()
            for taxon in taxa:
                regions.setdefault(taxon, []).append((start, end))
    return regions


def apply_masks(
    records: dict[str, str], regions: dict[str, list[tuple[int, int]]]
) -> dict[str, str]:
    """Set the recombinant regions to N in the taxa they were predicted for."""
    out: dict[str, str] = {}
    for name, seq in records.items():
        spans = regions.get(name)
        if not spans:
            out[name] = seq
            continue
        chars = list(seq)
        for start, end in spans:
            lo, hi = max(start - 1, 0), min(end, len(chars))
            chars[lo:hi] = ["N"] * (hi - lo)
        out[name] = "".join(chars)
    return out


def polymorphic_sites(records: dict[str, str]) -> dict[str, str]:
    """Columns with two or more distinct bases among ACGT (N and gaps ignored)."""
    names = list(records)
    seqs = [records[n].upper() for n in names]
    length = min(len(s) for s in seqs) if seqs else 0
    keep = [i for i in range(length) if len({s[i] for s in seqs if s[i] in "ACGT"}) >= 2]
    return {n: "".join(s[i] for i in keep) for n, s in zip(names, seqs, strict=True)}


# Gubbins asks for one of these whenever it runs with more than one thread;
# some conda builds (osx-arm64 among them) ship only the single-threaded
# ``raxmlHPC``, and Gubbins then exits without building a tree.
MULTITHREADED_RAXML = (
    "raxmlHPC-PTHREADS-AVX2",
    "raxmlHPC-PTHREADS-AVX",
    "raxmlHPC-PTHREADS-SSE3",
    "raxmlHPC-PTHREADS",
)
IQTREE_BINARIES = ("iqtree2", "iqtree")


def multithreaded_raxml_available() -> bool:
    return any(shutil.which(name) for name in MULTITHREADED_RAXML)


def resolve_tree_builder(
    requested: str | None, threads: int, logger: logging.Logger, *, on_host: bool
) -> tuple[str | None, int]:
    """Pick the Gubbins tree builder and thread count that can actually run.

    A requested builder is passed through unchanged. Otherwise Gubbins' own
    default (RAxML) stands unless the run is native, multi-threaded and the
    host has no multi-threaded RAxML build: then IQ-TREE is used when present,
    else Gubbins runs single-threaded. Both fallbacks are logged.
    """
    if requested:
        return requested, threads
    if threads <= 1 or not on_host or multithreaded_raxml_available():
        return None, threads
    if any(shutil.which(name) for name in IQTREE_BINARIES):
        logger.warning(
            "No multi-threaded RAxML build (%s) on PATH; Gubbins would exit with "
            "--threads %d. Using --tree-builder iqtree instead (set "
            "--tool-arg gubbins_tree_builder=... to choose).",
            "/".join(MULTITHREADED_RAXML),
            threads,
        )
        return "iqtree", threads
    logger.warning(
        "No multi-threaded RAxML build (%s) or IQ-TREE on PATH; running Gubbins "
        "with a single thread.",
        "/".join(MULTITHREADED_RAXML),
    )
    return None, 1


# Above this fraction of variable columns a set is far outside what a
# recombination scanner is built for: Gubbins expects isolates of one species.
DIVERGENCE_WARN_FRACTION = 0.10


def variable_fraction(records: dict[str, str], sample: int = 200_000) -> tuple[float, int]:
    """Fraction of columns with more than one base among ACGT, and the length.

    Estimated on up to ``sample`` evenly spaced columns so a genus-scale
    alignment (tens of genomes, millions of columns) costs a second or two.
    """
    seqs = [s.upper() for s in records.values()]
    length = min((len(s) for s in seqs), default=0)
    if length == 0 or len(seqs) < 2:
        return 0.0, length
    step = max(1, length // sample)
    columns = range(0, length, step)
    varying = sum(1 for i in columns if len({s[i] for s in seqs if s[i] in "ACGT"}) >= 2)
    return varying / len(list(columns)), length


class GubbinsMasker(Masker):
    capabilities = ToolCapabilities(
        name="gubbins",
        required_binaries=(BinarySpec("run_gubbins.py", version_args=("--version",)),),
        conda=("bioconda::gubbins",),
        # gubbins_tree_builder / gubbins_first_tree_builder name Gubbins'
        # --tree-builder / --first-tree-builder; gubbins_args is a quoted
        # string of further run_gubbins.py arguments.
        accepted_extras=frozenset(
            {"gubbins_tree_builder", "gubbins_first_tree_builder", "gubbins_args"}
        ),
    )

    def mask(
        self,
        full_alignment: Path,
        out_dir: Path,
        params: MaskParams,
        logger: logging.Logger,
    ) -> Path:
        """Run Gubbins on the whole-genome alignment; return the
        recombination-filtered polymorphic-sites FASTA."""
        out_dir.mkdir(parents=True, exist_ok=True)
        prefix = out_dir / "gubbins"
        cleaned = sanitise_alignment(full_alignment, out_dir / "input_alignment.fasta", logger)
        records = read_fasta(cleaned)
        excluded = {name for name in records if name in params.exclude}
        if excluded:
            # Gubbins sees the ingroup only; its predictions are applied to
            # every record below, so the excluded outgroup keeps its place.
            subset = out_dir / "ingroup_alignment.fasta"
            write_fasta(subset, {n: s for n, s in records.items() if n not in excluded})
            logger.info(
                "Gubbins runs on %d ingroup genome(s); %s stay(s) out of the scan.",
                len(records) - len(excluded),
                ", ".join(sorted(excluded)),
            )
            gubbins_input = subset
        else:
            gubbins_input = cleaned
        requested = params.extra.get("gubbins_tree_builder")
        tree_builder, threads = resolve_tree_builder(
            str(requested) if requested else None,
            params.threads,
            logger,
            on_host=runs_on_host(self.capabilities),
        )
        scan_records = read_fasta(gubbins_input) if excluded else records
        fraction, length = variable_fraction(scan_records)
        logger.info(
            "Gubbins input: %d genome(s), %d columns, about %.1f%% variable.",
            len(scan_records),
            length,
            100 * fraction,
        )
        if fraction > DIVERGENCE_WARN_FRACTION:
            logger.warning(
                "About %.0f%% of the alignment is variable. Gubbins is built for "
                "isolates of one species, where a few percent is usual; on a set "
                "this diverse it may exhaust its stack or return regions that mean "
                "little. Consider a species-level target or --mask none.",
                100 * fraction,
            )
        argv: list[str | Path] = ["run_gubbins.py", "--threads", str(threads)]
        if tree_builder:
            argv += ["--tree-builder", tree_builder]
        first = params.extra.get("gubbins_first_tree_builder")
        if first:
            argv += ["--first-tree-builder", str(first)]
        argv += shlex.split(str(params.extra.get("gubbins_args", "")))
        argv += ["--prefix", prefix, gubbins_input]
        try:
            run_tool(
                self.capabilities,
                argv,
                logger=logger,
                cwd=out_dir,
                log_prefix="gubbins",
            )
        except ToolExecutionError as exc:
            raise WorkdirError(
                f"Gubbins failed on an alignment of {len(scan_records)} genome(s) "
                f"whose columns are about {100 * fraction:.0f}% variable. Its "
                "recombination scan allocates per-SNP arrays on the thread stack "
                "and dies on very diverse input; a within-species set is what it "
                "expects. Run the stage with --mask none, or narrow the target."
            ) from exc
        if not excluded:
            filtered = Path(str(prefix) + ".filtered_polymorphic_sites.fasta")
            if not filtered.exists():
                raise WorkdirError("Gubbins did not produce a filtered polymorphic sites FASTA")
            return filtered
        gff = Path(str(prefix) + ".recombination_predictions.gff")
        if not gff.exists():
            raise WorkdirError("Gubbins did not produce a recombination predictions GFF")
        regions = read_recombination_gff(gff)
        masked = apply_masks(records, regions)
        n_sites = sum(len(r) for r in regions.values())
        logger.info(
            "Masked %d recombinant region(s) across %d genome(s); extracting variable sites.",
            n_sites,
            len(regions),
        )
        out = Path(str(prefix) + ".masked_polymorphic_sites.fasta")
        write_fasta(out, polymorphic_sites(masked))
        return out
