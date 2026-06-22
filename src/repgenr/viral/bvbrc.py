"""Legacy BV-BRC genome selection (the --source bvbrc back-end).

Selects genomes from the BV-BRC group FASTA plus the Entrez-derived taxonomy
tables (``metadata_base.tsv`` / ``metadata_ncbi.tsv``) written by the vmetadata
stage. Kept for back-compatibility; the default path is NCBI Virus
(:mod:`repgenr.viral.selection`). The BV-BRC header form is
``... | <bvbrc_id>] [<species> | ...``, and the per-taxid length stats come from
the metadata tables rather than from per-record fields.
"""

from __future__ import annotations

import logging
import shutil
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from statistics import mean, median
from typing import TYPE_CHECKING

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..core.binaries import check_binaries
from ..core.context import WorkdirContext
from ..core.errors import UserInputError, WorkdirError
from ..core.process import run as run_cmd
from ._common import MASHTREE, parse_targets, select_outgroup_from_matrix
from .entrez import TAXNAMES_ORDERED

if TYPE_CHECKING:
    from ..stages.vgenome import VgenomeParams


def run_select(
    ctx: WorkdirContext,
    params: VgenomeParams,
    fasta: Path,
    base_tsv: Path,
    ncbi_tsv: Path,
    logger: logging.Logger,
) -> int:
    targets = parse_targets(params)
    if not targets:
        raise UserInputError("Select a taxonomy: --target-genus/-species/-serotype/-custom.")

    base = _read_base(base_tsv)
    ncbi = _read_ncbi(ncbi_tsv)

    # One metadata-only pass over the (large) group FASTA; sequences are fetched
    # lazily for just the records actually written. The old path re-parsed the
    # whole FASTA four times.
    records = _scan_fasta(fasta)

    selected, headers = _select_by_taxonomy(records, ncbi, targets, logger)
    if not selected:
        raise UserInputError("No sequences matched the taxonomy selection.")

    length_range = _determine_length_range(selected, base, params, logger)
    kept = _filter_by_length(selected, headers, length_range, params, logger)
    if not kept:
        raise UserInputError("No sequences passed length/discard filtering.")

    if params.print_fasta_headers:
        for taxid in kept:
            for bvbrc_id in kept[taxid]:
                logger.info(headers[taxid][bvbrc_id])
    if params.glance:
        logger.info("--glance specified; stopping before writing genomes")
        return 0

    sequences = _open_sequences(fasta, params.ignore_duplicates)
    n_written = _write_genomes(ctx, records, sequences, kept, params, logger)

    tool_versions: dict[str, str] = {}
    if not params.no_outgroup:
        tool_versions = _determine_outgroup(
            ctx, records, sequences, base, kept, length_range, params, logger
        )

    ctx.config.record_stage(
        "vgenome",
        params={"selected": n_written, "no_outgroup": params.no_outgroup},
        tool_versions=tool_versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Wrote %d viral genomes", n_written)
    return n_written


def _read_base(path: Path) -> dict[str, dict]:
    base: dict[str, dict] = {}
    with open(path) as fo:
        next(fo)
        for line in fo:
            f = line.rstrip("\n").split("\t")
            taxid = f[0]
            base[taxid] = {
                "num": int(f[2]), "seq_min": int(f[3]), "seq_max": int(f[4]),
                "seq_med": int(f[5]), "seq_mean": int(f[6]),
            }
    return base


def _read_ncbi(path: Path) -> dict[str, list[dict]]:
    ncbi: dict[str, list[dict]] = {}
    n = len(TAXNAMES_ORDERED)
    with open(path) as fo:
        next(fo)
        for line in fo:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 + 2 * n:
                continue
            taxid = f[0]
            names = f[3 : 3 + n]
            taxids = f[3 + n : 3 + 2 * n]
            ncbi[taxid] = [
                {"taxlevelname": TAXNAMES_ORDERED[i], "taxname": names[i], "taxid": taxids[i]}
                for i in range(n)
            ]
    return ncbi


def _matches(ncbi: dict, taxid: str, level: str, values: list[str]) -> bool:
    for entry in ncbi.get(taxid, []):
        if entry["taxlevelname"] == level and (
            entry["taxname"].lower() in values or entry["taxid"] in values
        ):
            return True
    return False


@dataclass
class _Record:
    """Per-record BV-BRC metadata (no sequence held)."""

    name: str  # record id / first header token; also the output filename stem
    bvbrc_id: str
    taxid: str
    description: str
    length: int


def _scan_fasta(fasta: Path) -> list[_Record]:
    """One pass over the group FASTA -> per-record metadata (sequences not held)."""
    records: list[_Record] = []
    for rec in SeqIO.parse(str(fasta), "fasta"):
        bvbrc_id = rec.description.split("| ")[1].split("]")[0]
        records.append(
            _Record(
                name=rec.id,
                bvbrc_id=bvbrc_id,
                taxid=bvbrc_id.split(".")[0],
                description=rec.description,
                length=len(rec.seq),
            )
        )
    return records


def _open_sequences(fasta: Path, ignore_duplicates: bool) -> Mapping[str, SeqRecord]:
    """Lazy name -> SeqRecord index so only the written records load sequences.

    ``SeqIO.index`` rejects duplicate record ids; with --ignore-duplicates fall
    back to an in-memory dict (last record wins), matching the old overwrite
    behavior, otherwise raise a clear error.
    """
    try:
        return SeqIO.index(str(fasta), "fasta")
    except ValueError as exc:
        if not ignore_duplicates:
            raise WorkdirError(
                "Duplicate sequence id in the BV-BRC FASTA. Use --ignore-duplicates."
            ) from exc
        return {rec.id: rec for rec in SeqIO.parse(str(fasta), "fasta")}


def _select_by_taxonomy(records: list[_Record], ncbi: dict, targets: dict, logger):
    selected: dict[str, dict[str, int]] = {}
    headers: dict[str, dict[str, str]] = {}
    for rec in records:
        ok = True
        for level, values in targets.items():
            if level == "custom":
                hit = False
                for kv in values:
                    key, val = kv.split(":")
                    if _matches(ncbi, rec.taxid, key, [val]):
                        hit = True
                if not hit:
                    ok = False
            elif not _matches(ncbi, rec.taxid, level, values):
                ok = False
        if not ok:
            continue

        selected.setdefault(rec.taxid, {})[rec.bvbrc_id] = rec.length
        headers.setdefault(rec.taxid, {})[rec.bvbrc_id] = rec.description
    logger.info("Taxonomy selection matched %d taxids", len(selected))
    return selected, headers


def _determine_length_range(selected, base, params: VgenomeParams, logger) -> tuple[int, int]:
    if params.length_range:
        try:
            lo, hi = (int(x) for x in params.length_range.split("-"))
        except ValueError as exc:
            raise UserInputError("--length-range must be 'start-end', e.g. 25000-35000") from exc
        return lo, hi

    stat_values: list[int] = []
    extremes: list[int] = []
    for taxid in selected:
        if taxid not in base:
            continue
        meta = base[taxid]
        stat_values.append(
            meta["seq_med"] if params.length_method == "median_of_medians" else meta["seq_mean"]
        )
        if params.length_all:
            extremes += [meta["seq_min"], meta["seq_max"]]
    if not stat_values:
        raise UserInputError("Could not determine target length; try --length-range.")

    if params.length_all:
        return min(extremes), max(extremes)
    midpoint = (
        median(stat_values) if params.length_method == "median_of_medians" else mean(stat_values)
    )
    dev = params.length_deviation / 100
    rng = (int(midpoint * (1 - dev)), int(midpoint * (1 + dev)))
    logger.info("Target length range: %d-%d", rng[0], rng[1])
    return rng


def _filter_by_length(selected, headers, length_range, params: VgenomeParams, logger):
    lo, hi = length_range
    discard = [x.strip() for x in params.discard.split(",")] if params.discard else None
    kept: dict[str, dict[str, int]] = {}
    for taxid, by_id in selected.items():
        for bvbrc_id, seq_len in by_id.items():
            if not (lo <= seq_len <= hi):
                continue
            if discard and any(tag in headers[taxid][bvbrc_id] for tag in discard):
                continue
            kept.setdefault(taxid, {})[bvbrc_id] = seq_len
    logger.info("After length/discard filtering: %d taxids", len(kept))
    return kept


def _write_record(seq: SeqRecord, dest: Path) -> None:
    dest.write_text(f">{seq.description}\n{seq.seq}\n")


def _write_genomes(ctx, records, sequences, kept, params: VgenomeParams, logger) -> int:
    genomes_dir = ctx.genomes_dir
    if genomes_dir.exists():
        shutil.rmtree(genomes_dir)
    genomes_dir.mkdir(parents=True)

    written = 0
    for rec in records:
        if rec.taxid not in kept or rec.bvbrc_id not in kept[rec.taxid]:
            continue
        target = genomes_dir / f"{rec.name}.fasta"
        if target.exists():
            if not params.ignore_duplicates:
                raise WorkdirError(
                    f"Duplicate sequence id {rec.name}. Use --ignore-duplicates to proceed."
                )
            logger.warning("Duplicate sequence id %s; overwriting", rec.name)
        _write_record(sequences[rec.name], target)
        written += 1
    return written


def _determine_outgroup(
    ctx, records, sequences, base, kept, length_range, params: VgenomeParams, logger
) -> dict[str, str]:
    """Pick an outgroup via mashtree; return the resolved tool versions."""
    versions = check_binaries((MASHTREE,))
    outgroup_wd = ctx.workdir / "virus_outgroup_wd"
    if outgroup_wd.exists():
        shutil.rmtree(outgroup_wd)
    genomes_dir = outgroup_wd / "genomes"
    genomes_dir.mkdir(parents=True)

    candidates = {
        taxid for taxid, meta in base.items()
        if taxid not in kept and meta["num"] >= params.outgroup_candidates_taxid_min_genomes
    }
    if not candidates:
        raise WorkdirError(
            "No outgroup candidates found. Lower the minimum genomes per candidate "
            "taxid, or pass --no-outgroup."
        )

    max_per_taxid = 3
    written: dict[str, int] = {}
    for rec in records:
        if written.get(rec.taxid, 0) > max_per_taxid:
            continue
        if not (rec.taxid in kept or rec.taxid in candidates) or rec.taxid not in base:
            continue
        metric = base[rec.taxid]["seq_med"]
        if not (metric * 0.9 <= rec.length <= metric * 1.1):
            continue
        prefix = "S" if rec.taxid in kept else "O"
        _write_record(sequences[rec.name], genomes_dir / f"{prefix}_{rec.name}.fasta")
        written[rec.taxid] = written.get(rec.taxid, 0) + 1

    matrix = outgroup_wd / "distance_matrix.tsv"
    genome_files = sorted(genomes_dir.glob("*.fasta"))
    run_cmd(
        ["mashtree", "--genomesize", str(int(mean(length_range))), "--mindepth", "0",
         "--outmatrix", matrix, *genome_files],
        logger=logger, log_prefix="mashtree", stdout_path=outgroup_wd / "mashtree.dnd",
    )
    if not matrix.exists():
        raise WorkdirError("mashtree produced no distance matrix; check the installation.")

    outgroup_id = select_outgroup_from_matrix(matrix, logger)
    if outgroup_id is None:
        raise WorkdirError("Could not assign an outgroup; specify one manually.")
    if outgroup_id.startswith("O_"):
        outgroup_id = outgroup_id[2:]

    (ctx.workdir / "outgroup_accession.txt").write_text(outgroup_id + "\n")
    ctx.outgroup_dir.mkdir(parents=True, exist_ok=True)
    if outgroup_id in sequences:
        _write_record(sequences[outgroup_id], ctx.outgroup_dir / f"{outgroup_id}.fasta")
    logger.info("Selected outgroup: %s", outgroup_id)
    if not params.keep_files:
        shutil.rmtree(outgroup_wd, ignore_errors=True)
    return versions
