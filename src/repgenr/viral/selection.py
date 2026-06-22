"""NCBI Virus genome selection (the default viral back-end).

From NCBI Virus records (``virus_records.json``) plus the downloaded
``download.fa``, select sequences by taxonomy and length, write canonical
genomes and a ``selection.tsv``, optionally group a segmented virus's segments
per isolate, and pick an outgroup with mashtree. This is the records-based path
of the vgenome stage, kept separate from the legacy BV-BRC path.
"""

from __future__ import annotations

import logging
import re
import shutil
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path
from statistics import mean, median
from typing import TYPE_CHECKING

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..core.binaries import check_binaries
from ..core.context import WorkdirContext
from ..core.contracts import SelectionRow, genome_filename, write_selection
from ..core.errors import UserInputError
from ..core.process import run as run_cmd
from ._common import MASHTREE, parse_targets, select_outgroup_from_matrix
from .ncbi_virus import VirusRecord, read_records

if TYPE_CHECKING:
    from ..stages.vgenome import VgenomeParams


def run_records(
    ctx: WorkdirContext,
    params: VgenomeParams,
    download_wd: Path,
    fasta: Path,
    records_json: Path,
    logger: logging.Logger,
) -> int:
    targets = parse_targets(params)
    if not targets:
        raise UserInputError("Select a taxonomy: --target-genus/-species/-serotype/-custom.")
    records = read_records(records_json)
    selected = [r for r in records if _record_matches(r, targets)]
    if not selected:
        raise UserInputError("No sequences matched the taxonomy selection.")

    discard = [x.strip() for x in params.discard.split(",")] if params.discard else None
    seqs = _seq_map(fasta)
    lo, hi = 0, 0
    if params.group_segments:
        # Segmented viruses: skip the per-segment length filter (segments differ
        # in length); group whole isolates instead.
        kept = [r for r in selected if r.accession in seqs]
    else:
        lo, hi = _length_range_records(selected, params, logger)
        kept = [r for r in selected if _passes_length(r, lo, hi, discard, seqs)]
    if not kept:
        raise UserInputError("No sequences passed length/discard filtering.")

    if params.print_fasta_headers:
        for r in kept:
            logger.info(seqs[r.accession].description)
    if params.glance:
        logger.info("--glance specified; stopping before writing genomes")
        return 0

    genomes_dir = ctx.genomes_dir
    if genomes_dir.exists():
        shutil.rmtree(genomes_dir)
    genomes_dir.mkdir(parents=True)
    if params.group_segments:
        selection_rows = _write_isolate_groups(genomes_dir, kept, seqs, logger)
    else:
        selection_rows = []
        for r in kept:
            name = genome_filename(r.family, r.genus, r.species, r.accession)
            rec = seqs[r.accession]
            (genomes_dir / name).write_text(f">{rec.description}\n{rec.seq}\n")
            selection_rows.append(
                SelectionRow(r.accession, r.family, r.genus, r.species, False, name)
            )

    tool_versions: dict[str, str] = {}
    if not params.no_outgroup and not params.group_segments:
        og, tool_versions = _determine_outgroup_records(
            ctx, records, kept, (lo, hi), params, seqs, logger
        )
        if og is not None:
            selection_rows.append(
                SelectionRow(og.accession, og.family, og.genus, og.species, True,
                             genome_filename(og.family, og.genus, og.species, og.accession))
            )

    write_selection(ctx.workdir / "selection.tsv", selection_rows)
    n_written = sum(1 for r in selection_rows if not r.is_outgroup)
    ctx.config.record_stage(
        "vgenome",
        params={
            "source": "ncbi_virus", "selected": n_written,
            "group_segments": params.group_segments, "no_outgroup": params.no_outgroup,
        },
        tool_versions=tool_versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Wrote %d viral genomes (NCBI Virus)", n_written)
    return n_written


def _norm(value: str) -> str:
    return value.lower().replace("-", "").replace("_", "").replace(" ", "")


def _record_matches(rec, targets: dict[str, list[str]]) -> bool:
    """All requested taxonomy levels must match the record (normalized compare)."""
    for level, values in targets.items():
        norm_values = [_norm(v) for v in values]
        if level == "genus":
            ok = _norm(rec.genus) in norm_values
        elif level == "species":
            ok = _norm(rec.species) in norm_values
        elif level == "serotype":
            ok = any(_norm(v) in _norm(rec.organism) for v in values)
        elif level == "custom":
            def _custom_ok(kv: str) -> bool:
                key, _, val = kv.partition(":")
                return _norm(str(getattr(rec, key.strip(), ""))) == _norm(val)
            ok = all(_custom_ok(kv) for kv in values if ":" in kv)
        else:
            ok = False
        if not ok:
            return False
    return bool(targets)


def _length_range_records(selected, params: VgenomeParams, logger) -> tuple[int, int]:
    if params.length_range:
        try:
            lo, hi = (int(x) for x in params.length_range.split("-"))
        except ValueError as exc:
            raise UserInputError("--length-range must be 'start-end', e.g. 25000-35000") from exc
        return lo, hi
    by_species: dict[str, list[int]] = {}
    for r in selected:
        by_species.setdefault(r.species, []).append(r.length)
    all_lens = [r.length for r in selected]
    if params.length_all:
        return min(all_lens), max(all_lens)
    medians = [int(median(v)) for v in by_species.values()]
    midpoint = median(medians) if params.length_method == "median_of_medians" else mean(all_lens)
    dev = params.length_deviation / 100
    rng = (int(midpoint * (1 - dev)), int(midpoint * (1 + dev)))
    logger.info("Target length range: %d-%d", rng[0], rng[1])
    return rng


def _seq_map(fasta: Path) -> Mapping[str, SeqRecord]:
    """accession -> SeqRecord, from an NCBI Virus genomic.fna (headers '>acc ...').

    Uses an on-disk offset index (``SeqIO.index``) rather than loading every
    record into a dict: only the kept/outgroup accessions are materialized on
    access, so memory stays bounded for a large genus download.
    """
    return SeqIO.index(str(fasta), "fasta")


def _passes_length(rec: VirusRecord, lo: int, hi: int, discard, seqs) -> bool:
    if not (lo <= rec.length <= hi) or rec.accession not in seqs:
        return False
    if discard:
        desc = seqs[rec.accession].description
        if any(tag in desc for tag in discard):
            return False
    return True


def _isolate_token(isolate: str) -> str:
    """Make an isolate name a single safe accession-like token."""
    token = re.sub(r"[^A-Za-z0-9.-]", "", isolate.replace(" ", "-").replace("_", "-"))
    return f"iso-{token or 'NA'}"


def _write_isolate_groups(genomes_dir, records, seqs, logger):
    """Combine each isolate's segments into one genome; keep singletons as-is.

    Records sharing an ``isolate`` name (segmented viruses) are concatenated in
    descending length order (a deterministic, segment-number-free ordering) into a
    single canonical genome; isolates with one sequence (and records without an
    isolate) are written individually.
    """
    groups: dict[str, list] = {}
    singletons: list = []
    for r in records:
        (groups.setdefault(r.isolate, []) if r.isolate else singletons).append(r)
    for iso, recs in list(groups.items()):
        if iso and len(recs) <= 1:
            singletons.extend(recs)
            del groups[iso]

    rows: list[SelectionRow] = []
    grouped = 0
    for iso, recs in groups.items():
        rep = recs[0]
        acc = _isolate_token(iso)
        name = genome_filename(rep.family, rep.genus, rep.species, acc)
        ordered = sorted(recs, key=lambda r: -r.length)
        seq = "".join(str(seqs[r.accession].seq) for r in ordered)
        (genomes_dir / name).write_text(f">{acc} {iso} ({len(ordered)} segments)\n{seq}\n")
        rows.append(SelectionRow(acc, rep.family, rep.genus, rep.species, False, name))
        grouped += len(recs)
    for r in singletons:
        name = genome_filename(r.family, r.genus, r.species, r.accession)
        (genomes_dir / name).write_text(
            f">{seqs[r.accession].description}\n{seqs[r.accession].seq}\n"
        )
        rows.append(SelectionRow(r.accession, r.family, r.genus, r.species, False, name))
    logger.info(
        "Grouped %d segment sequences into %d isolate genomes (+ %d single-record genomes)",
        grouped, len(groups), len(singletons),
    )
    return rows


def _determine_outgroup_records(ctx, records, kept, length_range, params, seqs, logger):
    """Pick an outgroup from species outside the selection, via mashtree.

    Returns ``(outgroup_record_or_None, tool_versions)`` so the caller can record
    the resolved mashtree version in provenance.
    """
    versions = check_binaries((MASHTREE,))
    kept_species = {r.species for r in kept}
    kept_acc = {r.accession for r in kept}
    cand_by_species: dict[str, list] = {}
    for r in records:
        if r.species in kept_species or r.accession not in seqs:
            continue
        cand_by_species.setdefault(r.species, []).append(r)
    candidates = {
        sp: rs for sp, rs in cand_by_species.items()
        if len(rs) >= params.outgroup_candidates_taxid_min_genomes
    }
    if not candidates:
        logger.warning("No outgroup candidates found; proceeding without an outgroup.")
        return None, versions

    outgroup_wd = ctx.workdir / "virus_outgroup_wd"
    if outgroup_wd.exists():
        shutil.rmtree(outgroup_wd)
    gdir = outgroup_wd / "genomes"
    gdir.mkdir(parents=True)
    lo, hi = length_range
    mid = (lo + hi) / 2
    rec_by_acc: dict[str, VirusRecord] = {}

    def _emit(prefix: str, recs: list, cap: int) -> None:
        written = 0
        for r in recs:
            if written >= cap or not (mid * 0.85 <= r.length <= mid * 1.15):
                continue
            (gdir / f"{prefix}_{r.accession}.fasta").write_text(
                f">{r.accession}\n{seqs[r.accession].seq}\n"
            )
            rec_by_acc[r.accession] = r
            written += 1

    _emit("S", list(kept), 12)
    for rs in candidates.values():
        _emit("O", rs, 3)

    genome_files = sorted(gdir.glob("*.fasta"))
    if len([f for f in genome_files if f.name.startswith("O_")]) == 0:
        logger.warning("No length-compatible outgroup candidates; proceeding without one.")
        return None, versions
    matrix = outgroup_wd / "distance_matrix.tsv"
    run_cmd(
        ["mashtree", "--genomesize", str(int(mid)), "--mindepth", "0",
         "--outmatrix", matrix, *genome_files],
        logger=logger, log_prefix="mashtree", stdout_path=outgroup_wd / "mashtree.dnd",
    )
    label = select_outgroup_from_matrix(matrix, logger) if matrix.exists() else None
    if label is None:
        logger.warning("Could not assign an outgroup; proceeding without one.")
        return None, versions
    acc = label[2:] if label.startswith(("O_", "S_")) else label
    og = rec_by_acc.get(acc)
    if og is None or acc in kept_acc:
        return None, versions

    ctx.outgroup_dir.mkdir(parents=True, exist_ok=True)
    name = genome_filename(og.family, og.genus, og.species, og.accession)
    (ctx.outgroup_dir / name).write_text(f">{seqs[acc].description}\n{seqs[acc].seq}\n")
    (ctx.workdir / "outgroup_accession.txt").write_text(og.accession + "\n")
    logger.info("Selected outgroup: %s (%s)", og.accession, og.species)
    if not params.keep_files:
        shutil.rmtree(outgroup_wd, ignore_errors=True)
    return og, versions
