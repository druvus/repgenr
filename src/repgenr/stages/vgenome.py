"""viral genome stage.

Ports ``vgenome.py``: from the viral metadata tables, select sequences by
taxonomy (genus/species/serotype/custom) and length, write them to ``genomes/``,
and determine an outgroup with mashtree (writing ``outgroup/`` and
``outgroup_accession.txt``) so the bacterial-style derep/phylo/tree2tax stages
work unchanged.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from statistics import mean, median, stdev

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..core.binaries import BinarySpec, check_binaries
from ..core.context import WorkdirContext
from ..core.errors import UserInputError, WorkdirError
from ..core.process import run as run_cmd
from ..viral.entrez import TAXNAMES_ORDERED
from ..viral.ncbi_virus import VirusRecord

_MASHTREE = BinarySpec("mashtree", version_args=("--version",), min_version="1.2")


@dataclass
class VgenomeParams:
    target_genus: str | None = None
    target_species: str | None = None
    target_serotype: str | None = None
    target_custom: str | None = None
    length_all: bool = False
    length_deviation: int = 10
    length_method: str = "median_of_medians"
    length_range: str | None = None
    discard: str | None = None
    no_outgroup: bool = False
    outgroup_candidates_taxid_min_genomes: int = 5
    glance: bool = False
    print_fasta_headers: bool = False
    ignore_duplicates: bool = False
    keep_files: bool = False
    extra: dict = field(default_factory=dict)


def run(ctx: WorkdirContext, params: VgenomeParams) -> int:
    logger = ctx.logger
    download_wd = ctx.workdir / "virus_download_wd"
    fasta = download_wd / "download.fa"
    records_json = download_wd / "virus_records.json"
    if not fasta.exists():
        raise WorkdirError("Viral metadata missing. Run the vmetadata stage first.")
    if records_json.exists():
        return _run_records(ctx, params, download_wd, fasta, records_json, logger)

    base_tsv = download_wd / "metadata_base.tsv"
    ncbi_tsv = download_wd / "metadata_ncbi.tsv"
    if not base_tsv.exists() or not ncbi_tsv.exists():
        raise WorkdirError("Viral metadata missing. Run the vmetadata stage first.")

    targets = _parse_targets(params)
    if not targets:
        raise UserInputError("Select a taxonomy: --target-genus/-species/-serotype/-custom.")

    base = _read_base(base_tsv)
    ncbi = _read_ncbi(ncbi_tsv)

    selected, headers = _select_by_taxonomy(fasta, ncbi, targets, logger)
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

    n_written = _write_genomes(ctx, fasta, kept, params, logger)

    if not params.no_outgroup:
        _determine_outgroup(ctx, fasta, base, kept, length_range, params, logger)

    ctx.config.record_stage(
        "vgenome",
        params={"selected": n_written, "no_outgroup": params.no_outgroup},
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Wrote %d viral genomes", n_written)
    return n_written


# --- NCBI Virus (records-based) path ----------------------------------------

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


def _seq_map(fasta: Path) -> dict[str, SeqRecord]:
    """accession -> SeqRecord, from an NCBI Virus genomic.fna (headers '>acc ...')."""
    return {entry.id: entry for entry in SeqIO.parse(str(fasta), "fasta")}


def _passes_length(rec: VirusRecord, lo: int, hi: int, discard, seqs) -> bool:
    if not (lo <= rec.length <= hi) or rec.accession not in seqs:
        return False
    if discard:
        desc = seqs[rec.accession].description
        if any(tag in desc for tag in discard):
            return False
    return True


def _run_records(ctx, params, download_wd, fasta, records_json, logger) -> int:
    from ..core.contracts import SelectionRow, genome_filename, write_selection
    from ..viral.ncbi_virus import read_records

    targets = _parse_targets(params)
    if not targets:
        raise UserInputError("Select a taxonomy: --target-genus/-species/-serotype/-custom.")
    records = read_records(records_json)
    selected = [r for r in records if _record_matches(r, targets)]
    if not selected:
        raise UserInputError("No sequences matched the taxonomy selection.")

    lo, hi = _length_range_records(selected, params, logger)
    discard = [x.strip() for x in params.discard.split(",")] if params.discard else None
    seqs = _seq_map(fasta)
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
    selection_rows: list[SelectionRow] = []
    for r in kept:
        name = genome_filename(r.family, r.genus, r.species, r.accession)
        rec = seqs[r.accession]
        (genomes_dir / name).write_text(f">{rec.description}\n{rec.seq}\n")
        selection_rows.append(
            SelectionRow(r.accession, r.family, r.genus, r.species, False, name)
        )

    if not params.no_outgroup:
        og = _determine_outgroup_records(ctx, records, kept, (lo, hi), params, seqs, logger)
        if og is not None:
            selection_rows.append(
                SelectionRow(og.accession, og.family, og.genus, og.species, True,
                             genome_filename(og.family, og.genus, og.species, og.accession))
            )

    write_selection(ctx.workdir / "selection.tsv", selection_rows)
    ctx.config.record_stage(
        "vgenome",
        params={"source": "ncbi_virus", "selected": len(kept), "no_outgroup": params.no_outgroup},
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Wrote %d viral genomes (NCBI Virus)", len(kept))
    return len(kept)


def _determine_outgroup_records(ctx, records, kept, length_range, params, seqs, logger):
    """Pick an outgroup from species outside the selection, via mashtree."""
    check_binaries((_MASHTREE,))
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
        return None

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
        return None
    matrix = outgroup_wd / "distance_matrix.tsv"
    run_cmd(
        ["mashtree", "--genomesize", str(int(mid)), "--mindepth", "0",
         "--outmatrix", matrix, *genome_files],
        logger=logger, log_prefix="mashtree", stdout_path=outgroup_wd / "mashtree.dnd",
    )
    label = _select_outgroup_from_matrix(matrix, logger) if matrix.exists() else None
    if label is None:
        logger.warning("Could not assign an outgroup; proceeding without one.")
        return None
    acc = label[2:] if label.startswith(("O_", "S_")) else label
    og = rec_by_acc.get(acc)
    if og is None or acc in kept_acc:
        return None

    from ..core.contracts import genome_filename
    ctx.outgroup_dir.mkdir(parents=True, exist_ok=True)
    name = genome_filename(og.family, og.genus, og.species, og.accession)
    (ctx.outgroup_dir / name).write_text(f">{seqs[acc].description}\n{seqs[acc].seq}\n")
    (ctx.workdir / "outgroup_accession.txt").write_text(og.accession + "\n")
    logger.info("Selected outgroup: %s (%s)", og.accession, og.species)
    if not params.keep_files:
        shutil.rmtree(outgroup_wd, ignore_errors=True)
    return og


def _parse_targets(params: VgenomeParams) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    if params.target_genus:
        out["genus"] = [x.strip().lower() for x in params.target_genus.split(",")]
    if params.target_species:
        out["species"] = [x.strip().lower() for x in params.target_species.split(",")]
    if params.target_serotype:
        out["serotype"] = [x.strip().lower() for x in params.target_serotype.split(",")]
    if params.target_custom:
        out["custom"] = [x.strip().lower() for x in params.target_custom.split(",")]
    return out


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


def _select_by_taxonomy(fasta: Path, ncbi: dict, targets: dict, logger):
    selected: dict[str, dict[str, int]] = {}
    headers: dict[str, dict[str, str]] = {}
    for entry in SeqIO.parse(str(fasta), "fasta"):
        bvbrc_id = entry.description.split("| ")[1].split("]")[0]
        taxid = bvbrc_id.split(".")[0]

        ok = True
        for level, values in targets.items():
            if level == "custom":
                hit = False
                for kv in values:
                    key, val = kv.split(":")
                    if _matches(ncbi, taxid, key, [val]):
                        hit = True
                if not hit:
                    ok = False
            elif not _matches(ncbi, taxid, level, values):
                ok = False
        if not ok:
            continue

        selected.setdefault(taxid, {})[bvbrc_id] = len(entry.seq)
        headers.setdefault(taxid, {})[bvbrc_id] = entry.description
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


def _write_genomes(ctx, fasta, kept, params: VgenomeParams, logger) -> int:
    genomes_dir = ctx.genomes_dir
    if genomes_dir.exists():
        shutil.rmtree(genomes_dir)
    genomes_dir.mkdir(parents=True)

    written = 0
    for entry in SeqIO.parse(str(fasta), "fasta"):
        bvbrc_id = entry.description.split("| ")[1].split("]")[0]
        taxid = bvbrc_id.split(".")[0]
        if taxid not in kept or bvbrc_id not in kept[taxid]:
            continue
        name = entry.description.split()[0]
        target = genomes_dir / f"{name}.fasta"
        if target.exists():
            if not params.ignore_duplicates:
                raise WorkdirError(
                    f"Duplicate sequence id {name}. Use --ignore-duplicates to proceed."
                )
            logger.warning("Duplicate sequence id %s; overwriting", name)
        target.write_text(f">{entry.description}\n{entry.seq}\n")
        written += 1
    return written


def _determine_outgroup(
    ctx, fasta, base, kept, length_range, params: VgenomeParams, logger
) -> None:
    check_binaries((_MASHTREE,))
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
    for entry in SeqIO.parse(str(fasta), "fasta"):
        bvbrc_id = entry.description.split("| ")[1].split("]")[0]
        taxid = bvbrc_id.split(".")[0]
        if written.get(taxid, 0) > max_per_taxid:
            continue
        if not (taxid in kept or taxid in candidates) or taxid not in base:
            continue
        metric = base[taxid]["seq_med"]
        if not (metric * 0.9 <= len(entry.seq) <= metric * 1.1):
            continue
        prefix = "S" if taxid in kept else "O"
        name = entry.description.split()[0]
        (genomes_dir / f"{prefix}_{name}.fasta").write_text(f">{entry.description}\n{entry.seq}\n")
        written[taxid] = written.get(taxid, 0) + 1

    matrix = outgroup_wd / "distance_matrix.tsv"
    genome_files = sorted(genomes_dir.glob("*.fasta"))
    run_cmd(
        ["mashtree", "--genomesize", str(int(mean(length_range))), "--mindepth", "0",
         "--outmatrix", matrix, *genome_files],
        logger=logger, log_prefix="mashtree", stdout_path=outgroup_wd / "mashtree.dnd",
    )
    if not matrix.exists():
        raise WorkdirError("mashtree produced no distance matrix; check the installation.")

    outgroup_id = _select_outgroup_from_matrix(matrix, logger)
    if outgroup_id is None:
        raise WorkdirError("Could not assign an outgroup; specify one manually.")
    if outgroup_id.startswith("O_"):
        outgroup_id = outgroup_id[2:]

    (ctx.workdir / "outgroup_accession.txt").write_text(outgroup_id + "\n")
    ctx.outgroup_dir.mkdir(parents=True, exist_ok=True)
    for entry in SeqIO.parse(str(fasta), "fasta"):
        if entry.description.split()[0] == outgroup_id:
            (ctx.outgroup_dir / f"{outgroup_id}.fasta").write_text(
                f">{entry.description}\n{entry.seq}\n"
            )
            break
    logger.info("Selected outgroup: %s", outgroup_id)
    if not params.keep_files:
        shutil.rmtree(outgroup_wd, ignore_errors=True)


def _select_outgroup_from_matrix(matrix: Path, logger) -> str | None:
    header: list[str] = []
    rows: list[list[str]] = []
    with open(matrix) as fo:
        for enum, line in enumerate(fo):
            parts = line.rstrip("\n").split("\t")
            if enum == 0:
                header = parts
            else:
                rows.append(parts)

    # group statistics for an informational warning
    o_vs_s: dict[str, dict[str, float]] = {}
    s_vs_o_all: list[float] = []
    for idx, seq1 in enumerate(header):
        if idx == 0 or not seq1.startswith("S"):
            continue
        for row in rows:
            seq2 = row[0]
            if not seq2.startswith("O"):
                continue
            val = float(row[idx])
            o_vs_s.setdefault(seq2, {})[seq1] = val
            s_vs_o_all.append(val)
    if not o_vs_s or not s_vs_o_all:
        return None

    threshold = mean(s_vs_o_all) - (stdev(s_vs_o_all) if len(s_vs_o_all) > 1 else 0.0)
    threshold = max(threshold, median(s_vs_o_all))

    # candidates sorted by (lowest max distance, then lowest min distance)
    summary = {
        o: {"min": min(vals.values()), "max": max(vals.values())} for o, vals in o_vs_s.items()
    }
    ordered = sorted(summary.items(), key=lambda x: (-x[1]["max"], x[1]["min"]))
    for candidate, dists in ordered:
        if dists["min"] >= threshold:
            return candidate
    return ordered[0][0] if ordered else None
