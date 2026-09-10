"""ingest stage: populate a working directory from local genome FASTAs.

The offline counterpart of ``metadata`` + ``genome`` (or ``vmetadata`` +
``vgenome``): it stages a directory of FASTA files under ``genomes/``, writes
the same ``selection.tsv`` contract the metadata stage publishes, fills the
SQLite manifest, and optionally sets one genome aside as the outgroup. The
downstream chain (``dereplicate`` -> ``phylo`` -> ``tree2tax``) then runs
without any download, which is what the local Nextflow harness already does
for the data-channel path.

Taxonomy comes from ``--selection`` (an 8-column ``selection.tsv``, which also
carries CheckM-style quality for ``--keeper quality``) or, failing that, from
canonical ``Family_genus_species_ACCESSION.fasta`` filenames; any other name
yields an empty taxonomy and the file stem as accession.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import (
    FASTA_SUFFIXES,
    SELECTION_TSV,
    SelectionRow,
    list_fasta,
    parse_genome_filename,
    read_selection,
    strip_fasta_suffix,
    write_selection,
)
from ..core.errors import UserInputError
from ..core.manifest import GenomeRecord

OUTGROUP_ACCESSION_TXT = "outgroup_accession.txt"


@dataclass
class IngestParams:
    genomes_dir: str
    selection: str | None = None
    outgroup: str | None = None
    copy: bool = False


def run(ctx: WorkdirContext, params: IngestParams) -> int:
    logger = ctx.logger
    source = Path(params.genomes_dir).expanduser()
    if not source.is_dir():
        raise UserInputError(f"--genomes-dir {source} is not a directory.")
    files = list_fasta(source)
    if not files:
        raise UserInputError(
            f"No genome FASTA files under {source} (expected {', '.join(FASTA_SUFFIXES)})."
        )
    by_name = {f.name: f for f in files}

    if params.selection:
        rows = _rows_from_selection(Path(params.selection), by_name)
    else:
        rows = [_row_from_filename(f.name) for f in files]

    outgroup_row, outgroup_file = _resolve_outgroup(rows, by_name, params.outgroup)
    outgroup_name = outgroup_row.filename if outgroup_row is not None else None
    ingroup = [r for r in rows if not r.is_outgroup and r.filename != outgroup_name]

    ctx.genomes_dir.mkdir(parents=True, exist_ok=True)
    _prune(ctx.genomes_dir, {r.filename for r in ingroup}, logger)
    for row in ingroup:
        _stage(by_name[row.filename], ctx.genomes_dir / row.filename, params.copy)

    acc_file = ctx.workdir / OUTGROUP_ACCESSION_TXT
    if outgroup_row is not None and outgroup_file is not None:
        ctx.outgroup_dir.mkdir(parents=True, exist_ok=True)
        _prune(ctx.outgroup_dir, {outgroup_row.filename}, logger)
        _stage(outgroup_file, ctx.outgroup_dir / outgroup_row.filename, params.copy)
        acc_file.write_text(outgroup_row.accession + "\n", encoding="utf-8")
        logger.info("Outgroup: %s (%s)", outgroup_row.accession, outgroup_row.filename)
    else:
        if ctx.outgroup_dir.exists():
            _prune(ctx.outgroup_dir, set(), logger)
        acc_file.unlink(missing_ok=True)

    selection_rows = [*ingroup, *([outgroup_row] if outgroup_row is not None else [])]
    ctx.manifest.replace_genomes([_record(r) for r in selection_rows])
    write_selection(ctx.workdir / SELECTION_TSV, selection_rows)

    mode = "copied" if params.copy else "linked"
    logger.info("Ingested %d genomes from %s (%s)", len(ingroup), source, mode)
    ctx.config.record_stage(
        "ingest",
        params={
            "genomes_dir": str(source),
            "selection": params.selection,
            # The raw flag (doctor re-derives the resume inputs from it) and
            # the accession it resolved to.
            "outgroup": params.outgroup,
            "outgroup_accession": outgroup_row.accession if outgroup_row is not None else None,
            "copy": params.copy,
            "total": len(ingroup),
        },
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return len(ingroup)


def _rows_from_selection(path: Path, by_name: dict[str, Path]) -> list[SelectionRow]:
    if not path.is_file():
        raise UserInputError(f"--selection {path} is not a file.")
    rows = read_selection(path)
    missing = [r.filename for r in rows if r.filename not in by_name]
    if missing:
        raise UserInputError(
            f"{len(missing)} selection row(s) name a file not present under --genomes-dir "
            f"(e.g. {', '.join(missing[:3])})."
        )
    return rows


def _row_from_filename(name: str) -> SelectionRow:
    family, genus, species, accession = parse_genome_filename(name)
    return SelectionRow(accession, family, genus, species, False, name)


def _resolve_outgroup(
    rows: list[SelectionRow], by_name: dict[str, Path], flag: str | None
) -> tuple[SelectionRow | None, Path | None]:
    """The outgroup row and its source file, from the selection or ``--outgroup``.

    ``--outgroup`` names a source genome (filename, stem or accession) or points
    at a FASTA file anywhere; it takes precedence over an outgroup row in the
    selection.
    """
    from_selection = [r for r in rows if r.is_outgroup]
    if len(from_selection) > 1:
        raise UserInputError("The selection marks more than one genome as outgroup.")
    if flag is None:
        if not from_selection:
            return None, None
        return from_selection[0], by_name[from_selection[0].filename]

    candidate = Path(flag).expanduser()
    if candidate.is_file():
        row = replace(_row_from_filename(candidate.name), is_outgroup=True)
        return row, candidate
    for row in rows:
        if flag in (row.filename, strip_fasta_suffix(row.filename), row.accession):
            return replace(row, is_outgroup=True), by_name[row.filename]
    raise UserInputError(
        f"--outgroup {flag!r} is neither a FASTA file nor a genome under --genomes-dir."
    )


def _stage(src: Path, dst: Path, copy: bool) -> None:
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    if copy:
        shutil.copy2(src, dst)
    else:
        # abspath, not resolve(): a resolved macOS firmlink path is not what the
        # container backend bind-mounts (see write_fofn in core.process).
        os.symlink(os.path.abspath(src), dst)


def _prune(directory: Path, keep: set[str], logger) -> None:
    """Remove genome FASTA files (or links) no longer part of the selection."""
    for f in directory.iterdir():
        if f.name in keep or not f.name.endswith(FASTA_SUFFIXES):
            continue
        if f.is_symlink() or f.is_file():
            f.unlink()
            logger.info("Removed %s (no longer selected)", f.name)


def _record(row: SelectionRow) -> GenomeRecord:
    return GenomeRecord(
        accession=row.accession,
        filename=row.filename,
        source="local",
        family=row.family or None,
        genus=row.genus or None,
        species=row.species or None,
        is_outgroup=row.is_outgroup,
        completeness=row.completeness,
        contamination=row.contamination,
    )
