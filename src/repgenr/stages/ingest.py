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
from ..core.errors import UserInputError, WorkdirError
from ..core.integrity import refuse_foreign_rows
from ..core.manifest import GenomeRecord, record_from_selection

OUTGROUP_ACCESSION_TXT = "outgroup_accession.txt"


@dataclass
class IngestParams:
    genomes_dir: str
    selection: str | None = None
    outgroup: str | None = None
    copy: bool = False
    # Discard genomes another entry path appended (assemble --append) instead
    # of refusing to overwrite the selection that holds them.
    drop_foreign: bool = False


@dataclass
class _Plan:
    """What ingest will stage, resolved and checked before anything is written."""

    source: Path
    by_name: dict[str, Path]
    ingroup: list[SelectionRow]
    outgroup_row: SelectionRow | None
    outgroup_file: Path | None


def precheck(ctx: WorkdirContext, params: IngestParams) -> None:
    """Refuse an ingest that cannot proceed, before anything is changed.

    The CLI harness calls this before it marks the stage record incomplete, so
    a mistyped directory or a bad genome file does not leave a finished ingest
    (and the stages built on it) looking interrupted.
    """
    if not params.drop_foreign:
        refuse_foreign_rows(ctx, "ingest", drop_foreign=False, logger=ctx.logger)
    _plan(params, logger=None)


def _plan(params: IngestParams, logger) -> _Plan:
    """Resolve and check the genome set; warn through ``logger`` when given."""
    source = Path(params.genomes_dir).expanduser()
    if not source.is_dir():
        raise UserInputError(f"--genomes-dir {source} is not a directory.")
    files = list_fasta(source)
    others = _other_entries(source, files)
    if not files:
        subdirs = [p for p in others if p.is_dir()]
        where = (
            f"; it holds {len(subdirs)} subdirectories, which ingest does not search"
            if subdirs
            else ""
        )
        raise UserInputError(
            f"No genome FASTA files under {source} (expected {', '.join(FASTA_SUFFIXES)}){where}."
        )
    by_name = {f.name: f for f in files}

    if params.selection:
        rows = _rows_from_selection(Path(params.selection), by_name, others)
    else:
        rows = [_row_from_filename(f.name) for f in files]
        skipped = [p.name for p in others if p.is_file()]
        if skipped and logger is not None:
            logger.warning(
                "Skipped %d file(s) under %s without a FASTA suffix (%s): %s",
                len(skipped),
                source,
                ", ".join(FASTA_SUFFIXES),
                _examples(skipped),
            )

    outgroup_row, outgroup_file = _resolve_outgroup(rows, by_name, params.outgroup)
    outgroup_name = outgroup_row.filename if outgroup_row is not None else None
    ingroup = [r for r in rows if not r.is_outgroup and r.filename != outgroup_name]
    _refuse_duplicates([*ingroup, *([outgroup_row] if outgroup_row is not None else [])])

    _refuse_unusable(
        [by_name[r.filename] for r in ingroup]
        + ([outgroup_file] if outgroup_file is not None else [])
    )
    return _Plan(source, by_name, ingroup, outgroup_row, outgroup_file)


def run(ctx: WorkdirContext, params: IngestParams) -> int:
    logger = ctx.logger
    refuse_foreign_rows(ctx, "ingest", drop_foreign=params.drop_foreign, logger=logger)
    plan = _plan(params, logger)
    source, by_name, ingroup = plan.source, plan.by_name, plan.ingroup
    outgroup_row, outgroup_file = plan.outgroup_row, plan.outgroup_file

    ctx.genomes_dir.mkdir(parents=True, exist_ok=True)
    if params.copy:
        logger.info("Copying %d genomes from %s into %s", len(ingroup), source, ctx.genomes_dir)
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
            "drop_foreign": params.drop_foreign,
            "total": len(ingroup),
        },
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return len(ingroup)


def _rows_from_selection(
    path: Path, by_name: dict[str, Path], others: list[Path]
) -> list[SelectionRow]:
    if not path.is_file():
        raise UserInputError(f"--selection {path} is not a file.")
    try:
        rows = read_selection(path)
    except WorkdirError as exc:
        # The file is the user's input here, not workdir state.
        raise UserInputError(f"--selection: {exc}") from exc
    missing = [r.filename for r in rows if r.filename not in by_name]
    if missing:
        unsupported = {p.name for p in others} & set(missing)
        note = (
            f"; {len(unsupported)} of them exist but lack a FASTA suffix "
            f"({', '.join(FASTA_SUFFIXES)})"
            if unsupported
            else ""
        )
        raise UserInputError(
            f"{len(missing)} selection row(s) name a file not present under --genomes-dir "
            f"(e.g. {', '.join(missing[:3])}){note}."
        )
    return rows


def _other_entries(source: Path, files: list[Path]) -> list[Path]:
    """Entries under ``source`` that are not genome FASTA files (dotfiles excluded)."""
    taken = {f.name for f in files}
    return sorted(p for p in source.iterdir() if not p.name.startswith(".") and p.name not in taken)


def _examples(names: list[str], limit: int = 5) -> str:
    more = f" (+{len(names) - limit} more)" if len(names) > limit else ""
    return ", ".join(names[:limit]) + more


def _refuse_unusable(paths: list[Path]) -> None:
    """Fail on empty or unreadable genome files (a dangling link included) before
    anything is staged.

    One ``stat`` per file: opening each file to check its first bytes cost about
    50 s for 1000 genomes on an exFAT disk, and this check also runs on the
    resume path. ``doctor`` checks the content of the staged genomes.
    """
    bad: list[str] = []
    for path in paths:
        try:
            if path.stat().st_size == 0:
                bad.append(path.name)
        except OSError:
            bad.append(path.name)
    if bad:
        raise UserInputError(
            f"{len(bad)} genome file(s) are empty or unreadable: {_examples(bad)}. "
            "Remove them from --genomes-dir (or from --selection) and re-run."
        )


def _refuse_duplicates(rows: list[SelectionRow]) -> None:
    """Fail when two genomes share an accession or a selection names one file twice.

    The manifest is keyed by accession, so a shared accession kept one row of
    several while genomes/ and selection.tsv kept them all. Non-canonical names
    with four or more tokens (``sample_1_run_A.fasta``) are parsed as
    Family_genus_species_ACCESSION and can collide this way.
    """
    for attr, label in (("accession", "an accession"), ("filename", "a filename")):
        seen: dict[str, str] = {}
        clashes: list[str] = []
        for row in rows:
            key = getattr(row, attr)
            if key in seen:
                clashes.append(f"{key} ({seen[key]}, {row.filename})")
            else:
                seen[key] = row.filename
        if clashes:
            shown = "; ".join(clashes[:3]) + (
                f" (+{len(clashes) - 3} more)" if len(clashes) > 3 else ""
            )
            hint = (
                " Give each genome its own accession with --selection, or rename the files "
                "(Family_genus_species_ACCESSION.fasta)."
                if attr == "accession"
                else " List each file once in --selection."
            )
            raise UserInputError(f"{len(clashes)} genome(s) share {label}: {shown}.{hint}")


def _row_from_filename(name: str) -> SelectionRow:
    family, genus, species, accession = parse_genome_filename(name)
    return SelectionRow(accession, family, genus, species, False, name)


def _resolve_outgroup(
    rows: list[SelectionRow], by_name: dict[str, Path], flag: str | None
) -> tuple[SelectionRow | None, Path | None]:
    """The outgroup row and its source file, from the selection or ``--outgroup``.

    ``--outgroup`` names a source genome (filename, stem or accession) or points
    at a FASTA file anywhere. When the selection also marks an outgroup, both
    must name the same genome; a conflict is an error rather than a silent drop
    of the selection's outgroup row.
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
        # A path to a genome under --genomes-dir is that genome's row, so its
        # selection accession is kept rather than one parsed from the filename.
        resolved = candidate.resolve()
        for row in rows:
            if by_name[row.filename].resolve() == resolved:
                _refuse_conflict(from_selection, row, flag)
                return replace(row, is_outgroup=True), by_name[row.filename]
        row = replace(_row_from_filename(candidate.name), is_outgroup=True)
        _refuse_conflict(from_selection, row, flag)
        _refuse_external_clash(rows, row, flag)
        return row, candidate
    for row in rows:
        if flag in (row.filename, strip_fasta_suffix(row.filename), row.accession):
            _refuse_conflict(from_selection, row, flag)
            return replace(row, is_outgroup=True), by_name[row.filename]
    raise UserInputError(
        f"--outgroup {flag!r} is neither a FASTA file nor a genome under --genomes-dir."
    )


def _refuse_conflict(from_selection: list[SelectionRow], chosen: SelectionRow, flag: str) -> None:
    """Fail when ``--outgroup`` names another genome than the selection's outgroup row."""
    if not from_selection:
        return
    marked = from_selection[0]
    if (marked.filename, marked.accession) == (chosen.filename, chosen.accession):
        return
    raise UserInputError(
        f"The selection marks {marked.filename} ({marked.accession}) as the outgroup, "
        f"but --outgroup {flag!r} names {chosen.filename} ({chosen.accession}). "
        "Name the same genome in both, or drop one of them."
    )


def _refuse_external_clash(rows: list[SelectionRow], chosen: SelectionRow, flag: str) -> None:
    """Fail when an outgroup file from outside --genomes-dir shares an ingroup genome's
    filename or accession.

    Both name a genome in genomes/, the manifest and the tree, so the outgroup
    would silently replace that ingroup genome rather than sit beside it.
    """
    for row in rows:
        if row.is_outgroup:
            continue
        if row.filename == chosen.filename or row.accession == chosen.accession:
            raise UserInputError(
                f"--outgroup {flag!r} is a file outside --genomes-dir, but its name gives "
                f"{chosen.filename} ({chosen.accession}), which is also the ingroup genome "
                f"{row.filename} ({row.accession}). Rename the outgroup file, or name the "
                "genome under --genomes-dir to set that genome aside."
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
    return record_from_selection(row, "local")
