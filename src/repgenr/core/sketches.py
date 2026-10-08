"""Genome sketches: one sourmash signature collection per genome of a workdir.

The genome-writing stages (genome, vgenome, ingest, assemble) and the
``repgenr sketch`` command leave ``sketches/<name>.sig.zip`` beside the genome
set, where ``<name>`` is the genome's record name (:func:`record_name`, the
name the aligners and tree builders give it). Each file is a sourmash zip
collection of three DNA signatures, k=21, k=31 and k=51 at scaled=1000 without
abundance, named after the genome, so a consumer can map a signature to its
genome without parsing paths.

Whether a sketch is current is decided from the sha256 of the genome FASTA and
the parameter string recorded with it, not from file times: the manifest holds
both for each genome row (``sketch_file``, ``sketch_params``,
``sketch_digest``); a genome under ``outgroup/`` without a manifest row is
recorded in ``sketches/outgroup.json`` instead.

Sketches are not stage deliverables. A genome set without them is complete;
the stages that compare genomes build their own sketches when these are absent.
Every sketch is written to a hidden temporary file and renamed into place, so
an interrupted run leaves no partial sketch, and the records are written in
batches as the sketches finish, so a rerun sketches only what is missing.
"""

from __future__ import annotations

import json
import logging
import os
import sqlite3
from collections.abc import Collection, Mapping, Sequence
from concurrent.futures import CancelledError, Future, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from . import process
from .containers import run_tool
from .contracts import (
    SELECTION_TSV,
    accession_from_filename,
    atomic_replace,
    list_fasta,
    read_selection,
    record_name,
)
from .errors import MissingBinaryError, RepGenRError, ToolExecutionError, WorkdirError
from .inputs import file_digest
from .manifest import SketchRecord
from .plugins import preflight, tool_available
from .resources import usable_cpus
from .sourmash import SOURMASH_TOOL

if TYPE_CHECKING:
    from .context import WorkdirContext

KSIZES: tuple[int, ...] = (21, 31, 51)
SCALED = 1000
# The parameter string passed to ``sourmash sketch dna -p`` and recorded with
# each sketch; a sketch recorded with another string is stale.
SKETCH_PARAMS = ",".join(f"k={k}" for k in KSIZES) + f",scaled={SCALED}"
SKETCHES_DIR = "sketches"
SKETCH_SUFFIX = ".sig.zip"
# Sketch records of outgroup genomes that have no manifest row.
OUTGROUP_JSON = "outgroup.json"
# sourmash picks the output format from the file name, so the temporary file
# keeps the .sig.zip suffix; the leading dot hides it from every listing.
_PARTIAL = ".partial"
# The threads a writer stage without -t/--threads gives the sketch step.
_DEFAULT_THREADS = 16

ACTION_PRESENT = "present"
ACTION_WRITTEN = "written"
ACTION_REPLACED = "replaced"
ACTION_FORCED = "forced"
ACTION_COPIED = "copied"


@dataclass(frozen=True)
class SketchTarget:
    """A genome to sketch: its FASTA and where its sketch is recorded.

    ``accession`` names its manifest row; None records it in
    ``sketches/outgroup.json`` under the file name.
    """

    path: Path
    accession: str | None

    @property
    def name(self) -> str:
        return record_name(self.path)


@dataclass(frozen=True)
class SketchSource:
    """A sketch another working directory holds, reusable when its record matches."""

    sketch: Path
    record: SketchRecord


@dataclass
class SketchSummary:
    present: int = 0  # up to date, left as they were
    written: int = 0  # sketched for the first time
    replaced: int = 0  # a stale sketch replaced
    forced: int = 0  # a current sketch written again under --force
    copied: int = 0  # reused from another working directory
    removed: int = 0  # sketches of genomes no longer in the set

    def as_dict(self) -> dict[str, int]:
        return {
            "present": self.present,
            "written": self.written,
            "replaced": self.replaced,
            "forced": self.forced,
            "copied": self.copied,
            "removed": self.removed,
        }

    def line(self) -> str:
        forced = f", {self.forced} replaced (--force)" if self.forced else ""
        copied = f", {self.copied} copied" if self.copied else ""
        return (
            f"Sketches: {self.present} present, {self.written} written, "
            f"{self.replaced} stale replaced{forced}{copied}, {self.removed} removed"
        )


@dataclass
class SketchStatus:
    """Record names of the genomes by the state of their sketch."""

    present: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    stale: list[str] = field(default_factory=list)

    @property
    def total(self) -> int:
        return len(self.present) + len(self.missing) + len(self.stale)


def sketches_dir(workdir: Path) -> Path:
    return workdir / SKETCHES_DIR


def sketch_path(ctx: WorkdirContext, genome_path: Path) -> Path:
    """Where the sketch of ``genome_path`` lives in the workdir of ``ctx``."""
    return sketches_dir(ctx.workdir) / (record_name(genome_path) + SKETCH_SUFFIX)


def default_threads() -> int:
    """Sketch threads for a writer stage that takes no -t/--threads."""
    return max(1, min(_DEFAULT_THREADS, usable_cpus()))


def _relative(ctx: WorkdirContext, path: Path) -> str:
    return path.relative_to(ctx.workdir).as_posix()


def _row_paths(ctx: WorkdirContext) -> dict[str, Path]:
    """accession -> genome path of each manifest row, whether or not the file exists.

    The bacterial genome stage records no filename on its outgroup row; that
    row is matched to the FASTA under ``outgroup/`` that carries its accession.
    """
    outgroup_files = list_fasta(ctx.outgroup_dir)
    paths: dict[str, Path] = {}
    for row in ctx.manifest.all_genomes(include_outgroup=True):
        if row.filename:
            directory = ctx.outgroup_dir if row.is_outgroup else ctx.genomes_dir
            paths[row.accession] = directory / row.filename
        elif row.is_outgroup:
            match = [p for p in outgroup_files if accession_from_filename(p.name) == row.accession]
            if len(match) == 1:
                paths[row.accession] = match[0]
    return paths


def sketch_targets(ctx: WorkdirContext) -> list[SketchTarget]:
    """Every genome of the workdir whose FASTA is present, outgroup included.

    The genomes are the manifest rows with a file; a FASTA under
    ``outgroup/`` without a row is added with no accession.
    """
    targets = [
        SketchTarget(path, accession)
        for accession, path in _row_paths(ctx).items()
        if path.is_file()
    ]
    seen = {t.path for t in targets}
    targets += [SketchTarget(p, None) for p in list_fasta(ctx.outgroup_dir) if p not in seen]
    return sorted(targets, key=lambda t: (t.accession is None, t.path.name))


def _expected_names(ctx: WorkdirContext) -> set[str]:
    """Record names of every genome of the set, whether or not its file is present."""
    names = {record_name(path) for path in _row_paths(ctx).values()}
    names |= {record_name(p) for p in list_fasta(ctx.outgroup_dir)}
    return names


def _read_outgroup_json(directory: Path) -> dict[str, SketchRecord]:
    path = directory / OUTGROUP_JSON
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}  # an unreadable record: the outgroup is sketched again
    if not isinstance(data, dict):
        return {}
    out: dict[str, SketchRecord] = {}
    for name, rec in data.items():
        if isinstance(rec, dict):
            out[name] = SketchRecord(rec.get("file"), rec.get("params"), rec.get("digest"))
    return out


def _write_outgroup_json(directory: Path, records: Mapping[str, SketchRecord]) -> None:
    path = directory / OUTGROUP_JSON
    if not records:
        path.unlink(missing_ok=True)
        return
    payload = {
        name: {"file": rec.file, "params": rec.params, "digest": rec.digest}
        for name, rec in sorted(records.items())
    }
    with atomic_replace(path) as fo:
        json.dump(payload, fo, indent=2, sort_keys=True)
        fo.write("\n")


def _record_key(target: SketchTarget) -> str:
    return target.accession if target.accession is not None else target.path.name


def _recorded(ctx: WorkdirContext) -> dict[str, SketchRecord]:
    """Recorded sketches by manifest accession, and by file name for outgroup.json."""
    records = dict(ctx.manifest.sketch_records())
    records.update(_read_outgroup_json(sketches_dir(ctx.workdir)))
    return records


def _is_current(rec: SketchRecord | None, rel: str, digest: str, out: Path) -> bool:
    return (
        rec is not None
        and rec.file == rel
        and rec.params == SKETCH_PARAMS
        and rec.digest == digest
        and out.is_file()
    )


def _partial(out: Path) -> Path:
    stem = out.name[: -len(SKETCH_SUFFIX)]
    return out.with_name(f".{stem}{_PARTIAL}{SKETCH_SUFFIX}")


def _clear_partials(directory: Path) -> None:
    """Remove temporary sketches an interrupted run left behind."""
    if not directory.is_dir():
        return
    for entry in directory.iterdir():
        if entry.name.startswith(".") and entry.name.endswith(_PARTIAL + SKETCH_SUFFIX):
            entry.unlink(missing_ok=True)


def sketch_command(genome: Path, name: str, out: Path) -> list[str | os.PathLike[str]]:
    """The one sourmash call that sketches a genome (three k-mer sizes, one file)."""
    return ["sourmash", "sketch", "dna", "-p", SKETCH_PARAMS, "--name", name, "-o", out, genome]


def _sketch_one(target: SketchTarget, out: Path, logger: logging.Logger) -> None:
    tmp = _partial(out)
    tmp.unlink(missing_ok=True)
    try:
        run_tool(
            SOURMASH_TOOL,
            sketch_command(target.path, target.name, tmp),
            logger=logger,
            log_prefix="sourmash",
            # A linked genome (ingest) lives elsewhere; its directory must be
            # visible inside a container too.
            extra_mounts=[str(target.path.resolve().parent)],
        )
        if not tmp.is_file():
            raise WorkdirError(f"sourmash wrote no sketch for {target.path.name}")
        os.replace(tmp, out)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def _copy_one(source: Path, out: Path) -> None:
    tmp = _partial(out)
    tmp.unlink(missing_ok=True)
    try:
        process.link_or_copy(source, tmp)
        os.replace(tmp, out)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def sketch_genomes(
    ctx: WorkdirContext,
    genomes: Sequence[SketchTarget],
    threads: int,
    logger: logging.Logger,
    *,
    force: bool = False,
    reuse: Mapping[str, SketchSource] | None = None,
) -> SketchSummary:
    """Sketch the genomes whose sketch is missing or stale, ``threads`` at a time.

    One sourmash process per genome. An up-to-date sketch (same FASTA sha256,
    same parameters) is left as it is unless ``force``. ``reuse`` maps record
    names to sketches of another workdir; one whose recorded digest and
    parameters match the genome is copied (or linked) instead of sketched.
    Each sketch is recorded on this thread as soon as it is in place, so an
    interrupted run keeps the records of every sketch it finished. After a
    failure no further genome is started; the running ones finish and are
    recorded, then the first failure is raised.
    """
    directory = sketches_dir(ctx.workdir)
    directory.mkdir(parents=True, exist_ok=True)
    _clear_partials(directory)
    records = _recorded(ctx)
    outgroup_records = _read_outgroup_json(directory)
    reuse = reuse or {}
    summary = SketchSummary()

    def work(target: SketchTarget) -> tuple[SketchTarget, str, SketchRecord]:
        out = sketch_path(ctx, target.path)
        rel = _relative(ctx, out)
        digest = file_digest(target.path)
        fresh = SketchRecord(rel, SKETCH_PARAMS, digest)
        current = _is_current(records.get(_record_key(target)), rel, digest, out)
        if current and not force:
            return target, ACTION_PRESENT, fresh
        existed = out.is_file()
        source = None if force else reuse.get(target.name)
        if (
            source is not None
            and source.record.params == SKETCH_PARAMS
            and source.record.digest == digest
            and source.sketch.is_file()
        ):
            _copy_one(source.sketch, out)
            return target, ACTION_COPIED, fresh
        _sketch_one(target, out, logger)
        if current:
            return target, ACTION_FORCED, fresh
        return target, (ACTION_REPLACED if existed else ACTION_WRITTEN), fresh

    def record(target: SketchTarget, action: str, rec: SketchRecord) -> None:
        setattr(summary, action, getattr(summary, action) + 1)
        if action == ACTION_PRESENT and records.get(_record_key(target)) == rec:
            return
        if target.accession is not None:
            ctx.manifest.set_sketches([(target.accession, rec)])
        else:
            # Records of genomes that left the set are pruned by remove_stale.
            outgroup_records[target.path.name] = rec
            _write_outgroup_json(directory, outgroup_records)

    workers = max(1, min(threads, len(genomes)))
    failures: list[Exception] = []
    pool = ThreadPoolExecutor(max_workers=workers)
    try:
        futures: list[Future[tuple[SketchTarget, str, SketchRecord]]] = [
            pool.submit(work, t) for t in genomes
        ]
        for future in as_completed(futures):
            try:
                target, action, rec = future.result()
            except CancelledError:
                continue  # not started after an earlier failure
            except Exception as exc:
                if not failures:
                    for other in futures:
                        other.cancel()  # queued genomes do not start
                failures.append(exc)
                continue
            record(target, action, rec)
    except BaseException:
        # A termination signal: queued genomes must not start; the running
        # ones settle (their tools are stopped by the handler).
        pool.shutdown(wait=True, cancel_futures=True)
        raise
    pool.shutdown(wait=True)
    if failures:
        raise failures[0]
    return summary


def remove_stale(ctx: WorkdirContext, logger: logging.Logger | None = None) -> int:
    """Remove sketches of genomes no longer in the set; return how many.

    Also removes temporary files of an interrupted run and outgroup.json
    records of genomes that left the set or gained a manifest row, and the
    directory itself when nothing is left in it.
    """
    directory = sketches_dir(ctx.workdir)
    if not directory.is_dir():
        return 0
    _clear_partials(directory)
    expected = _expected_names(ctx)
    removed = 0
    for entry in directory.iterdir():
        if entry.name.startswith(".") or not entry.name.endswith(SKETCH_SUFFIX):
            continue
        if entry.is_dir() or entry.name[: -len(SKETCH_SUFFIX)] in expected:
            continue
        entry.unlink()
        removed += 1
        if logger is not None:
            logger.debug("Removed sketch %s (genome no longer in the set)", entry.name)
    unrowed = {p.name for p in list_fasta(ctx.outgroup_dir)} - {
        path.name for path in _row_paths(ctx).values()
    }
    kept = {n: r for n, r in _read_outgroup_json(directory).items() if n in unrowed}
    _write_outgroup_json(directory, kept)
    if not any(not e.name.startswith(".") for e in directory.iterdir()):
        for entry in directory.iterdir():  # dotfiles only (exFAT, Finder)
            if not entry.is_dir():
                entry.unlink(missing_ok=True)
        if not any(directory.iterdir()):
            directory.rmdir()
    if removed and logger is not None:
        logger.info("Removed %d sketch(es) of genomes no longer in the set", removed)
    return removed


def sketch_status(ctx: WorkdirContext, *, verify: bool = True) -> SketchStatus:
    """Present, missing and stale sketches of the genomes of the workdir.

    With ``verify`` each genome FASTA is hashed and compared with its record;
    without it a sketch counts as present when its file exists and its
    record carries the current parameters (no genome is read).
    """
    records = _recorded(ctx)
    status = SketchStatus()
    for target in sketch_targets(ctx):
        out = sketch_path(ctx, target.path)
        rel = _relative(ctx, out)
        rec = records.get(_record_key(target))
        if not out.is_file():
            status.missing.append(target.name)
        elif verify:
            current = _is_current(rec, rel, file_digest(target.path), out)
            (status.present if current else status.stale).append(target.name)
        elif rec is not None and rec.file == rel and rec.params == SKETCH_PARAMS:
            status.present.append(target.name)
        else:
            status.stale.append(target.name)
    return status


def sketch_counts(workdir: Path) -> tuple[int, int] | None:
    """(sketches present, genomes present) from selection.tsv; None without sketches/.

    Reads no genome and opens no manifest, for ``repgenr status``.
    """
    directory = sketches_dir(workdir)
    if not directory.is_dir():
        return None
    selection = workdir / SELECTION_TSV
    try:
        rows = read_selection(selection) if selection.is_file() else []
    except (OSError, ValueError, RepGenRError):
        return None
    genomes = [
        workdir / ("outgroup" if row.is_outgroup else "genomes") / row.filename for row in rows
    ]
    present = [g for g in genomes if g.is_file()]
    sketched = sum(1 for g in present if (directory / (record_name(g) + SKETCH_SUFFIX)).is_file())
    return sketched, len(present)


def expected_sketch_files(workdir: Path) -> list[Path]:
    """The sketch of every genome selection.tsv lists and the workdir holds."""
    selection = workdir / SELECTION_TSV
    if not selection.is_file():
        return []
    try:
        rows = read_selection(selection)
    except (OSError, ValueError, RepGenRError):
        return []
    out = []
    for row in rows:
        genome = workdir / ("outgroup" if row.is_outgroup else "genomes") / row.filename
        if genome.is_file():
            out.append(sketches_dir(workdir) / (record_name(genome) + SKETCH_SUFFIX))
    return out


def require_sourmash_if_requested(flag: bool | None) -> None:
    """Refuse an explicit --sketch without sourmash before a stage writes anything."""
    if flag:
        preflight(SOURMASH_TOOL)


def sourmash_for_sketches(flag: bool | None, stage: str, logger: logging.Logger) -> dict[str, str]:
    """Resolve --sketch/--no-sketch of a writer stage to sourmash's versions, or {}.

    ``flag`` None (the default) sketches when sourmash can run and otherwise
    logs one INFO line naming the flag; True requires sourmash
    (MissingBinaryError otherwise); False never sketches.
    """
    if flag is False:
        return {}
    if flag is None and not tool_available(SOURMASH_TOOL):
        logger.info(
            "%s: sourmash not found; genome sketches not written (--sketch requires them; "
            "'repgenr sketch' adds them later).",
            stage,
        )
        return {}
    try:
        return preflight(SOURMASH_TOOL)
    except MissingBinaryError as exc:
        if flag:
            raise
        logger.info("%s: genome sketches not written: %s", stage, exc)
        return {}


def sketch_stage_genomes(
    ctx: WorkdirContext,
    flag: bool | None,
    stage: str,
    logger: logging.Logger,
    *,
    threads: int | None = None,
    only: Collection[Path] | None = None,
    reuse: Mapping[str, SketchSource] | None = None,
) -> tuple[dict[str, Any], dict[str, str]]:
    """The sketch step at the end of a genome-writing stage.

    Prunes sketches of genomes that left the set (always, even with
    ``--no-sketch``), then sketches the missing and stale ones when sourmash
    can run: every genome of the set, or only the genome files in ``only``
    (assemble --append). Returns the summary for the stage record and
    sourmash's versions. A sketch failure under the default fails nothing:
    sketches are not deliverables, so it is logged and the stage finishes;
    with an explicit ``--sketch`` it is raised.
    """
    try:
        removed = remove_stale(ctx, logger)
    except OSError as exc:
        if flag:
            raise
        _warn_stopped(stage, exc, ctx, logger)
        removed = 0
    versions = sourmash_for_sketches(flag, stage, logger)
    if not versions:
        state = "off" if flag is False else "unavailable"
        return {"state": state, "removed": removed}, {}
    targets = sketch_targets(ctx)
    if only is not None:
        wanted = {Path(p).name for p in only}
        targets = [t for t in targets if t.path.name in wanted]
    try:
        summary = sketch_genomes(ctx, targets, threads or default_threads(), logger, reuse=reuse)
    except (ToolExecutionError, WorkdirError, OSError) as exc:
        if flag or process.stop_requested.is_set():
            raise
        _warn_stopped(stage, exc, ctx, logger)
        return {"state": "failed", "removed": removed}, versions
    summary.removed = removed
    logger.info("%s", summary.line())
    return {"state": "done", **summary.as_dict()}, versions


def _warn_stopped(
    stage: str, exc: BaseException, ctx: WorkdirContext, logger: logging.Logger
) -> None:
    logger.warning(
        "%s: genome sketching stopped (%s); the genome set is complete without "
        "sketches, and 'repgenr sketch -wd %s' builds the missing ones.",
        stage,
        exc,
        ctx.workdir,
    )


def sources_from_workdir(workdir: Path) -> dict[str, SketchSource]:
    """Reusable sketches of an earlier workdir, by record name (for ingest --from-workdir).

    Read from its manifest (read-only); a workdir without sketches, without a
    manifest, or with a pre-v4 manifest offers none.
    """
    from .manifest import MANIFEST_FILENAME, Manifest

    directory = sketches_dir(workdir)
    path = workdir / MANIFEST_FILENAME
    if not directory.is_dir() or not path.is_file():
        return {}
    try:
        manifest = Manifest.open_readonly(path)
    except WorkdirError:
        return {}
    try:
        records = manifest.sketch_records()
    except (sqlite3.Error, WorkdirError):  # an unreadable source manifest: sketch instead
        records = {}
    finally:
        manifest.close()
    out: dict[str, SketchSource] = {}
    for rec in records.values():
        if not rec.file:
            continue
        sketch = workdir / rec.file
        name = Path(rec.file).name
        if name.endswith(SKETCH_SUFFIX):
            out[name[: -len(SKETCH_SUFFIX)]] = SketchSource(sketch, rec)
    return out
