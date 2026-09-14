"""Input-completeness guards between pipeline stages.

The resume fingerprint detects that inputs *changed*; these checks detect that
inputs are *incomplete* -- the state a crashed upstream stage leaves behind (a
prefix of ``genomes/``, or ``derep/representatives/`` out of step with
``clusters.tsv``). A consuming stage refuses to run on such inputs unless the
user passes ``--allow-incomplete``, which downgrades the refusal to a warning
(for deliberately hand-pruned genome sets).
"""

from __future__ import annotations

import logging
from pathlib import Path

from .contracts import (
    EXCUSED_RUNS_TSV,
    MISSING_ACCESSIONS_TXT,
    SELECTION_TSV,
    list_fasta,
    read_clusters,
    read_excused_runs,
    read_selection,
)
from .errors import UserInputError, WorkdirError


def looks_like_fasta(path: Path) -> bool:
    """Cheap first-bytes check: does this file plausibly hold FASTA data?

    Catches HTML error pages served with HTTP 200 and empty/absent files
    without parsing whole (potentially multi-GB) genomes.
    """
    try:
        with open(path, "rb") as fo:
            head = fo.read(4096)
    except OSError:
        return False
    return head.lstrip().startswith(b">")


def excused_accessions(workdir: Path) -> set[str]:
    """Selected accessions no genome is expected for.

    The union of ``missing_accessions.txt`` (the genome stage asked NCBI and
    got nothing) and the runs in ``excused_runs.tsv`` (the assemble stage could
    not fetch, assemble or pass a run). Both are recorded deliberately, so the
    completeness guard reports the rest of a shortfall without them.
    """
    excused: set[str] = set()
    path = workdir / MISSING_ACCESSIONS_TXT
    if path.exists():
        excused |= {
            line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()
        }
    runs = workdir / EXCUSED_RUNS_TSV
    if runs.exists():
        excused |= {row.run_accession for row in read_excused_runs(runs)}
    return excused


def check_genome_completeness(
    genomes_dir: Path,
    workdir: Path,
    *,
    logger: logging.Logger,
    allow_incomplete: bool = False,
) -> list[str]:
    """Verify ``genomes_dir`` holds every genome the selection promised.

    Expected filenames come from ``selection.tsv`` (outgroup rows excluded);
    accessions NCBI legitimately returned nothing for (recorded by the genome
    stage in ``missing_accessions.txt``) are excused. Returns the shortfall;
    raises :class:`WorkdirError` on a shortfall unless ``allow_incomplete``.
    """
    selection = workdir / SELECTION_TSV
    if not selection.exists():
        logger.debug("No %s; skipping genome completeness check", SELECTION_TSV)
        return []
    excused = excused_accessions(workdir)
    expected = {
        row.filename
        for row in read_selection(selection)
        if not row.is_outgroup and row.accession not in excused
    }
    present = {p.name for p in list_fasta(genomes_dir)}
    shortfall = sorted(expected - present)
    if not shortfall:
        return []
    message = (
        f"{genomes_dir} is missing {len(shortfall)} of {len(expected)} selected "
        f"genomes (e.g. {', '.join(shortfall[:3])}). The genome stage likely "
        "crashed or was interrupted; re-run it (or pass --allow-incomplete to "
        "proceed with the partial set)."
    )
    if allow_incomplete:
        logger.warning("%s", message)
        return shortfall
    raise WorkdirError(message)


def check_representatives_consistency(
    representatives_dir: Path,
    clusters_tsv: Path,
    *,
    logger: logging.Logger,
    allow_incomplete: bool = False,
) -> list[str]:
    """Verify the representatives directory matches ``clusters.tsv``.

    A crash inside the dereplicate stage's output write leaves a prefix of the
    representatives next to the previous run's cluster table; the mismatch is
    detected here before a tree is built over it. Returns the missing
    representative filenames; raises unless ``allow_incomplete``.
    """
    if not clusters_tsv.exists():
        logger.debug("No %s; skipping representatives consistency check", clusters_tsv)
        return []
    expected = set(read_clusters(clusters_tsv))
    present = {p.name for p in list_fasta(representatives_dir)}
    shortfall = sorted(expected - present)
    if not shortfall:
        return []
    message = (
        f"{representatives_dir} is out of step with {clusters_tsv}: "
        f"{len(shortfall)} representative(s) listed but absent "
        f"(e.g. {', '.join(shortfall[:3])}). The dereplicate stage likely "
        "crashed mid-write; re-run it (or pass --allow-incomplete)."
    )
    if allow_incomplete:
        logger.warning("%s", message)
        return shortfall
    raise WorkdirError(message)


# Manifest sources written by stages that add to a working directory rather
# than select it; a re-selection would drop their rows and prune their files.
FOREIGN_SOURCES = frozenset({"sra"})


def refuse_foreign_rows(ctx, stage: str, *, drop_foreign: bool, logger: logging.Logger) -> None:
    """Stop ``stage`` from silently discarding genomes another entry path appended.

    ``metadata`` and ``ingest`` replace the manifest and the genome stage prunes
    what the manifest no longer lists, so reads-derived genomes added with
    ``assemble --append`` would vanish. Refuse unless the user passed
    ``--drop-foreign``, in which case say what is being dropped.
    """
    try:
        foreign = [
            g.accession
            for g in ctx.manifest.all_genomes(include_outgroup=True)
            if g.source in FOREIGN_SOURCES
        ]
    except Exception:  # no manifest yet, or an unreadable one: nothing to protect
        return
    if not foreign:
        return
    listing = ", ".join(foreign[:5]) + (" ..." if len(foreign) > 5 else "")
    if not drop_foreign:
        raise UserInputError(
            f"The working directory holds {len(foreign)} genome(s) appended from sequencing "
            f"runs ({listing}); re-running {stage} would drop them and the genome stage "
            "would prune their files. Pass --drop-foreign to discard them, or run "
            "`assemble --append` again afterwards to restore them."
        )
    logger.warning(
        "Dropping %d appended genome(s) (%s) as --drop-foreign asked", len(foreign), listing
    )
