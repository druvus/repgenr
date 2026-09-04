"""Input digests for the stage-resume fingerprint.

A stage's skip decision must reflect its actual inputs, not only its CLI
parameters, so re-running an upstream stage invalidates downstream stages.
Directories of genomes can hold thousands of files and tens of gigabytes, so
they are digested from file metadata (name, size, mtime_ns) rather than
content: one ``stat`` per file, catching adds, removes, renames, and rewrites
(every stage that regenerates a file bumps its mtime). Small contract files
(selection.tsv, clusters.tsv, tree.nwk) are digested by content. A same-size,
mtime-preserving in-place edit is therefore not detected; ``--force`` covers
that deliberate case.
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .manifest import Manifest

# Stable sentinel for a missing input path: a stage whose inputs do not exist
# yet degrades to params-only fingerprinting instead of erroring.
ABSENT = "absent"

_CHUNK = 1 << 20  # 1 MiB streaming chunks, matching core.http


def file_digest(path: Path) -> str:
    """Chunked sha256 of a file's content; ``ABSENT`` when the file is missing."""
    if not path.is_file():
        return ABSENT
    digest = hashlib.sha256()
    with open(path, "rb") as fo:
        for chunk in iter(lambda: fo.read(_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dir_stat_digest(path: Path) -> str:
    """sha256 over sorted (name, size, mtime_ns) of a directory's files.

    Flat (non-recursive: the contract directories are flat) and dotfile-filtered,
    matching :func:`repgenr.core.contracts.list_fasta` semantics. ``ABSENT`` when
    the directory is missing.
    """
    if not path.is_dir():
        return ABSENT
    digest = hashlib.sha256()
    for entry in sorted(path.iterdir()):
        if entry.name.startswith(".") or not entry.is_file():
            continue
        st = entry.stat()
        digest.update(f"{entry.name}\0{st.st_size}\0{st.st_mtime_ns}\n".encode())
    return digest.hexdigest()


def path_digest(path: Path) -> str:
    """Digest a path: directories by metadata, files by content."""
    if path.is_dir():
        return dir_stat_digest(path)
    return file_digest(path)


def inputs_digest(workdir: Path, paths: Iterable[Path]) -> dict[str, str]:
    """Digest each input path, keyed by its workdir-relative posix path."""
    out: dict[str, str] = {}
    for path in paths:
        try:
            key = path.relative_to(workdir).as_posix()
        except ValueError:
            key = path.as_posix()
        out[key] = path_digest(path)
    return out


def manifest_digest(manifest: Manifest, *, include_derep: bool = True) -> str:
    """sha256 of the manifest's genome rows, independent of insertion order.

    Hashes ordered query results, never the SQLite file bytes (which are not
    byte-stable across identical logical content). Includes CheckM completeness
    and contamination (consumed by the dereplicate keeper) alongside taxonomy.

    ``include_derep`` adds derep_status/representative to the hash; it must be
    False for a stage that WRITES those columns itself (dereplicate), or the
    stage would invalidate its own resume on every first re-run: the digest
    computed before the stage runs would never match the one computed after,
    since the run just changed the very columns being hashed. tree2tax (which
    only READS derep status) keeps the default True.
    """
    digest = hashlib.sha256()
    rows = sorted(
        (
            r.accession,
            r.filename or "",
            r.family or "",
            r.genus or "",
            r.species or "",
            str(r.is_outgroup),
            *((r.derep_status or "", r.representative or "") if include_derep else ()),
            "" if r.completeness is None else repr(r.completeness),
            "" if r.contamination is None else repr(r.contamination),
        )
        for r in manifest.all_genomes(include_outgroup=True)
    )
    for row in rows:
        digest.update(("\0".join(row) + "\n").encode())
    return digest.hexdigest()


# Whether a stage's manifest digest includes derep_status/representative.
# "dereplicate" WRITES those columns itself (_update_manifest), so including
# them would invalidate its own resume on every first re-run (see
# manifest_digest's docstring). Every other manifest-input stage (tree2tax,
# which only READS derep status) keeps the default True.
_MANIFEST_DIGEST_INCLUDE_DEREP: dict[str, bool] = {"dereplicate": False}


def manifest_digest_for_stage(stage: str, manifest: Manifest) -> str:
    """manifest_digest() with the per-stage include_derep decision applied.

    This is the single source of truth for that decision: both the resume
    fingerprint (cli/base.py's _stage_input_digests, which stamps a stage's
    completion record) and `repgenr doctor`'s stale-input check
    (core/doctor.py's _check_stale_inputs, which recomputes the same digest to
    compare against that record) must call this instead of manifest_digest
    directly. Calling manifest_digest with a hardcoded or differing
    include_derep in only one of the two places means the two digests can
    never agree, so a completed stage (or a healthy one, on every future
    invocation) would falsely and permanently report as stale.
    """
    return manifest_digest(manifest, include_derep=_MANIFEST_DIGEST_INCLUDE_DEREP.get(stage, True))
