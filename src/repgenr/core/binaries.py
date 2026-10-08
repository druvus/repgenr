"""External binary discovery and version preflight.

Every adapter declares the binaries it needs via :class:`BinarySpec` in its
``ToolCapabilities``. Before a tool runs, :func:`check_binaries` confirms each
binary is on PATH and, when a minimum version is declared, that it is new
enough. Resolved versions are returned so the stage can record them in
``repgenr.yaml`` provenance.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from .errors import MissingBinaryError

# A version is not part of a longer token: "python3.12" in an interpreter path
# or "GLIBC_2.17" in a loader error is not the tool's version. A leading "v"
# ("dRep v3.4.5", "hal2maf v2.2") is allowed.
_VERSION_RE = re.compile(r"(?<![0-9A-UW-Za-uw-z_./])(\d+)\.(\d+)(?:\.(\d+))?")


@dataclass(frozen=True)
class BinarySpec:
    """A required external executable.

    ``version_args`` is the argument vector that prints a version (e.g.
    ``("--version",)``). ``min_version`` is an optional ``(major, minor, patch)``
    or dotted string requirement.
    """

    name: str
    version_args: tuple[str, ...] = field(default=("--version",))
    min_version: str | None = None
    # When the version cannot be parsed, lenient specs pass (avoid false rejects);
    # strict specs fail -- use this for tools where a legacy build that does not
    # answer its version flag can shadow a modern one on PATH (samtools/bcftools).
    strict_version: bool = False


def _parse_version(text: str) -> tuple[int, int, int] | None:
    match = _VERSION_RE.search(text)
    if not match:
        return None
    major, minor, patch = match.group(1), match.group(2), match.group(3)
    return (int(major), int(minor), int(patch or 0))


def _query_version(name: str, version_args: tuple[str, ...]) -> str | None:
    try:
        proc = subprocess.run(
            [name, *version_args],
            capture_output=True,
            text=True,
            timeout=30,
            # A tool queried without arguments (FastTree) may otherwise wait
            # on an inherited terminal for input until the timeout.
            stdin=subprocess.DEVNULL,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    blob = (proc.stdout or "") + (proc.stderr or "")
    # A tool that crashed while answering (a Python traceback, exit != 0) has
    # not reported a version; numbers in its traceback are file paths.
    if proc.returncode != 0 and "Traceback (most recent call last)" in blob:
        return None
    parsed = _parse_version(blob)
    if parsed:
        return ".".join(map(str, parsed))
    # Without a version number, keep the first line only when the tool accepted
    # the flag: a non-zero exit means the line is an error message (sibeliaz
    # has no version flag and answers `-v` with "illegal option"), not a version.
    if proc.returncode != 0:
        return None
    return blob.strip().splitlines()[0] if blob.strip() else None


def _owns_binary(data: object, name: str) -> bool:
    """True when a conda-meta record lists ``bin/<name>`` among its files."""
    return isinstance(data, dict) and f"bin/{name}" in (data.get("files") or ())


def _metadata_version(name: str, path: str | None = None) -> str | None:
    """The version of the conda package that installed ``name``, if any.

    Used when the binary does not report a version itself (SibeliaZ has no
    version flag). The package records sit in ``<prefix>/conda-meta`` next to
    the ``bin`` directory that holds the binary; the path is not resolved,
    since conda links files into the prefix. The record named after the binary
    is tried first, then any record whose file list contains ``bin/<name>``.
    """
    found = path or shutil.which(name)
    if found is None:
        return None
    meta = Path(found).parent.parent / "conda-meta"
    try:
        if not meta.is_dir():
            return None
        named = sorted(meta.glob(f"{name}-*.json"))
        others = [p for p in sorted(meta.glob("*.json")) if p not in set(named)]
        for record in (*named, *others):
            try:
                text = record.read_text(encoding="utf-8")
                if f'"bin/{name}"' not in text:
                    continue
                data = json.loads(text)
            except (OSError, ValueError):
                continue
            if _owns_binary(data, name):
                version = data.get("version")
                if isinstance(version, str) and version.strip():
                    return version.strip()
    except OSError:
        return None
    return None


def check_binaries(specs: tuple[BinarySpec, ...]) -> dict[str, str]:
    """Confirm all ``specs`` are present (and new enough). Return name -> version.

    Raises :class:`MissingBinaryError` listing every missing or too-old binary.
    """
    versions: dict[str, str] = {}
    problems: list[str] = []

    for spec in specs:
        found = shutil.which(spec.name)
        if found is None:
            problems.append(f"{spec.name}: not found on PATH")
            continue

        reported = _query_version(spec.name, spec.version_args)
        if _parse_version(reported or "") is None:
            # No version number from the tool itself: the conda package
            # record, when there is one, still names the installed version.
            reported = _metadata_version(spec.name, found) or reported
        versions[spec.name] = reported or "unknown"

        if spec.min_version is not None:
            have = _parse_version(reported or "")
            want = _parse_version(spec.min_version)
            if have is None and want is not None and spec.strict_version:
                got = f"got {reported!r}" if reported else "no version reported"
                problems.append(
                    f"{spec.name}: could not read a version ({got}); need "
                    f">= {spec.min_version}. A legacy build that does not answer its "
                    "version flag may be shadowing it on PATH -- use a modern environment."
                )
            elif have is not None and want is not None and have < want:
                problems.append(f"{spec.name}: version {reported} < required {spec.min_version}")

    if problems:
        raise MissingBinaryError(
            "Required external tools are missing or outdated:\n  " + "\n  ".join(problems)
        )
    return versions
