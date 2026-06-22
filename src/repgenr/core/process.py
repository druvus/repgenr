"""Safe subprocess execution.

Replaces the old ``subprocess.call(' '.join(cmd), shell=True)`` pattern and the
shell globs (``genomes/*.fasta``) that break past ``ARG_MAX`` for large genome
sets. Commands are always passed as argument vectors (no shell), output is
streamed to the logger, and a non-zero exit raises :class:`ToolExecutionError`.

For tools that genuinely need a large list of input files, write the list to a
file-of-filenames with :func:`write_fofn` and pass that path instead of a glob.
"""

from __future__ import annotations

import logging
import os
import shutil
import signal
import subprocess
import threading
import zipfile
from collections import deque
from collections.abc import Mapping, Sequence
from pathlib import Path

from .errors import ToolExecutionError, WorkdirError

_DEFAULT_TAIL = 50


def _default_timeout() -> float | None:
    """Global subprocess timeout (seconds) from ``REPGENR_SUBPROCESS_TIMEOUT``.

    Unset (the default) means no timeout, preserving prior behavior; operators
    can cap every external tool with one environment variable.
    """
    raw = os.environ.get("REPGENR_SUBPROCESS_TIMEOUT")
    if not raw:
        return None
    try:
        value = float(raw)
    except ValueError:
        return None
    return value if value > 0 else None


def run(
    command: Sequence[str | os.PathLike[str]],
    *,
    logger: logging.Logger,
    cwd: str | os.PathLike[str] | None = None,
    env: Mapping[str, str] | None = None,
    check: bool = True,
    stdout_path: str | os.PathLike[str] | None = None,
    log_prefix: str | None = None,
    timeout: float | None = None,
) -> int:
    """Run ``command`` (an argument vector) without a shell.

    Output is line-streamed to ``logger`` at INFO. When ``stdout_path`` is given,
    stdout is written there instead (stderr still goes to the logger) -- use this
    for tools that emit their result on stdout (e.g. ``mashtree``).

    ``timeout`` (seconds) caps the run: on expiry the whole process group is
    killed and :class:`ToolExecutionError` is raised, so a hung tool (or a stuck
    ``docker pull``) cannot wedge the pipeline. It defaults to the
    ``REPGENR_SUBPROCESS_TIMEOUT`` environment variable (unset = no timeout).

    Returns the process exit code. Raises :class:`ToolExecutionError` on a
    non-zero exit when ``check`` is True.
    """
    cmd = [str(part) for part in command]
    prefix = f"[{log_prefix}] " if log_prefix else ""
    logger.info("%s$ %s", prefix, " ".join(cmd))

    full_env = {**os.environ, **env} if env else None
    tail: deque[str] = deque(maxlen=_DEFAULT_TAIL)
    limit = timeout if timeout is not None else _default_timeout()

    out_handle = open(stdout_path, "w") if stdout_path is not None else None
    timer: threading.Timer | None = None
    timed_out = False
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd) if cwd is not None else None,
            env=full_env,
            stdout=(out_handle if out_handle is not None else subprocess.PIPE),
            stderr=subprocess.STDOUT if out_handle is None else subprocess.PIPE,
            text=True,
            bufsize=1,
            # Own process group so a timeout can kill the tool and its children.
            start_new_session=limit is not None,
        )

        if limit is not None:
            def _kill() -> None:
                nonlocal timed_out
                timed_out = True
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                except (ProcessLookupError, PermissionError, OSError):
                    pass
            timer = threading.Timer(limit, _kill)
            timer.start()

        stream = proc.stdout if out_handle is None else proc.stderr
        if stream is not None:
            for raw in stream:
                line = raw.rstrip("\n")
                if not line:
                    continue
                tail.append(line)
                logger.info("%s%s", prefix, line)
        returncode = proc.wait()
    finally:
        if timer is not None:
            timer.cancel()
        if out_handle is not None:
            out_handle.close()

    if timed_out:
        tail.append(f"[killed after {limit}s timeout]")
        raise ToolExecutionError(cmd, returncode, output="\n".join(tail))
    if check and returncode != 0:
        raise ToolExecutionError(cmd, returncode, output="\n".join(tail))
    return returncode


def unzip(zip_path: str | os.PathLike[str], dest: str | os.PathLike[str]) -> None:
    """Extract a zip, turning a truncated/corrupt archive into a clear error.

    A `datasets` download that was cut short (network reset, full disk) leaves a
    bad zip; ``zipfile`` then raises ``BadZipFile`` which would surface as a raw
    traceback. Map it to :class:`WorkdirError` naming the file so a re-run retries.
    """
    try:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(dest)
    except zipfile.BadZipFile as exc:
        raise WorkdirError(
            f"Corrupt or truncated download: {zip_path} ({exc}). Re-run to retry."
        ) from exc


def write_fofn(paths: Sequence[str | os.PathLike[str]], dest: str | os.PathLike[str]) -> Path:
    """Write a file-of-filenames (one absolute path per line) and return its path.

    Use this instead of shell globs when handing a large genome set to a tool.

    Paths are absolute but NOT symlink-resolved: the container backend binds
    un-resolved abspaths (macOS firmlinks resolve /Users -> /System/Volumes/Data,
    which is outside Docker's shared directories), so a tool reading this fofn
    inside a container must see the same un-resolved paths.
    """
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(dest_path, "w") as fo:
        for p in paths:
            fo.write(f"{os.path.abspath(os.fspath(p))}\n")
    return dest_path


def link_or_copy(src: str | os.PathLike[str], dst: str | os.PathLike[str]) -> None:
    """Stage ``src`` at ``dst`` cheaply: hardlink it, copying only as a fallback.

    Staging genomes into representatives/cluster dirs copies tens of GB at 1000s
    of genomes. A hardlink is instant and uses no extra disk; it shares the inode
    with the source, which is safe because these staged files are only read by
    downstream stages, never modified in place. Falls back to a real copy when
    the filesystem can't hardlink (cross-device, or exFAT/NTFS on the dev box).
    """
    # Resolve symlinks to the real file first. Tools such as skDER emit their
    # representative genomes as symlinks (often into a Nextflow-/container-staged
    # input tree); hardlinking the symlink itself -- which os.link does on macOS --
    # produces a broken, 0-byte staged file. Linking the real target instead keeps
    # the content and still shares the inode (no extra disk).
    src_s, dst_s = os.path.realpath(os.fspath(src)), os.fspath(dst)
    try:
        os.link(src_s, dst_s)
    except OSError:
        shutil.copy2(src_s, dst_s)
