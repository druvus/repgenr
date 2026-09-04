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

_module_logger = logging.getLogger(__name__)

# A malformed REPGENR_SUBPROCESS_TIMEOUT is reported once, not on every
# subprocess launch.
_warned_bad_timeout = False


def _default_timeout() -> float | None:
    """Global subprocess timeout (seconds) from ``REPGENR_SUBPROCESS_TIMEOUT``.

    Unset (the default) means no timeout, preserving prior behavior; operators
    can cap every external tool with one environment variable. A value that is
    not a positive number is ignored (no timeout) with a single warning.
    """
    global _warned_bad_timeout
    raw = os.environ.get("REPGENR_SUBPROCESS_TIMEOUT")
    if not raw:
        return None
    try:
        value = float(raw)
    except ValueError:
        value = -1.0
    if value > 0:
        return value
    if not _warned_bad_timeout:
        _warned_bad_timeout = True
        _module_logger.warning(
            "Ignoring REPGENR_SUBPROCESS_TIMEOUT=%r: expected a positive "
            "number of seconds; running without a subprocess timeout.",
            raw,
        )
    return None


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

    The command line is logged at INFO; the tool's own output is line-streamed
    at DEBUG (progress bars and per-file chatter would otherwise dominate the
    log), and the last lines are kept for the error message on failure. A line
    a tool redraws in place with carriage returns is logged once, in its final
    state. When ``stdout_path`` is given, stdout is written there instead
    (stderr still goes to the logger) -- use this for tools that emit their
    result on stdout (e.g. ``mashtree``).

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

    # stdout goes to a temp sibling first and is published only when the tool
    # did not fail: several adapters point stdout_path at a final deliverable
    # (e.g. tree/tree.nwk), and opening that in "w" mode would truncate the
    # previous good file the moment a doomed rebuild starts. os.devnull is
    # written directly (it has no meaningful siblings).
    out_target: Path | None = None
    out_tmp: Path | None = None
    if stdout_path is not None:
        if str(stdout_path) == os.devnull:
            out_tmp = Path(os.devnull)
        else:
            out_target = Path(stdout_path)
            out_tmp = out_target.with_name(out_target.name + ".part")
    out_handle = open(out_tmp, "w", encoding="utf-8") if out_tmp is not None else None
    timer: threading.Timer | None = None
    timed_out = False
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd) if cwd is not None else None,
            env=full_env,
            stdout=(out_handle if out_handle is not None else subprocess.PIPE),
            stderr=subprocess.STDOUT if out_handle is None else subprocess.PIPE,
            # Binary so that a "\r" inside a line is ours to interpret; text
            # mode's universal newlines would split every progress-bar redraw
            # into its own line.
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
                # A progress bar redraws one line with "\r"; keep its last frame.
                text = raw.decode("utf-8", errors="replace")
                line = text.rstrip("\r\n").rsplit("\r", 1)[-1]
                if not line:
                    continue
                tail.append(line)
                logger.debug("%s%s", prefix, line)
        returncode = proc.wait()
    except BaseException:
        if out_target is not None and out_tmp is not None:
            out_tmp.unlink(missing_ok=True)
        raise
    finally:
        if timer is not None:
            timer.cancel()
        if out_handle is not None:
            out_handle.close()

    failed = timed_out or (check and returncode != 0)
    if out_target is not None and out_tmp is not None:
        if failed:
            # Discard the partial capture; a previous good file stays intact.
            out_tmp.unlink(missing_ok=True)
        else:
            os.replace(out_tmp, out_target)
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


# Close to the common ARG_MAX of 1 MiB (macOS and most Linux configurations);
# execve also counts the environment, so warn with headroom.
_ARGV_BYTES_WARN = 900_000


def warn_argv_bytes(
    tool: str, argv: Sequence[str | os.PathLike[str]], logger: logging.Logger
) -> None:
    """Warn when an argv's byte size approaches the OS ARG_MAX limit.

    Some tools take every genome path on argv (mashtree, sibeliaz,
    snippy-core); at thousands of genomes the exec can fail with E2BIG.
    Use before run_tool for such invocations. Reduce the set (dereplicate
    first, or --process-size for dereplication) when the warning fires.
    """
    total = sum(len(os.fsencode(os.fspath(a))) + 1 for a in argv)
    if total >= _ARGV_BYTES_WARN:
        logger.warning(
            "%s receives ~%d KB across %d command-line arguments; the OS "
            "ARG_MAX limit (commonly 1 MB) may be exceeded and the tool may "
            "fail to launch. Reduce the genome set per call.",
            tool,
            total // 1000,
            len(argv),
        )


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
    with open(dest_path, "w", encoding="utf-8") as fo:
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
    # the content and still shares the inode (no extra disk). realpath() walks the
    # whole path with lstat, so only pay it when src is actually a symlink (the
    # common case -- staging plain genome files -- skips it).
    src_path = os.fspath(src)
    src_s = os.path.realpath(src_path) if os.path.islink(src_path) else src_path
    dst_s = os.fspath(dst)
    try:
        os.link(src_s, dst_s)
    except OSError:
        shutil.copy2(src_s, dst_s)
