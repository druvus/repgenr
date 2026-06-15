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
import subprocess
from collections import deque
from collections.abc import Mapping, Sequence
from pathlib import Path

from .errors import ToolExecutionError

_DEFAULT_TAIL = 50


def run(
    command: Sequence[str | os.PathLike[str]],
    *,
    logger: logging.Logger,
    cwd: str | os.PathLike[str] | None = None,
    env: Mapping[str, str] | None = None,
    check: bool = True,
    stdout_path: str | os.PathLike[str] | None = None,
    log_prefix: str | None = None,
) -> int:
    """Run ``command`` (an argument vector) without a shell.

    Output is line-streamed to ``logger`` at INFO. When ``stdout_path`` is given,
    stdout is written there instead (stderr still goes to the logger) -- use this
    for tools that emit their result on stdout (e.g. ``mashtree``).

    Returns the process exit code. Raises :class:`ToolExecutionError` on a
    non-zero exit when ``check`` is True.
    """
    cmd = [str(part) for part in command]
    prefix = f"[{log_prefix}] " if log_prefix else ""
    logger.info("%s$ %s", prefix, " ".join(cmd))

    full_env = {**os.environ, **env} if env else None
    tail: deque[str] = deque(maxlen=_DEFAULT_TAIL)

    out_handle = open(stdout_path, "w") if stdout_path is not None else None
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd) if cwd is not None else None,
            env=full_env,
            stdout=(out_handle if out_handle is not None else subprocess.PIPE),
            stderr=subprocess.STDOUT if out_handle is None else subprocess.PIPE,
            text=True,
            bufsize=1,
        )
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
        if out_handle is not None:
            out_handle.close()

    if check and returncode != 0:
        raise ToolExecutionError(cmd, returncode, output="\n".join(tail))
    return returncode


def write_fofn(paths: Sequence[str | os.PathLike[str]], dest: str | os.PathLike[str]) -> Path:
    """Write a file-of-filenames (one absolute path per line) and return its path.

    Use this instead of shell globs when handing a large genome set to a tool.
    """
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(dest_path, "w") as fo:
        for p in paths:
            fo.write(f"{os.fspath(Path(p).resolve())}\n")
    return dest_path
