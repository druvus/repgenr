"""Logging configuration.

Replaces the old shell ``tee`` plumbing in ``repgenr.py``. Each run logs to the
console and appends to ``<workdir>/repgenr.log`` with timestamps, so the full
output of every subprocess is captured in one place.

A message can carry a short console form (:func:`console_extra`): the console
prints that form with a pointer to the run log, and the run log keeps the full
text. Tool command lines with thousands of paths and long name lists use it,
so they do not fill the terminal.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence
from pathlib import Path
from typing import Any

_LOGGER_NAME = "repgenr"
_FORMAT = "%(asctime)s %(levelname)s %(message)s"
_DATEFMT = "%Y-%m-%d %H:%M:%S"
# LogRecord attribute that holds the console form of a message.
_CONSOLE_ATTR = "repgenr_console"
_DEFAULT_LOG = "repgenr.log"
# Terminal width a shortened console line is fitted to, header and suffix included.
CONSOLE_COLUMNS = 120


def console_suffix(log_name: str = _DEFAULT_LOG) -> str:
    """The pointer appended to a shortened console line."""
    return f" (full text in {log_name})"


def console_room(levelname: str = "INFO", log_name: str = _DEFAULT_LOG) -> int:
    """Characters left for a short console form so the line fits ``CONSOLE_COLUMNS``.

    The line is the timestamp, the level name, the short form and
    :func:`console_suffix`.
    """
    header = len("2026-01-01 00:00:00") + 1 + len(levelname) + 1  # matches _DATEFMT
    return CONSOLE_COLUMNS - header - len(console_suffix(log_name))


def console_extra(short: str) -> dict[str, Any]:
    """``extra=`` for a log call whose console form is ``short``.

    The run log keeps the full message. Without a run log (no working
    directory, or the file could not be opened) the console prints the full
    message, so the text is never lost.
    """
    return {_CONSOLE_ATTR: short}


def capped_names(names: Sequence[str], limit: int = 5) -> str:
    """The first ``limit`` names, comma-separated, then "and N more"."""
    shown = ", ".join(names[:limit])
    if len(names) > limit:
        return f"{shown} and {len(names) - limit} more"
    return shown


class _ConsoleFormatter(logging.Formatter):
    """Prints the short console form of a record when it carries one."""

    def __init__(self, log_name: str) -> None:
        super().__init__(_FORMAT, datefmt=_DATEFMT)
        self._suffix = console_suffix(log_name)

    def format(self, record: logging.LogRecord) -> str:
        short = getattr(record, _CONSOLE_ATTR, None)
        if short is None:
            return super().format(record)
        # A copy: the file handler formats the same record with the full text.
        shown = logging.makeLogRecord(record.__dict__)
        shown.msg = f"{short}{self._suffix}"
        shown.args = None
        return super().format(shown)


def configure_logging(
    workdir: str | os.PathLike[str] | None = None,
    *,
    level: int = logging.INFO,
    log_filename: str = _DEFAULT_LOG,
) -> logging.Logger:
    """Configure and return the ``repgenr`` logger.

    Idempotent: repeated calls replace handlers rather than stacking them. When
    ``workdir`` is given (and exists or can be created), a file handler appends
    to ``<workdir>/<log_filename>``.
    """
    logger = logging.getLogger(_LOGGER_NAME)
    # The logger passes everything to the handlers; each handler filters: the
    # console respects the requested level (INFO/--quiet/--verbose) while the file
    # always keeps full DEBUG detail, so an unexpected-error traceback is captured
    # in the run log without cluttering the console.
    logger.setLevel(logging.DEBUG)
    for handler in logger.handlers:
        handler.close()
    logger.handlers.clear()
    logger.propagate = False

    formatter = logging.Formatter(_FORMAT, datefmt=_DATEFMT)

    console = logging.StreamHandler()
    console.setFormatter(formatter)
    console.setLevel(level)
    logger.addHandler(console)

    if workdir is not None:
        wd = Path(workdir)
        try:
            wd.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(wd / log_filename)
            file_handler.setFormatter(formatter)
            file_handler.setLevel(logging.DEBUG)
            logger.addHandler(file_handler)
        except OSError:
            logger.warning("Could not open log file under %s; console only", wd)
        else:
            # The full text of a shortened line is in the run log.
            console.setFormatter(_ConsoleFormatter(log_filename))

    return logger
