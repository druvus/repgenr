"""Logging configuration.

Replaces the old shell ``tee`` plumbing in ``repgenr.py``. Each run logs to the
console and appends to ``<workdir>/repgenr.log`` with timestamps, so the full
output of every subprocess is captured in one place.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

_LOGGER_NAME = "repgenr"
_FORMAT = "%(asctime)s %(levelname)s %(message)s"
_DATEFMT = "%Y-%m-%d %H:%M:%S"


def configure_logging(
    workdir: str | os.PathLike[str] | None = None,
    *,
    level: int = logging.INFO,
    log_filename: str = "repgenr.log",
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

    return logger


def get_logger() -> logging.Logger:
    """Return the shared ``repgenr`` logger (configure_logging if not yet set up)."""
    logger = logging.getLogger(_LOGGER_NAME)
    if not logger.handlers:
        return configure_logging()
    return logger
