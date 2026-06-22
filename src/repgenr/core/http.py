"""Shared HTTP client with retry/backoff for external services.

RepGenR talks to flaky, externally-operated services (the GTDB API and metadata
download, NCBI Entrez). Routing them through one process-wide ``requests.Session``
with a urllib3 retry policy gives every caller the same transient-failure
resilience -- retrying 429/5xx with exponential backoff and honoring
``Retry-After`` -- and the same clean error surface: a failed request raises
:class:`WorkdirError` naming the URL, never a bare ``requests`` traceback. The
streaming :func:`download` additionally verifies the byte count against
``Content-Length`` and writes through a ``.part`` file, so an interrupted
transfer never leaves a truncated file that later parses as if complete.
"""

from __future__ import annotations

import logging
from functools import lru_cache
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .errors import WorkdirError

_DEFAULT_TIMEOUT = 120
_CHUNK = 1 << 20  # 1 MiB streaming chunks

_RETRY = Retry(
    total=5,
    backoff_factor=1.0,  # 0s, 1s, 2s, 4s, 8s between attempts
    status_forcelist=(429, 500, 502, 503, 504),
    allowed_methods=frozenset({"GET", "HEAD"}),
    respect_retry_after_header=True,
    raise_on_status=False,
)


@lru_cache(maxsize=1)
def session() -> requests.Session:
    """Process-wide session with the retry policy mounted on http(s)."""
    s = requests.Session()
    adapter = HTTPAdapter(max_retries=_RETRY)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers.update({"User-Agent": "repgenr"})
    return s


def _get(url: str, *, params: dict | None, timeout: int) -> requests.Response:
    try:
        resp = session().get(url, params=params, timeout=timeout)
        resp.raise_for_status()
        return resp
    except requests.RequestException as exc:
        raise WorkdirError(f"HTTP request failed: {url} ({exc})") from exc


def get_json(url: str, *, params: dict | None = None, timeout: int = _DEFAULT_TIMEOUT) -> dict:
    """GET ``url`` and parse JSON; raises :class:`WorkdirError` on failure."""
    resp = _get(url, params=params, timeout=timeout)
    try:
        return resp.json()
    except ValueError as exc:
        raise WorkdirError(f"Expected JSON from {url} but could not parse it ({exc})") from exc


def get_text(url: str, *, params: dict | None = None, timeout: int = _DEFAULT_TIMEOUT) -> str:
    """GET ``url`` and return the response body as text (status-checked)."""
    return _get(url, params=params, timeout=timeout).text


def download(
    url: str,
    dest: Path,
    *,
    timeout: int = _DEFAULT_TIMEOUT,
    logger: logging.Logger | None = None,
) -> Path:
    """Stream ``url`` to ``dest``, verifying the size and writing atomically.

    Writes to ``<dest>.part`` and renames on success; a transfer that drops
    short of the server's ``Content-Length`` is deleted and raised as a
    :class:`WorkdirError` rather than left as a silently-truncated file.
    """
    dest = Path(dest)
    tmp = dest.with_name(dest.name + ".part")
    written = 0
    expected = 0
    try:
        with session().get(url, stream=True, timeout=timeout) as resp:
            resp.raise_for_status()
            expected = int(resp.headers.get("Content-Length", 0) or 0)
            with open(tmp, "wb") as fo:
                for chunk in resp.iter_content(chunk_size=_CHUNK):
                    fo.write(chunk)
                    written += len(chunk)
    except requests.RequestException as exc:
        tmp.unlink(missing_ok=True)
        raise WorkdirError(f"Download failed: {url} ({exc})") from exc

    if expected and written != expected:
        tmp.unlink(missing_ok=True)
        raise WorkdirError(
            f"Incomplete download: {url} got {written} of {expected} bytes."
        )
    tmp.replace(dest)
    if logger is not None:
        logger.info("Downloaded %s (%d bytes)", dest.name, written)
    return dest
