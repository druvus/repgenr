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

import hashlib
import logging
import re
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


_MD5_RE = re.compile(r"^([0-9a-fA-F]{32})\s+\.?/?(.+)$")


def verify_md5_manifest(
    path: Path, manifest_url: str, *, logger: logging.Logger
) -> bool:
    """Verify ``path`` against an md5sum-format manifest published beside it.

    Returns True when the checksum matches. An unavailable manifest or a file
    not listed in it is logged and skipped (False) -- upstream layouts vary --
    but an actual mismatch raises :class:`WorkdirError`: the download is
    corrupt and must not be parsed.
    """
    try:
        text = get_text(manifest_url)
    except WorkdirError as exc:
        logger.warning("Checksum manifest unavailable (%s); skipping verification", exc)
        return False
    expected = None
    for line in text.splitlines():
        m = _MD5_RE.match(line.strip())
        if m and Path(m.group(2)).name == path.name:
            expected = m.group(1).lower()
            break
    if expected is None:
        logger.warning("%s not listed in %s; skipping verification", path.name, manifest_url)
        return False
    digest = hashlib.md5()
    with open(path, "rb") as fo:
        for chunk in iter(lambda: fo.read(_CHUNK), b""):
            digest.update(chunk)
    actual = digest.hexdigest()
    if actual != expected:
        raise WorkdirError(
            f"Checksum mismatch for {path.name}: expected {expected}, got {actual}. "
            "The download is corrupt; delete it and re-run."
        )
    logger.info("Checksum verified for %s", path.name)
    return True
