"""Tests for the shared HTTP client (no network).

The session is replaced with a fake so the JSON/text/download helpers and their
error surfaces are exercised deterministically: status errors and unparseable
bodies become WorkdirError, and a short download is rejected with no leftover
file.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest
import requests

from repgenr.core import http
from repgenr.core.errors import WorkdirError

_LOG = logging.getLogger("test")
# Taken at import, before the autouse fixture in tests/conftest.py stubs it.
_REAL_PROBE_SESSION = http._probe_session


class _FakeResp:
    def __init__(self, *, json_data=None, text="", content=b"", headers=None, status_exc=None):
        self._json = json_data
        self.text = text
        self._content = content
        self.headers = headers or {}
        self._status_exc = status_exc

    def raise_for_status(self):
        if self._status_exc is not None:
            raise self._status_exc

    def json(self):
        if isinstance(self._json, ValueError):
            raise self._json
        return self._json

    def iter_content(self, chunk_size=1):
        for i in range(0, len(self._content), chunk_size):
            yield self._content[i : i + chunk_size]

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class _FakeSession:
    def __init__(self, resp):
        self._resp = resp

    def get(self, url, **kw):
        return self._resp


def _patch(monkeypatch, resp):
    monkeypatch.setattr(http, "session", lambda: _FakeSession(resp))


def test_get_json_ok(monkeypatch) -> None:
    _patch(monkeypatch, _FakeResp(json_data={"ok": 1}))
    assert http.get_json("https://x/y") == {"ok": 1}


def test_get_json_bad_body(monkeypatch) -> None:
    _patch(monkeypatch, _FakeResp(json_data=ValueError("nope")))
    with pytest.raises(WorkdirError, match="could not parse"):
        http.get_json("https://x/y")


def test_get_text_ok(monkeypatch) -> None:
    _patch(monkeypatch, _FakeResp(text="hello"))
    assert http.get_text("https://x/y") == "hello"


def test_status_error_becomes_workdir_error(monkeypatch) -> None:
    _patch(monkeypatch, _FakeResp(status_exc=requests.HTTPError("500")))
    with pytest.raises(WorkdirError, match="HTTP request failed"):
        http.get_text("https://x/y")


def test_download_ok(monkeypatch, tmp_path: Path) -> None:
    body = b"A" * 2048
    _patch(monkeypatch, _FakeResp(content=body, headers={"Content-Length": str(len(body))}))
    dest = tmp_path / "f.gz"
    http.download("https://x/f.gz", dest)
    assert dest.read_bytes() == body
    assert not (tmp_path / "f.gz.part").exists()  # temp renamed away


def test_download_truncated_is_rejected(monkeypatch, tmp_path: Path) -> None:
    # server promises 4096 bytes but the body is short -> reject, leave nothing
    _patch(monkeypatch, _FakeResp(content=b"A" * 10, headers={"Content-Length": "4096"}))
    dest = tmp_path / "f.gz"
    with pytest.raises(WorkdirError, match="Incomplete download"):
        http.download("https://x/f.gz", dest)
    assert not dest.exists()
    assert not (tmp_path / "f.gz.part").exists()


def test_download_request_error(monkeypatch, tmp_path: Path) -> None:
    _patch(monkeypatch, _FakeResp(status_exc=requests.ConnectionError("reset")))
    with pytest.raises(WorkdirError, match="Download failed"):
        http.download("https://x/f.gz", tmp_path / "f.gz")
    assert not (tmp_path / "f.gz.part").exists()


# --- MD5 manifest verification ------------------------------------------------


def test_verify_md5_manifest_ok(monkeypatch, tmp_path: Path) -> None:
    f = tmp_path / "data.tsv.gz"
    f.write_bytes(b"payload")
    import hashlib

    digest = hashlib.md5(b"payload").hexdigest()
    manifest = f"{digest}  ./data.tsv.gz\nffff  ./other.txt\n"
    monkeypatch.setattr(http, "get_text", lambda url, **k: manifest)
    assert http.verify_md5_manifest(f, "https://x/MD5SUM.txt", logger=_LOG) is True


def test_verify_md5_manifest_mismatch_raises(monkeypatch, tmp_path: Path) -> None:
    f = tmp_path / "data.tsv.gz"
    f.write_bytes(b"payload")
    manifest = "0" * 32 + "  ./data.tsv.gz\n"
    monkeypatch.setattr(http, "get_text", lambda url, **k: manifest)
    with pytest.raises(WorkdirError, match="hecksum"):
        http.verify_md5_manifest(f, "https://x/MD5SUM.txt", logger=_LOG)


def test_verify_md5_manifest_unavailable_skips(monkeypatch, tmp_path: Path) -> None:
    f = tmp_path / "data.tsv.gz"
    f.write_bytes(b"payload")

    def boom(url, **k):
        raise WorkdirError("404")

    monkeypatch.setattr(http, "get_text", boom)
    assert http.verify_md5_manifest(f, "https://x/MD5SUM.txt", logger=_LOG) is False


def test_verify_md5_manifest_unlisted_file_skips(monkeypatch, tmp_path: Path) -> None:
    f = tmp_path / "data.tsv.gz"
    f.write_bytes(b"payload")
    monkeypatch.setattr(http, "get_text", lambda url, **k: "aa" * 16 + "  ./other.txt\n")
    assert http.verify_md5_manifest(f, "https://x/MD5SUM.txt", logger=_LOG) is False


def test_verify_md5_accepts_the_right_digest_and_rejects_a_wrong_one(tmp_path) -> None:
    import hashlib

    from repgenr.core.errors import WorkdirError
    from repgenr.core.http import verify_md5

    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(b"@r1\nACGT\n+\nIIII\n")
    good = hashlib.md5(path.read_bytes()).hexdigest()
    verify_md5(path, good.upper())  # case-insensitive, no return value needed
    with pytest.raises(WorkdirError, match="Checksum mismatch"):
        verify_md5(path, "0" * 32)


def _http_error(status: int) -> requests.HTTPError:
    response = requests.Response()
    response.status_code = status
    return requests.HTTPError(f"{status} error", response=response)


@pytest.mark.parametrize("status", [404, 503])
def test_status_error_carries_the_status(monkeypatch, tmp_path: Path, status: int) -> None:
    """A caller can tell an unknown resource (404) from a failing server."""
    _patch(monkeypatch, _FakeResp(status_exc=_http_error(status)))
    with pytest.raises(http.HTTPStatusError) as info:
        http.get_json("https://x/y")
    assert info.value.status == status
    with pytest.raises(http.HTTPStatusError) as info:
        http.download("https://x/y", tmp_path / "f")
    assert info.value.status == status


def test_connection_error_is_not_a_status_error(monkeypatch) -> None:
    _patch(monkeypatch, _FakeResp(status_exc=requests.ConnectionError("refused")))
    with pytest.raises(WorkdirError) as info:
        http.get_json("https://x/y")
    assert not isinstance(info.value, http.HTTPStatusError)


def test_requests_use_a_short_connect_timeout(monkeypatch, tmp_path: Path) -> None:
    """A blocked network is reported in minutes: through an unreachable proxy
    one GTDB API request waited 481 s (six 120 s connect attempts)."""
    seen = []

    class _Recording(_FakeSession):
        def get(self, url, **kw):
            seen.append(kw["timeout"])
            return super().get(url, **kw)

    body = b"x"
    resp = _FakeResp(json_data={}, content=body, headers={"Content-Length": "1"})
    monkeypatch.setattr(http, "session", lambda: _Recording(resp))
    http.get_json("https://x/y")
    http.download("https://x/y", tmp_path / "f")
    assert seen and all(isinstance(t, tuple) and t[0] <= 30 and t[1] >= 120 for t in seen)


class _ProbeSession:
    def __init__(self, exc: Exception | None = None):
        self.exc = exc
        self.calls: list[tuple[str, dict]] = []

    def get(self, url, **kw):
        self.calls.append((url, kw))
        if self.exc is not None:
            raise self.exc
        return None

    def close(self) -> None:
        pass


def test_require_reachable_passes_when_the_host_answers(monkeypatch) -> None:
    probe = _ProbeSession()
    monkeypatch.setattr(http, "_probe_session", lambda: probe)
    http.require_reachable(http.NCBI_DATASETS_URL, what="NCBI datasets")
    [(url, kw)] = probe.calls
    assert url == http.NCBI_DATASETS_URL
    # One attempt with the shared 15 s connect timeout.
    assert kw["timeout"][0] == 15


@pytest.mark.parametrize(
    "exc",
    [requests.ConnectTimeout("timed out"), requests.ConnectionError("refused")],
)
def test_require_reachable_names_the_host_when_it_does_not_answer(monkeypatch, exc) -> None:
    """datasets on a blocked network made three attempts of about 8.5 minutes each."""
    monkeypatch.setattr(http, "_probe_session", lambda: _ProbeSession(exc))
    with pytest.raises(WorkdirError, match="api.ncbi.nlm.nih.gov") as info:
        http.require_reachable(http.NCBI_DATASETS_URL, what="NCBI datasets")
    assert info.value.exit_code == 3
    assert "HTTPS_PROXY" in str(info.value)


def test_probe_session_reads_the_proxy_environment(monkeypatch) -> None:
    """The probe goes through the same proxies as datasets and core.http."""
    monkeypatch.setenv("HTTPS_PROXY", "http://proxy.example:3128")
    monkeypatch.delenv("NO_PROXY", raising=False)
    monkeypatch.delenv("no_proxy", raising=False)
    s = _REAL_PROBE_SESSION()
    try:
        settings = s.merge_environment_settings(http.NCBI_DATASETS_URL, {}, None, None, None)
    finally:
        s.close()
    assert settings["proxies"]["https"] == "http://proxy.example:3128"
