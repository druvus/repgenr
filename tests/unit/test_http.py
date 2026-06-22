"""Tests for the shared HTTP client (no network).

The session is replaced with a fake so the JSON/text/download helpers and their
error surfaces are exercised deterministically: status errors and unparseable
bodies become WorkdirError, and a short download is rejected with no leftover
file.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import requests

from repgenr.core import http
from repgenr.core.errors import WorkdirError


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
