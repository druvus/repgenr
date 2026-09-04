"""The v2 resume fingerprint: params (minus scheduling), input digests, env."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

from repgenr.cli.base import _env_fragment, _stage_fingerprint
from repgenr.core.containers import configure_container


@dataclass
class _P:
    secondary_ani: float = 0.99
    threads: int = 16
    num_processes: int = 1


def test_threads_do_not_change_fingerprint() -> None:
    base = _stage_fingerprint("dereplicate", _P(threads=16, num_processes=1), {}, {})
    more = _stage_fingerprint("dereplicate", _P(threads=64, num_processes=8), {}, {})
    assert base == more  # only scheduling differs -> same fingerprint, resume skips


def test_result_param_changes_fingerprint() -> None:
    a = _stage_fingerprint("dereplicate", _P(secondary_ani=0.99), {}, {})
    b = _stage_fingerprint("dereplicate", _P(secondary_ani=0.95), {}, {})
    assert a != b


def test_stage_name_changes_fingerprint() -> None:
    assert _stage_fingerprint("phylo", _P(), {}, {}) != _stage_fingerprint(
        "dereplicate", _P(), {}, {}
    )


def test_input_digest_changes_fingerprint() -> None:
    a = _stage_fingerprint("phylo", _P(), {"derep/representatives": "d1"}, {})
    b = _stage_fingerprint("phylo", _P(), {"derep/representatives": "d2"}, {})
    assert a != b


def test_env_changes_fingerprint() -> None:
    native = _stage_fingerprint("phylo", _P(), {}, {"container": ["none", None, False]})
    docker = _stage_fingerprint("phylo", _P(), {}, {"container": ["docker", None, False]})
    assert native != docker


def test_env_fragment_tracks_backend_platform_wave_only(tmp_path) -> None:
    try:
        configure_container(backend="docker", platform="linux/amd64", wave_enabled=True)
        with_extras = _env_fragment()
        # engine and cache_dir are plumbing, not result-affecting
        configure_container(
            backend="docker",
            platform="linux/amd64",
            wave_enabled=True,
            engine="podman",
            cache_dir=str(tmp_path),
        )
        assert _env_fragment() == with_extras
        configure_container(backend="singularity", platform="linux/amd64", wave_enabled=True)
        assert _env_fragment() != with_extras
    finally:
        configure_container("none")


def test_v2_hash_never_matches_v1_shape() -> None:
    """The fpv marker guarantees pre-redesign fingerprints can never false-match."""
    params = _P()
    v1_blob = json.dumps(
        {"stage": "dereplicate", "params": {"secondary_ani": 0.99}},
        sort_keys=True,
        default=str,
    )
    v1 = hashlib.sha256(v1_blob.encode("utf-8")).hexdigest()
    v2 = _stage_fingerprint("dereplicate", params, {}, {})
    assert v1 != v2
