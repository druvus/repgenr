"""Containers of a stopped repgenr are stopped too.

Three busybox containers that ignore SIGTERM run on parallel_map worker
threads. After one signal (the escalation timer kills the docker clients) or
two (repgenr ends at once and only the final stop hook runs), no container
labelled with repgenr's PID may remain. Skipped unless the Docker daemon
answers and the busybox image is present or can be pulled.
"""

from __future__ import annotations

import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

pytestmark = [pytest.mark.live, pytest.mark.container]

IMAGE = "busybox:1.36"


@pytest.fixture(scope="module", autouse=True)
def docker_preflight() -> None:
    if shutil.which("docker") is None:
        pytest.skip("docker binary not found")
    if subprocess.run(["docker", "info"], capture_output=True).returncode != 0:
        pytest.skip("docker daemon not running")
    if subprocess.run(["docker", "image", "inspect", IMAGE], capture_output=True).returncode:
        if subprocess.run(["docker", "pull", IMAGE], capture_output=True).returncode:
            pytest.skip(f"{IMAGE} not available")


_DRIVER = f"""
import logging, sys
from repgenr.core import containers, executors, process
from repgenr.core.plugins import ToolCapabilities
logging.basicConfig(level=logging.INFO, stream=sys.stderr)
log = logging.getLogger("signals")
containers.configure_container("docker")
process.install_termination_handler()
caps = ToolCapabilities(name="busy", container={IMAGE!r})
def task(i):
    containers.run_tool(caps, ["sh", "-c", "trap '' TERM; sleep 300"], logger=log, cwd=sys.argv[1])
executors.parallel_map(task, [0, 1, 2], 3)
"""


def _labelled(pid: int) -> list[str]:
    out = subprocess.run(
        ["docker", "ps", "-q", "--filter", f"label=repgenr.pid={pid}"],
        capture_output=True,
        text=True,
    ).stdout
    return out.split()


@pytest.mark.parametrize("signals", [1, 2])
def test_worker_containers_are_stopped_with_repgenr(tmp_path: Path, signals: int) -> None:
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(sys.path)}
    proc = subprocess.Popen(
        [sys.executable, "-c", _DRIVER, str(tmp_path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        env=env,
        text=True,
    )
    try:
        deadline = time.monotonic() + 60
        while len(_labelled(proc.pid)) < 3:
            assert proc.poll() is None, proc.communicate()[1]
            assert time.monotonic() < deadline, "the containers never started"
            time.sleep(0.5)
        proc.send_signal(signal.SIGTERM)
        if signals == 2:
            time.sleep(1)
            proc.send_signal(signal.SIGTERM)
        proc.communicate(timeout=30)
        deadline = time.monotonic() + 15
        while _labelled(proc.pid):
            assert time.monotonic() < deadline, f"left running: {_labelled(proc.pid)}"
            time.sleep(0.5)
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.communicate()
        left = _labelled(proc.pid)
        if left:
            subprocess.run(["docker", "stop", "--time", "0", *left], capture_output=True)
