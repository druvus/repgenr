"""A failure inside the output read loop must not leave the tool running."""

from __future__ import annotations

import logging
import sys
import time
from pathlib import Path

import pytest

from repgenr.core import process


class _RaisingLogger(logging.Logger):
    """Raises on the first streamed line, standing in for any consumer failure."""

    def debug(self, msg, *args, **kwargs):  # noqa: ANN001
        if args and args[-1] == "started":
            raise RuntimeError("consumer failed")
        super().debug(msg, *args, **kwargs)


def test_read_loop_failure_kills_the_child(tmp_path: Path) -> None:
    marker = tmp_path / "child_finished"
    script = (
        "import sys, time, pathlib; print('started', flush=True); "
        f"time.sleep(2); pathlib.Path({str(marker)!r}).write_text('done')"
    )
    logger = _RaisingLogger("test.child.cleanup")
    logger.setLevel(logging.DEBUG)
    started = time.monotonic()
    with pytest.raises(RuntimeError, match="consumer failed"):
        process.run([sys.executable, "-c", script], logger=logger)
    # The exception propagated promptly and the orphan never got to finish.
    assert time.monotonic() - started < 1.5
    time.sleep(2.5)
    assert not marker.exists(), "child kept running after the read loop failed"


def test_sigterm_stops_the_running_tool(tmp_path: Path) -> None:
    """A repgenr process ended by SIGTERM stops its tool and removes the partial
    output (live: FastTree kept running after its phylo run was terminated)."""
    import os
    import signal
    import subprocess

    started = tmp_path / "started"
    finished = tmp_path / "finished"
    tool = (
        "import pathlib, time; "
        f"pathlib.Path({str(started)!r}).write_text('x'); "
        f"time.sleep(3); pathlib.Path({str(finished)!r}).write_text('x')"
    )
    out = tmp_path / "tree.nwk"
    driver = (
        "import logging, sys\n"
        "from repgenr.core import process\n"
        "process.install_termination_handler()\n"
        f"process.run([sys.executable, '-c', {tool!r}], logger=logging.getLogger('t'),"
        f" stdout_path={str(out)!r})\n"
    )
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(sys.path)}
    proc = subprocess.Popen(
        [sys.executable, "-c", driver], stderr=subprocess.PIPE, env=env, text=True
    )
    deadline = time.monotonic() + 10
    while not started.exists():
        assert time.monotonic() < deadline, "the tool never started"
        time.sleep(0.05)
    proc.send_signal(signal.SIGTERM)
    _, stderr = proc.communicate(timeout=10)
    assert proc.returncode == 128 + signal.SIGTERM
    assert "stopped 1 running tool" in stderr
    time.sleep(3.5)
    assert not finished.exists(), "the tool kept running after repgenr was terminated"
    assert not out.exists() and not (tmp_path / "tree.nwk.part").exists()


def test_sigterm_during_parallel_map_starts_no_queued_tool(tmp_path: Path) -> None:
    """Queued pool items do not launch their tools once SIGTERM arrived (review
    of #223: 8 tasks on 2 workers all started before repgenr exited)."""
    import os
    import signal
    import subprocess

    marks = tmp_path / "marks"
    marks.mkdir()
    tool = (
        "import pathlib, sys, time; "
        f"pathlib.Path({str(marks)!r}, sys.argv[1]).write_text('x'); time.sleep(3)"
    )
    driver = (
        "import logging, sys\n"
        "from repgenr.core import process\n"
        "from repgenr.core.executors import parallel_map\n"
        "process.install_termination_handler()\n"
        "log = logging.getLogger('t')\n"
        f"parallel_map(lambda i: process.run([sys.executable, '-c', {tool!r}, str(i)],"
        " logger=log), range(8), 2)\n"
    )
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(sys.path)}
    proc = subprocess.Popen(
        [sys.executable, "-c", driver], stderr=subprocess.PIPE, env=env, text=True
    )
    deadline = time.monotonic() + 10
    while len(list(marks.iterdir())) < 2:
        assert time.monotonic() < deadline, "the first tools never started"
        time.sleep(0.05)
    sent = time.monotonic()
    proc.send_signal(signal.SIGTERM)
    proc.communicate(timeout=20)
    assert proc.returncode == 128 + signal.SIGTERM
    assert time.monotonic() - sent < 2.5, "repgenr waited for queued work"
    time.sleep(1)
    assert len(list(marks.iterdir())) == 2, sorted(p.name for p in marks.iterdir())


# A tool that starts a helper: the shell runs ``sleep`` in the background and
# waits for it, as run_gubbins.py runs RAxML or skDER runs skani. Each process
# records its PID so the test can check that both are gone.
def _helper_tool(tmp_path: Path, *, trap_term: bool = False) -> list[str]:
    child = tmp_path / "child.pid"
    grandchild = tmp_path / "grandchild.pid"
    trap = "trap '' TERM INT; " if trap_term else ""
    script = (
        f"{trap}sleep 60 & echo $! > {grandchild}.tmp; mv {grandchild}.tmp {grandchild}; "
        f"echo $$ > {child}; wait"
    )
    return ["/bin/sh", "-c", script]


def _pids(tmp_path: Path) -> list[int]:
    return [int((tmp_path / f"{n}.pid").read_text()) for n in ("child", "grandchild")]


def _wait_started(tmp_path: Path, proc) -> None:  # noqa: ANN001
    deadline = time.monotonic() + 10
    while not ((tmp_path / "child.pid").exists() and (tmp_path / "grandchild.pid").exists()):
        assert proc.poll() is None, proc.communicate()
        assert time.monotonic() < deadline, "the tool never started"
        time.sleep(0.05)


def _alive(pid: int) -> bool:
    import os

    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    # A zombie still answers kill(0); it is not running.
    import subprocess

    state = subprocess.run(
        ["ps", "-o", "stat=", "-p", str(pid)], capture_output=True, text=True
    ).stdout.strip()
    return bool(state) and not state.startswith("Z")


def _assert_gone_within(pids: list[int], seconds: float) -> None:
    deadline = time.monotonic() + seconds
    while any(_alive(p) for p in pids):
        if time.monotonic() > deadline:
            import os
            import signal

            survivors = [p for p in pids if _alive(p)]
            for p in survivors:
                os.kill(p, signal.SIGKILL)
            raise AssertionError(f"still running after {seconds} s: {survivors}")
        time.sleep(0.05)


def _driver(*commands: list[str], grace: float | None = None) -> str:
    lines = [
        "import logging, sys",
        "from repgenr.core import process",
        "from repgenr.core.executors import parallel_map",
        "process.install_termination_handler()",
        "log = logging.getLogger('t')",
    ]
    if grace is not None:
        lines.append(f"process.STOP_GRACE_SECONDS = {grace}")
    if len(commands) > 1:
        # The tools run on pool threads while the main thread waits.
        call = "parallel_map(lambda c: process.run(c, logger=log)"
        lines.append(f"{call}, {list(commands)!r}, {len(commands)})")
    else:
        lines.append(f"process.run({commands[0]!r}, logger=log)")
    return "\n".join(lines) + "\n"


def _start_driver(driver: str):  # noqa: ANN202
    import os
    import subprocess

    env = {**os.environ, "PYTHONPATH": os.pathsep.join(sys.path)}
    return subprocess.Popen(
        [sys.executable, "-c", driver], stderr=subprocess.PIPE, env=env, text=True
    )


def _finish(proc, timeout: float = 10) -> str:  # noqa: ANN001
    """Wait for the driver; kill it on a hang so no test leaves it behind."""
    import subprocess

    try:
        return proc.communicate(timeout=timeout)[1]
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.communicate()
        raise


@pytest.mark.parametrize("signame", ["SIGTERM", "SIGHUP", "SIGINT"])
def test_signal_stops_the_tool_and_its_helpers(tmp_path: Path, signame: str) -> None:
    """The tool runs in its own process group; a signal to repgenr stops the
    tool and the helpers it started (RAxML under run_gubbins.py, live), and
    Ctrl-C still ends repgenr although the tool no longer shares its terminal."""
    import signal

    sig = getattr(signal, signame)
    proc = _start_driver(_driver(_helper_tool(tmp_path)))
    _wait_started(tmp_path, proc)
    pids = _pids(tmp_path)
    sent = time.monotonic()
    proc.send_signal(sig)
    stderr = _finish(proc)
    elapsed = time.monotonic() - sent
    _assert_gone_within(pids, 2)
    assert elapsed < 2, "repgenr did not end promptly"
    if sig == signal.SIGINT:
        assert proc.returncode != 0
        assert "KeyboardInterrupt" in stderr
    else:
        assert proc.returncode == 128 + sig, stderr


def test_tool_in_own_session_ignores_terminal_signals(tmp_path: Path) -> None:
    """The tool leads its own session, so a Ctrl-C sent to repgenr's process
    group (what a terminal does) reaches repgenr only, which then stops it."""
    import os

    proc = _start_driver(_driver(_helper_tool(tmp_path)))
    _wait_started(tmp_path, proc)
    child, grandchild = _pids(tmp_path)
    try:
        assert os.getsid(child) == child
        assert os.getpgid(grandchild) == child
        assert os.getsid(child) != os.getsid(proc.pid)
    finally:
        proc.terminate()
        _finish(proc)
    _assert_gone_within([child, grandchild], 2)


@pytest.mark.parametrize("signame", ["SIGTERM", "SIGINT"])
def test_signal_stops_helpers_of_a_tool_on_a_pool_thread(tmp_path: Path, signame: str) -> None:
    """With the tool on a parallel_map thread the main thread is not inside
    run(); the handler still stops the tool group (Ctrl-C used to reach the
    tool through the terminal; it no longer does)."""
    import signal

    dirs = [tmp_path / "a", tmp_path / "b"]
    for d in dirs:
        d.mkdir()
    proc = _start_driver(_driver(*(_helper_tool(d) for d in dirs)))
    for d in dirs:
        _wait_started(d, proc)
    pids = [pid for d in dirs for pid in _pids(d)]
    sent = time.monotonic()
    proc.send_signal(getattr(signal, signame))
    _finish(proc)
    elapsed = time.monotonic() - sent
    _assert_gone_within(pids, 2)
    assert elapsed < 2, "repgenr did not end promptly"


def test_helpers_that_ignore_sigterm_are_killed_after_the_grace_period(tmp_path: Path) -> None:
    """A helper that ignores SIGTERM is killed (SIGKILL) once the grace period ends."""
    import signal

    proc = _start_driver(_driver(_helper_tool(tmp_path, trap_term=True), grace=0.5))
    _wait_started(tmp_path, proc)
    pids = _pids(tmp_path)
    proc.send_signal(signal.SIGTERM)
    _finish(proc)
    _assert_gone_within(pids, 2)
    assert proc.returncode == 128 + signal.SIGTERM


def test_read_loop_failure_kills_the_helpers(tmp_path: Path) -> None:
    """Without a signal handler (library use), a failure inside run() stops
    the whole tool group, not only the tool."""
    logger = _RaisingLogger("test.child.helpers")
    logger.setLevel(logging.DEBUG)
    command = _helper_tool(tmp_path)
    command[-1] = command[-1].replace("wait", "echo started; wait")
    with pytest.raises(RuntimeError, match="consumer failed"):
        process.run(command, logger=logger)
    _assert_gone_within(_pids(tmp_path), 2)


def _state(pid: int) -> str:
    import subprocess

    return subprocess.run(
        ["ps", "-o", "stat=", "-p", str(pid)], capture_output=True, text=True
    ).stdout.strip()


def _wait_state(pid: int, stopped: bool) -> None:
    deadline = time.monotonic() + 3
    while _state(pid).startswith("T") != stopped:
        assert time.monotonic() < deadline, f"pid {pid} state {_state(pid)!r}"
        time.sleep(0.05)


def test_ctrl_z_suspends_and_resumes_the_tool_group(tmp_path: Path) -> None:
    """Ctrl-Z (SIGTSTP) no longer reaches a tool in its own session; repgenr
    suspends the tool group with itself and resumes it on SIGCONT."""
    import signal

    proc = _start_driver(_driver(_helper_tool(tmp_path)))
    _wait_started(tmp_path, proc)
    pids = _pids(tmp_path)
    try:
        proc.send_signal(signal.SIGTSTP)
        _wait_state(proc.pid, stopped=True)
        for pid in pids:
            _wait_state(pid, stopped=True)
        proc.send_signal(signal.SIGCONT)
        _wait_state(proc.pid, stopped=False)
        for pid in pids:
            _wait_state(pid, stopped=False)
    finally:
        proc.send_signal(signal.SIGCONT)
        proc.terminate()
        _finish(proc)
    _assert_gone_within(pids, 2)
