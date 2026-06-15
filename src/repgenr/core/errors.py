"""Exception hierarchy for RepGenR.

Stages and plugins raise these instead of calling ``sys.exit``. The CLI layer
(:mod:`repgenr.cli.main`) catches :class:`RepGenRError`, logs it cleanly and
sets a non-zero exit code, so library callers and tests can handle failures
without a process exit.
"""

from __future__ import annotations


class RepGenRError(Exception):
    """Base class for all expected, user-facing RepGenR failures."""


class UserInputError(RepGenRError):
    """Invalid or missing user input (bad arguments, missing target, etc.)."""


class WorkdirError(RepGenRError):
    """A working directory is missing required files or is in a bad state."""


class MissingBinaryError(RepGenRError):
    """A required external tool was not found on PATH or failed its version check."""


class ToolExecutionError(RepGenRError):
    """An external tool exited with a non-zero status."""

    def __init__(self, command: list[str], returncode: int, output: str | None = None):
        self.command = command
        self.returncode = returncode
        self.output = output
        rendered = " ".join(command)
        msg = f"command failed (exit {returncode}): {rendered}"
        if output:
            msg += f"\n--- output tail ---\n{output}"
        super().__init__(msg)


class PluginError(RepGenRError):
    """A requested plugin/tool could not be found or loaded."""
