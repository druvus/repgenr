"""Exception hierarchy for RepGenR.

Stages and plugins raise these instead of calling ``sys.exit``. The CLI layer
(:mod:`repgenr.cli.main`) catches :class:`RepGenRError`, logs it cleanly and
sets a non-zero exit code, so library callers and tests can handle failures
without a process exit.
"""

from __future__ import annotations


class RepGenRError(Exception):
    """Base class for all expected, user-facing RepGenR failures.

    ``exit_code`` is the process status the CLI exits with for the class, so
    a caller can tell a usage error from a missing tool or a failed run
    without parsing the log. 1 is reserved for unexpected exceptions.
    """

    exit_code: int = 1


class UserInputError(RepGenRError):
    """Invalid or missing user input (bad arguments, missing target, etc.)."""

    exit_code = 2


class WorkdirError(RepGenRError):
    """A working directory is missing required files or is in a bad state."""

    exit_code = 3


class MissingBinaryError(RepGenRError):
    """A required external tool was not found on PATH or failed its version check."""

    exit_code = 4


class ToolExecutionError(RepGenRError):
    """An external tool exited with a non-zero status.

    The CLI exits 6, or with the tool's own status under
    ``REPGENR_PROPAGATE_TOOL_EXIT=1`` (see ``cli.base``).
    """

    exit_code = 6

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

    exit_code = 5
