"""Process and command execution helpers."""

from __future__ import annotations

import os
import shlex
import subprocess
import sys
from collections.abc import Iterator, Mapping, Sequence
from typing import BinaryIO, TextIO, TypeVar, cast

from SpliceGrapher.shared.format_utils import time_string
from SpliceGrapher.shared.logging_utils import get_logger

LOGGER = get_logger(__name__)
AttrT = TypeVar("AttrT")
SubprocessStream = int | TextIO | BinaryIO | None
CompletedProcessResult = subprocess.CompletedProcess[str] | subprocess.CompletedProcess[bytes]


def get_attribute(mapping: Mapping[str, object], key: str, default: AttrT) -> AttrT:
    """Return ``mapping[key]`` when present, else ``default``."""
    return cast(AttrT, mapping[key]) if key in mapping else default


def id_factory(prefix: str = "", initial: int = 1) -> Iterator[str]:
    """Generate unique ids using the given prefix."""
    counter = initial
    while True:
        yield f"{prefix}{counter}"
        counter += 1


def log_message(message: str, logstream: TextIO | None = None) -> None:
    """Emit a structured process message and optionally mirror it to a stream."""
    LOGGER.info("process_message", message=message.rstrip("\n"))
    if logstream is not None:
        logstream.write(message)


def run_command(
    command: str | Sequence[str],
    *,
    shell: bool = False,
    check: bool = False,
    stdout: SubprocessStream = None,
    stderr: SubprocessStream = None,
    text: bool = False,
) -> CompletedProcessResult:
    """Run a command with explicit shell behavior."""
    command_args: str | list[str]
    if isinstance(command, str):
        command_args = command if shell else shlex.split(command)
    else:
        command_args = list(command)

    return subprocess.run(
        command_args,
        shell=shell,
        check=check,
        stdout=stdout,
        stderr=stderr,
        text=text,
    )


def run_logged_command(
    command: str,
    *,
    logstream: TextIO | None = None,
    debug: bool = False,
    exit_on_error: bool = True,
    stderr: SubprocessStream = None,
    stdout: SubprocessStream = None,
) -> None:
    """Log and run a shell command with SGN defaults."""
    LOGGER.info("command_started", command=command, debug=debug)
    message = "    " + time_string(f"{command}\n")
    if logstream is not None:
        logstream.write(message)

    retcode = 0
    if not debug:
        stderr_stream = stderr if stderr is not None else subprocess.DEVNULL
        stdout_stream = stdout if stdout is not None else subprocess.DEVNULL
        completed = run_command(
            command,
            shell=True,
            check=False,
            stderr=stderr_stream,
            stdout=stdout_stream,
        )
        retcode = completed.returncode

    if exit_on_error and retcode != 0:
        LOGGER.error("command_failed", command=command, return_code=retcode)
        code_type = "signal" if retcode < 0 else "code"
        raise RuntimeError(f"Error running command: returned {retcode} {code_type}\n{command}")


def write_startup_message() -> None:
    """Emit the standardized startup event for scripts."""
    script_name = os.path.basename(sys.argv[0])
    LOGGER.info("startup", script=script_name)


__all__ = [
    "CompletedProcessResult",
    "SubprocessStream",
    "get_attribute",
    "id_factory",
    "log_message",
    "run_command",
    "run_logged_command",
    "write_startup_message",
]
