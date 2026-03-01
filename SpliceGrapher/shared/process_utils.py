"""Process and command execution helpers extracted from shared.utils."""

import os
import shlex
import subprocess
import sys
from collections.abc import Sequence
from typing import Any

from SpliceGrapher.shared.format_utils import time_string
from SpliceGrapher.shared.logging_utils import get_logger

LOGGER = get_logger(__name__)


def getAttribute(key, default, **args):
    """Returns the value for the given key in the arguments dict, if found;
    otherwise returns default value."""
    return default if key not in args else args[key]


def idFactory(pfx="", initial=1):
    """Generates unique ids using the given prefix.  For example,
    idFactory('ev_') will generate 'ev_1', 'ev_2', ..."""
    prefix = pfx if pfx else ""
    counter = initial
    while True:
        yield "%s%d" % (prefix, counter)
        counter += 1


def logMessage(s, logstream=None):
    """Allows log messages to be output both to stderr and a log file."""
    LOGGER.info("process_message", message=s.rstrip("\n"))
    sys.stderr.write(s)
    if logstream:
        logstream.write(s)


def run_command(
    command: str | Sequence[str],
    *,
    shell: bool = False,
    check: bool = False,
    stdout: Any = None,
    stderr: Any = None,
    text: bool = False,
) -> subprocess.CompletedProcess[Any]:
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


def runCommand(s, **args):
    """Announces a command runs it.."""
    logstream = getAttribute("logstream", None, **args)
    debug = getAttribute("debug", False, **args)
    exitOnError = getAttribute("exitOnError", True, **args)
    stderr = getAttribute("stderr", None, **args)
    stdout = getAttribute("stdout", None, **args)
    LOGGER.info("command_started", command=s, debug=debug)
    message = "    " + time_string(f"{s}\n")
    sys.stderr.write(message)
    if logstream:
        logstream.write(message)

    retcode = 0
    if not debug:
        stderr_stream = stderr if stderr is not None else subprocess.DEVNULL
        stdout_stream = stdout if stdout is not None else subprocess.DEVNULL
        completed = run_command(
            s,
            shell=True,
            check=False,
            stderr=stderr_stream,
            stdout=stdout_stream,
        )
        retcode = completed.returncode

    if exitOnError and retcode != 0:
        LOGGER.error("command_failed", command=s, return_code=retcode)
        code_type = "signal" if retcode < 0 else "code"
        raise Exception("Error running command: returned %d %s\n%s" % (retcode, code_type, s))


def writeStartupMessage():
    """Standardized startup message for all scripts."""
    base = os.path.basename(sys.argv[0])
    LOGGER.info("startup", script=base)
    sys.stderr.write(time_string(f"{base} Started\n"))
