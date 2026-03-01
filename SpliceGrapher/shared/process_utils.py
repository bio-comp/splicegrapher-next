"""Process and command execution helpers extracted from shared.utils."""

import os
import subprocess
import sys

from SpliceGrapher.shared.format_utils import time_string


def getAttribute(key, default, **args):
    """Returns the value for the given key in the arguments dict, if found;
    otherwise returns default value."""
    return default if key not in args else args[key]


def idFactory(pfx="", initial=1):
    """Generates unique ids using the given prefix.  For example,
    idFactory('ev_') will generate 'ev_1', 'ev_2', ..."""
    prefix = pfx if pfx else ""
    counter = initial
    while 1:
        yield "%s%d" % (prefix, counter)
        counter += 1


def logMessage(s, logstream=None):
    """Allows log messages to be output both to stderr and a log file."""
    sys.stderr.write(s)
    if logstream:
        logstream.write(s)


def runCommand(s, **args):
    """Announces a command runs it.."""
    logstream = getAttribute("logstream", None, **args)
    debug = getAttribute("debug", False, **args)
    exitOnError = getAttribute("exitOnError", True, **args)
    stderr = getAttribute("stderr", None, **args)
    stdout = getAttribute("stdout", None, **args)
    message = "    " + time_string(f"{s}\n")
    sys.stderr.write(message)
    if logstream:
        logstream.write(message)

    retcode = 0
    if not debug:
        stderr_stream = stderr if stderr is not None else subprocess.DEVNULL
        stdout_stream = stdout if stdout is not None else subprocess.DEVNULL
        retcode = subprocess.call(
            s,
            shell=True,
            stderr=stderr_stream,
            stdout=stdout_stream,
        )

    if exitOnError and retcode < 0:
        raise Exception("Error running command: returned %d signal\n%s" % (retcode, s))


def writeStartupMessage():
    """Standardized startup message for all scripts."""
    base = os.path.basename(sys.argv[0])
    sys.stderr.write(time_string(f"{base} Started\n"))
