"""
Module that provides control over the sys.stderr and sys.out streams.
This is particularly useful for Python-wrapped C code that may not
provide control over these streams.
"""

import os
import sys

SAVE_STDOUT = os.dup(sys.stdout.fileno())
SAVE_STDERR = os.dup(sys.stderr.fileno())
NULL_STREAM = open(os.devnull, "a")


def hideAll(verbose=False):
    """Hides output to both stdout and stderr."""
    hideStderr(verbose=verbose)
    hideStdout(verbose=verbose)


def hideStderr(verbose=False):
    """Hides output to sys.stderr until the next call to showStderr()."""
    if verbose:
        sys.stderr.write("Hiding stderr\n")
    os.dup2(NULL_STREAM.fileno(), 2)


def hideStdout(verbose=False):
    """Hides output to sys.stdout until the next call to showStdout()."""
    if verbose:
        sys.stdout.write("Hiding stdout\n")
    os.dup2(NULL_STREAM.fileno(), 1)


def showAll(verbose=False):
    """Reinstates output to both stdout and stderr."""
    showStderr(verbose=verbose)
    showStdout(verbose=verbose)


def showStderr(verbose=False):
    """Reinstates output to sys.stderr."""
    sys.stderr.flush()
    os.dup2(SAVE_STDERR, 2)
    if verbose:
        sys.stderr.write("Reinstated stderr\n")


def showStdout(verbose=False):
    """Reinstates output to sys.stdout."""
    sys.stdout.flush()
    os.dup2(SAVE_STDOUT, 1)
    if verbose:
        sys.stdout.write("Reinstated stdout\n")
