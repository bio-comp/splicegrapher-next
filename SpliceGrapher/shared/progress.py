"""Progress and random-iterator helpers extracted from shared.utils."""

import random
import sys

from SpliceGrapher.shared.format_utils import commaFormat


class ProgressIndicator:
    """A simple progress indicator."""

    def __init__(self, increment, description="", verbose=True):
        self.limit = increment
        self.barlim = int(self.limit / 10)
        self.dotlim = int(self.barlim / 5)
        self.descr = description
        self.started = False
        self.verbose = verbose
        self.ctr = 0

    def count(self):
        """Returns the current count."""
        return self.ctr

    def finish(self):
        """Finishes the progress output by appending a newline,
        if anything has been written."""
        if self.started:
            sys.stderr.write("\n")
        self.started = False

    def reset(self):
        """Resets the indicator to be used again."""
        self.ctr = 0
        self.finish()

    def update(self):
        """Updates the indicator."""
        self.ctr += 1
        if not self.verbose:
            return
        if self.ctr % self.limit == 0:
            sys.stderr.write("%s %s\n" % (commaFormat(self.ctr), self.descr))
        elif self.ctr % self.barlim == 0:
            sys.stderr.write("|")
        elif self.ctr % self.dotlim == 0:
            self.started = True
            sys.stderr.write(".")


class RandomListIterator:
    """
    Returns items at random, with replacement, from a list of values.
    For lists with more than ~1000 elements, this is about 30% more
    efficient than using random.sample() repeatedly.
    """

    def __init__(self, values, **args):
        self.values = values
        self.rand = random.Random()
        self.limit = len(self.values) - 1
        if "seed" in args:
            self.rand.seed = args["seed"]

    def __iter__(self):
        return self

    def __next__(self):
        """Iterator implementation that returns a random value, with
        replacement, from a list."""
        i = self.rand.randint(0, self.limit)
        return self.values[i]

    # Python 2 compatibility alias
    next = __next__
