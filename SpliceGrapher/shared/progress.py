"""Progress and random-iterator helpers extracted from shared.utils."""

from __future__ import annotations

import random
import sys
from collections.abc import Iterator, Sequence
from typing import Generic, TypeVar

from tqdm.auto import tqdm

T = TypeVar("T")


def _is_tty(stream: object) -> bool:
    """Return True when a stream is interactive."""
    try:
        isatty = getattr(stream, "isatty", None)
        return bool(isatty and isatty())
    except Exception:
        return False


class ProgressIndicator:
    """Compatibility wrapper around a tqdm progress bar."""

    def __init__(
        self,
        increment: int,
        description: str = "",
        verbose: bool = True,
    ) -> None:
        self.limit = max(1, int(increment))
        self.descr = description
        self.verbose = verbose
        self.ctr = 0
        self._enabled = self.verbose and _is_tty(sys.stderr)
        self._bar = tqdm(
            total=None,
            desc=self.descr or None,
            disable=not self._enabled,
            leave=False,
            unit="records",
            miniters=self.limit,
            dynamic_ncols=True,
            file=sys.stderr,
        )

    def count(self) -> int:
        """Returns the current count."""
        return self.ctr

    def finish(self) -> None:
        """Finishes progress output."""
        self._bar.close()

    def reset(self) -> None:
        """Resets the indicator to be used again."""
        self.finish()
        self.ctr = 0
        self._bar = tqdm(
            total=None,
            desc=self.descr or None,
            disable=not self._enabled,
            leave=False,
            unit="records",
            miniters=self.limit,
            dynamic_ncols=True,
            file=sys.stderr,
        )

    def update(self) -> None:
        """Updates the indicator."""
        self.ctr += 1
        if not self.verbose:
            return
        self._bar.update(1)


class RandomListIterator(Generic[T], Iterator[T]):
    """
    Returns items at random, with replacement, from a list of values.
    For lists with more than ~1000 elements, this is about 30% more
    efficient than using random.sample() repeatedly.
    """

    def __init__(self, values: Sequence[T], *, seed: int | None = None) -> None:
        self.values = values
        self.rand = random.Random()
        self.limit = len(self.values) - 1
        if seed is not None:
            self.rand.seed(seed)

    def __iter__(self) -> RandomListIterator[T]:
        return self

    def __next__(self) -> T:
        """Iterator implementation that returns a random value, with
        replacement, from a list."""
        i = self.rand.randint(0, self.limit)
        return self.values[i]

    # Python 2 compatibility alias
    next = __next__
