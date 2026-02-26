"""Collection and searching helpers extracted from shared.utils."""

from SpliceGrapher.shared.process_utils import getAttribute


def asList(value, delim=","):
    """Generic method for creating a list from a given value.
    The value may be a string of delimiter-separated values,
    a list or a set."""
    if isinstance(value, str):
        return value.split(delim)
    elif isinstance(value, list):
        return value
    elif isinstance(value, set):
        return list(value)
    else:
        raise ValueError("Expected a string, a list or a set; received %s" % type(value))


def asSet(value, delim=","):
    """Generic method for creating a set from a given value.
    The value may be a string of delimiter-separated values,
    a list or a set."""
    if isinstance(value, str):
        return set(value.split(delim))
    elif isinstance(value, list):
        return set(value)
    elif isinstance(value, set):
        return value
    else:
        raise ValueError("Expected a string, a list or a set; received %s" % type(value))


def bsearch(X, target, getValue=lambda a: a, **args):
    """
    Generic binary search returns the index of the value in X that
    is closest to the target.  'getValue' provides access to the correct
    attribute of each element in X.  Assumes X is sorted such that
    getValue(a) < getValue(b) iff a < b.
    """
    if not X:
        raise ValueError("Cannot perform binary search on an empty list")
    lo = getAttribute("low", 0, **args)
    hi = getAttribute("high", len(X) - 1, **args)
    mid = lo
    while lo < hi:
        mid = (lo + hi) // 2
        midval = getValue(X[mid])
        if midval < target:
            lo = mid + 1
        elif midval > target:
            hi = mid - 1
        else:
            return mid
    return lo
