"""Formatting and parsing helpers extracted from shared.utils."""

import time


def commaFormat(d):
    """Formats integer values using commas.  For example, 123456789 becomes '123,456,789'"""
    return f"{int(d):,}"


def dictString(valDict, delim=","):
    """Returns a simple string representation of a list of values."""
    return delim.join(["%s -> %s" % (x, valDict[x]) for x in valDict])


def listString(vals, delim=","):
    """Returns a simple string representation of a list of values."""
    return delim.join([str(x) for x in vals])


def substringAfter(s, tag):
    """Returns the substring after the given tag, if found; None otherwise."""
    pos = s.find(tag)
    if pos >= 0:
        return s[pos + len(tag) :]


def substringBefore(s, tag):
    """Returns the substring before the given tag, if found; None otherwise."""
    pos = s.find(tag)
    if pos >= 0:
        return s[:pos]


def substringBetween(s, tag1, tag2):
    """Returns the substring between the two tags, if found; None otherwise."""
    p1 = s.find(tag1)
    if p1 < 0:
        return None
    p1 += len(tag1)
    p2 = s.find(tag2, p1)
    if p2 > 0:
        return s[p1:p2]


def timestamp(format_string="%Y%m%d%H%M%S"):
    """Returns a timestamp unique to the current second."""
    return time.strftime(format_string, time.localtime())


def timeString(s, format_string="%X", LF=False):
    """Returns the input string with user-readable a timestamp prefix."""
    timestamp = time.strftime(format_string, time.localtime())
    result = "%s %s" % (timestamp, s)
    if LF:
        result += "\n"
    return result


def to_numeric(s):
    """Attempts to return the given value as a float or an int, else returns the original string."""
    try:
        return int(s)
    except ValueError:
        # (fails on numbers with decimals like '12.34')
        pass

    try:
        return float(s)
    except ValueError:
        pass

    return s
