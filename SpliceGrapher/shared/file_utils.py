"""Filesystem and file-format helpers extracted from shared.utils."""

import gzip
import os
from glob import glob


def ezopen(fileName):
    """Allows clients to open files without regard for whether they're gzipped."""
    if not (os.path.exists(fileName) and os.path.isfile(fileName)):
        raise ValueError("file does not exist at %s" % fileName)

    fileHandle = gzip.GzipFile(fileName)
    try:
        fileHandle.readline()
        return gzip.GzipFile(fileName)
    except Exception:
        return open(fileName)
    finally:
        fileHandle.close()


def fileLen(path):
    """Simple function to get an exact file length."""
    line_count = 0
    with open(path) as instream:
        for line_count, _line in enumerate(instream, start=1):
            pass
    return line_count


def filePrefix(f):
    """Returns the filename prefix for a file.  For example:
    /my/dir/myfile.ext --> myfile"""
    _head, tail = os.path.split(f)
    prefix, _suffix = os.path.splitext(tail)
    return prefix


def findFile(name, path, delim=":"):
    """Finds the first instance of a file name in the given path string."""
    paths = path.split(delim)
    for p in paths:
        filePath = os.path.join(p, name)
        if os.path.exists(filePath) and os.path.isfile(filePath):
            return filePath


def makeGraphListFile(spliceGraphDir):
    """Given a path to a top-level directory of splice graphs, returns
    a file that contains the paths to all graphs under that directory.
    Follows the standard SpliceGrapher directory structure:
        top-level-dir/chromosome-dir/splice-graph-file"""
    subdirs = os.path.join(spliceGraphDir, "*")
    target = os.path.join(subdirs, "*.gff")
    graphList = glob(target)
    if not graphList:
        raise ValueError("No splice graphs found in %s\n" % spliceGraphDir)

    graphList.sort()
    result = "%s.lis" % spliceGraphDir
    with open(result, "w") as graphStream:
        graphStream.write("\n".join(graphList))
    return result


def validateDir(path):
    """Standard method for validating directory paths."""
    validateFile(path)
    if not os.path.isdir(path):
        raise Exception("'%s' is not a directory; exiting." % path)


def validateFile(path):
    """Standard method for validating file paths."""
    if not path:
        raise Exception("'%s' is not a valid file path; exiting." % path)

    if not os.path.exists(path):
        raise Exception("File '%s' not found; exiting." % path)
