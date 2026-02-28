"""A parser for FASTA files."""

from __future__ import annotations

import os
import random
import sys
from os import PathLike
from pathlib import Path
from typing import IO, Iterator

from pyfaidx import Fasta

from SpliceGrapher.shared.file_utils import ezopen


class MalformedInput(Exception):
    """Exception raised when the input file does not look like a FASTA file."""

    def __init__(self, value: str):
        self.param = value
        super().__init__(value)

    def __repr__(self) -> str:
        return repr(self.param)


class FastaRecord:
    """A FASTA record."""

    def __init__(self, header: str, sequence: str):
        """Create a record with the given header and sequence."""
        self.header = header
        self.sequence = sequence

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, FastaRecord):
            return False
        return self.header == other.header and self.sequence == other.sequence

    def __str__(self) -> str:
        return ">" + self.header + "\n" + self.sequence + "\n"


def _to_text(line: str | bytes) -> str:
    if isinstance(line, bytes):
        return line.decode("utf-8")
    return line


def _fasta_itr_from_file(file_handle: IO[str] | IO[bytes]) -> Iterator[FastaRecord]:
    """Provide an iteration through FASTA records in an open file handle."""
    first_line = file_handle.readline()
    if not first_line:
        raise MalformedInput("FASTA file is empty")

    header = _to_text(first_line).strip()
    if not header.startswith(">"):
        raise MalformedInput(f"First header in FASTA file must start with '>' (found {header})")
    header = header[1:]

    sequence_chunks: list[str] = []
    line_number = 1
    for raw_line in file_handle:
        line_number += 1
        line = _to_text(raw_line).strip()
        if not line:
            raise MalformedInput(f"Blank line in FASTA file (line {line_number})")

        if line.startswith(">"):
            yield FastaRecord(header, "".join(sequence_chunks))
            header = line[1:]
            sequence_chunks = []
            continue
        sequence_chunks.append(line)

    yield FastaRecord(header, "".join(sequence_chunks))


def _fasta_itr_from_pyfaidx(path: Path) -> Iterator[FastaRecord]:
    """Provide an iteration through FASTA records using pyfaidx."""
    try:
        fasta_file = Fasta(
            str(path), as_raw=True, read_long_names=True, sequence_always_upper=False
        )
    except Exception as exc:  # pragma: no cover - raised via tests for malformed files
        raise MalformedInput(str(exc)) from exc

    try:
        for header in fasta_file.keys():
            yield FastaRecord(header, str(fasta_file[header]))
    except Exception as exc:
        raise MalformedInput(str(exc)) from exc
    finally:
        fasta_file.close()


def _fasta_itr_from_name(fname: str) -> Iterator[FastaRecord]:
    """Provide an iteration through FASTA records in a file path."""
    path = Path(fname)
    if path.suffix.lower() in {".gz", ".bgz"}:
        handle = ezopen(fname)
        try:
            yield from _fasta_itr_from_file(handle)
        finally:
            handle.close()
        return

    yield from _fasta_itr_from_pyfaidx(path)


def _fasta_itr(src: str | PathLike[str] | IO[str] | IO[bytes]) -> Iterator[FastaRecord]:
    """Provide an iterator over FASTA records from paths or open file handles."""
    if isinstance(src, (str, PathLike)):
        return _fasta_itr_from_name(os.fspath(src))
    if hasattr(src, "read"):
        return _fasta_itr_from_file(src)
    raise TypeError(f"Unsupported FASTA source type: {type(src)!r}")


def fasta_get_by_name(
    itr: Iterator[FastaRecord], name: str, byLength: bool = False
) -> FastaRecord | None:
    """Return the record in itr with the given name."""
    target = name.strip()
    target_len = len(target)
    for rec in itr:
        current = rec.header.strip()
        if byLength and target_len < len(current):
            current = current[:target_len]
        if current == target:
            return rec
    return None


class fasta_itr(object):
    """An iterator through a sequence of FASTA records."""

    def __init__(self, src: str | PathLike[str] | IO[str] | IO[bytes]):
        """Create an iterator through the records in src."""
        self.__itr = _fasta_itr(src)

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.__itr)

    # Python 2 compatibility alias
    next = __next__

    def __getitem__(self, name):
        return fasta_get_by_name(iter(self), name)


class fasta_slice(object):
    """Provide an iterator through FASTA records between first and last bounds."""

    def __init__(self, src, first, last=None):
        self.__itr = _fasta_itr(src)
        self.__first = first
        self.__last = last
        if type(first) == int:
            self.__current = 0
        elif type(first) == type(""):
            self.__current = None
        else:
            raise ValueError("bad first")

        self.__foundFirst = False
        if self.__first == 0 or self.__first == "":
            self.__foundFirst = True

    def __iter__(self):
        return self

    def __next__(self):
        if not self.__foundFirst:
            for rec in self.__itr:
                if type(self.__first) == int:
                    if self.__first == self.__current:
                        self.__foundFirst = True
                        break
                    self.__current += 1
                else:
                    if rec.header == self.__first:
                        self.__foundFirst = True
                        break
                    self.__current = rec.header
            if not self.__foundFirst:
                raise ValueError("did not find first record")
            return rec
        rec = next(self.__itr)

        if self.__last is not None:
            if type(self.__first) == int:
                self.__current += 1
                if self.__current == self.__last:
                    raise StopIteration
            else:
                if rec.header == self.__last:
                    raise StopIteration
                self.__current = rec.header
        return rec

    # Python 2 compatibility alias
    next = __next__

    def __getitem__(self, name):
        return fasta_get_by_name(iter(self), name)

    def save(self, fileName):
        with open(fileName, "w", encoding="utf-8") as outfile:
            for record in self:
                outfile.write(str(record))


def get_sequence(src, name):
    """Return the record in src with the given name."""
    return fasta_itr(src)[name]


def fasta_count(src):
    """Count the number of records in a FASTA source."""
    return sum(1 for _rec in fasta_itr(src))


def fasta_split(fileName, num_files, directory=None):
    """Split a FASTA file into a given number of output files."""
    num_records = fasta_count(fileName)

    if directory is None:
        base, ext = os.path.splitext(fileName)
    else:
        _dir, name = os.path.split(fileName)
        base, ext = os.path.splitext(name)
        base = os.path.join(directory, base)

    rec_num = 1
    file_num = 0
    recs_per_file = float(num_records) / float(num_files)
    rec_limit = 0
    outfile = None

    for rec in fasta_itr(fileName):
        if rec_num > round(rec_limit):
            if outfile is not None:
                outfile.close()
            file_num += 1
            outfile = open(base + "." + str(file_num) + ext, "w", encoding="utf-8")
            rec_limit += recs_per_file
        outfile.write(str(rec))
        rec_num += 1

    if outfile is not None:
        outfile.close()


class FastaRandomizer(object):
    """Load all FASTA sequences and return random subsets."""

    def __init__(self, fastaFile):
        self.records = [rec for rec in fasta_itr(fastaFile)]
        self.recRange = range(len(self.records))

    def randomRecords(self, n):
        if n > len(self.records):
            raise ValueError("%d FASTA records requested; only %d stored." % (n, len(self.records)))
        idList = random.sample(self.recRange, n)
        return [self.records[i] for i in idList]

    def __str__(self):
        return "FastaRandomizer - %d records" % len(self.records)


def truncateSequences(fastaFile, exonSize, intronSize, acceptor=False, outFile=None, verbose=False):
    """Truncate FASTA sequences based on given intron and exon sizes."""
    if verbose:
        sys.stderr.write("Loading sequences from %s\n" % fastaFile)
    fiter = fasta_itr(fastaFile)

    outStream = sys.stdout
    if outFile:
        if verbose:
            sys.stderr.write("Writing output to %s\n" % outFile)
        outStream = open(outFile, "w", encoding="utf-8")

    try:
        for rec in fiter:
            midpt = len(rec.sequence) // 2
            if acceptor:
                newRec = FastaRecord(
                    rec.header, rec.sequence[midpt - intronSize : midpt + exonSize]
                )
            else:
                newRec = FastaRecord(
                    rec.header, rec.sequence[midpt - exonSize : midpt + intronSize]
                )
            outStream.write(str(newRec))
    finally:
        if outFile:
            outStream.close()
