"""Low-level FASTA readers and typed iterator utilities."""

from __future__ import annotations

import os
from os import PathLike
from pathlib import Path
from typing import IO, Iterator

from pyfaidx import Fasta

from SpliceGrapher.formats.fasta.records import FastaRecord, MalformedInput
from SpliceGrapher.shared.file_utils import ez_open

FastaSource = str | PathLike[str] | IO[str] | IO[bytes]
SliceBound = int | str


def _to_text(line: str | bytes) -> str:
    if isinstance(line, bytes):
        return line.decode("utf-8")
    return line


def _iter_fasta_records_from_file(file_handle: IO[str] | IO[bytes]) -> Iterator[FastaRecord]:
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


def _iter_fasta_records_from_pyfaidx(path: Path) -> Iterator[FastaRecord]:
    try:
        fasta_file = Fasta(
            str(path),
            as_raw=True,
            read_long_names=True,
            sequence_always_upper=False,
        )
    except Exception as exc:  # pragma: no cover
        raise MalformedInput(str(exc)) from exc

    try:
        for header in fasta_file.keys():
            yield FastaRecord(header, str(fasta_file[header]))
    except Exception as exc:
        raise MalformedInput(str(exc)) from exc
    finally:
        fasta_file.close()


def _iter_fasta_records_from_name(source_name: str) -> Iterator[FastaRecord]:
    path = Path(source_name)
    if path.suffix.lower() in {".gz", ".bgz"}:
        handle = ez_open(source_name)
        try:
            yield from _iter_fasta_records_from_file(handle)
        finally:
            handle.close()
        return

    yield from _iter_fasta_records_from_pyfaidx(path)


def iter_fasta_records(source: FastaSource) -> Iterator[FastaRecord]:
    if isinstance(source, (str, PathLike)):
        return _iter_fasta_records_from_name(os.fspath(source))
    if hasattr(source, "read"):
        return _iter_fasta_records_from_file(source)
    raise TypeError(f"Unsupported FASTA source type: {type(source)!r}")


def fasta_get_by_name(
    records: Iterator[FastaRecord],
    name: str,
    *,
    by_length: bool = False,
) -> FastaRecord | None:
    target = name.strip()
    target_length = len(target)
    for record in records:
        current = record.header.strip()
        if by_length and target_length < len(current):
            current = current[:target_length]
        if current == target:
            return record
    return None


class FastaIterator:
    """Iterate through FASTA records from a supported source."""

    def __init__(self, source: FastaSource):
        self._records = iter_fasta_records(source)

    def __iter__(self) -> FastaIterator:
        return self

    def __next__(self) -> FastaRecord:
        return next(self._records)

    def __getitem__(self, name: str) -> FastaRecord | None:
        return fasta_get_by_name(iter(self), name)


class FastaSlice:
    """Iterate through FASTA records between the given first/last bounds."""

    def __init__(self, source: FastaSource, first: SliceBound, last: SliceBound | None = None):
        self._records = iter_fasta_records(source)
        self._first = first
        self._last = last
        self._current: int | str | None = 0 if isinstance(first, int) else None
        self._found_first = first in {0, ""}

    def __iter__(self) -> FastaSlice:
        return self

    def __next__(self) -> FastaRecord:
        if not self._found_first:
            record = self._find_first_record()
            self._found_first = True
            return record

        record = next(self._records)
        if self._last is not None:
            if isinstance(self._first, int):
                assert isinstance(self._current, int)
                self._current += 1
                if self._current == self._last:
                    raise StopIteration
            else:
                if record.header == self._last:
                    raise StopIteration
                self._current = record.header
        return record

    def __getitem__(self, name: str) -> FastaRecord | None:
        return fasta_get_by_name(iter(self), name)

    def save(self, file_name: str | PathLike[str]) -> None:
        with open(file_name, "w", encoding="utf-8") as output_handle:
            for record in self:
                output_handle.write(str(record))

    def _find_first_record(self) -> FastaRecord:
        for record in self._records:
            if isinstance(self._first, int):
                assert isinstance(self._current, int)
                if self._first == self._current:
                    return record
                self._current += 1
            else:
                if record.header == self._first:
                    return record
                self._current = record.header
        raise ValueError("did not find first record")
