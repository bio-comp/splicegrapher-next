# FASTA Package Hard Cut Design

## Summary

Hard-clean `SpliceGrapher/formats/fasta.py` by replacing the flat legacy module with a modern `SpliceGrapher.formats.fasta` package. The cut removes explicit Python 2 residue, deletes legacy camelCase/public alias debt, and separates record, reader, and operation concerns into a clearer package boundary.

## Goals

- Remove explicit Python 2 patterns from the FASTA module.
- Replace the flat module with a package directory that fits the current `formats/` layout.
- Expose a snake_case-first public API at `SpliceGrapher.formats.fasta`.
- Split cohesive concerns so FASTA parsing logic is easier to maintain and test.
- Preserve the top-level import namespace `SpliceGrapher.formats.fasta`.

## Non-Goals

- No broader reorganization of `formats/` beyond the FASTA package boundary.
- No compatibility wrappers for deleted legacy names like `fasta_itr`, `fasta_slice`, or `truncateSequences`.
- No changes to unrelated sequence or alignment modules in this slice.

## Current Problems

The existing `SpliceGrapher/formats/fasta.py` still carries direct legacy debt:

- `class fasta_itr(object)` and `class fasta_slice(object)`
- `next = __next__` aliases
- mixed camelCase and snake_case APIs
- `%`-style string interpolation
- direct `sys.stderr.write(...)` calls
- multiple responsibilities in one file: data records, iterators, file adapters, randomization, and transformation helpers

This is one of the clearest remaining “repo pollution” targets because the smell is explicit and the internal SGN call surface is currently limited to tests.

## Proposed Package Layout

Replace the flat file with a package:

```text
SpliceGrapher/formats/fasta/
├── __init__.py
├── records.py
├── readers.py
└── operations.py
```

Responsibilities:

- `records.py`
  - `MalformedInput`
  - `FastaRecord`
- `readers.py`
  - low-level FASTA iteration/adapters
  - `FastaIterator`
  - `FastaSlice`
- `operations.py`
  - lookup/count/split/truncate/randomize operations
  - `fasta_get_by_name`
  - `get_sequence`
  - `fasta_count`
  - `fasta_split`
  - `truncate_sequences`
  - `FastaRandomizer`
- `__init__.py`
  - public import boundary and explicit `__all__`

## Public API Contract

The public namespace remains `SpliceGrapher.formats.fasta`, but the exported names become modern-only:

- `MalformedInput`
- `FastaRecord`
- `FastaIterator`
- `FastaSlice`
- `FastaRandomizer`
- `fasta_get_by_name`
- `get_sequence`
- `fasta_count`
- `fasta_split`
- `truncate_sequences`

Deleted names:

- `fasta_itr`
- `fasta_slice`
- `truncateSequences`

This is an intentional hard cut. Callers must move to the new names.

## Design Details

### Records

`FastaRecord` stays as the immutable value object for FASTA headers and sequences. `MalformedInput` stays as the format-specific exception.

### Readers

Reader code stays responsible for turning file paths, streams, and pyfaidx-backed sources into `FastaRecord` iterators. The iterator classes become properly named Python 3 classes with typed constructor inputs and no Python 2 aliases.

### Operations

Higher-level utilities move out of the reader classes and into explicit functions/classes in `operations.py`. `truncate_sequences` replaces the old `truncateSequences` helper and logs through `structlog` instead of writing directly to stderr.

## Error Handling

- Continue raising `MalformedInput` for invalid FASTA input.
- Fail fast for bad source types rather than keeping permissive legacy fallback behavior.
- Preserve semantic behavior where it is already covered by tests.

## Testing Strategy

Update `tests/test_fasta_io.py` to the new package surface and add a module-layout guard that asserts:

- `SpliceGrapher.formats.fasta` resolves as a package
- the expected public names are exported
- deleted legacy aliases are absent

Verification for the slice:

- targeted FASTA tests
- Ruff on touched files
- MyPy on touched files
- full pytest
- `uv build`

## Migration Impact

Internal SGN impact should be minimal because the current runtime search only found test call sites. Downstream repos using deleted legacy FASTA names will break and need explicit updates; this is acceptable for this slice because the goal is to remove legacy pollution rather than preserve ghost APIs.
