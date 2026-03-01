# Parity Harness (T1)

Issue: `#23` Testing tranche T1 parity harness expansion.

## Purpose

Keep cleanup and extraction slices gated by explicit parity coverage for core/formats behavior.

## Required Targeted Test Matrix

The canonical matrix is tracked in `tests/fixtures/parity_harness.toml`.

Required areas:

- `core_shared_import_smoke` -> `tests/test_core_shared_import_smoke.py`
- `annotation_io` -> `tests/test_annotation_io.py`
- `alignment_io_parity` -> `tests/test_alignment_io_parity.py`
- `splicegrapher_alignment_io` -> `tests/test_splicegrapher_alignment_io.py`

## Fixture Conventions

Fixture conventions and builder pointers are documented in:

- `tests/fixtures/README.md`

## PR Expectations

For cleanup/extraction PRs that touch core or formats paths:

1. Run full gates:
   - `uv run ruff check .`
   - `uv run ruff format --check .`
   - `uv run pytest -q`
2. Ensure any changed parity area still maps in `tests/fixtures/parity_harness.toml`.
3. If new parity areas are added, update:
   - `tests/fixtures/parity_harness.toml`
   - `tests/test_parity_harness_matrix.py`
   - this document.
