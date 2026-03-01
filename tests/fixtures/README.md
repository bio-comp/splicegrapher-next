# Test Fixture Conventions

This directory contains fixture metadata and tiny deterministic test assets used by parity-focused test modules.

## Rules

- Keep fixture files small and deterministic.
- Prefer synthetic or reduced fixture inputs over large upstream datasets.
- Track parity matrix requirements in `tests/fixtures/parity_harness.toml`.
- Keep all paths in the manifest repository-relative.

## Required Parity Areas

- `core_shared_import_smoke`
- `annotation_io`
- `alignment_io_parity`
- `splicegrapher_alignment_io`

## Fixture Builders

- `tests/helpers/alignment_fixture_builder.py`
- `tests/helpers/idiffir_fixture_builder.py`
