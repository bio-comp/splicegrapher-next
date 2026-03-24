# Alignment Collect Lean Cleanup Design

## Goal

Clean up `SpliceGrapher/formats/alignment_io/collect.py` without adding new internal files or changing the public `alignment_io` package API.

## Problem

`collect.py` currently mixes two different responsibilities:
- the real pysam-backed collection engine (`_collect_pysam_data`, junction building, depth accumulation)
- cross-boundary glue that does not belong in a pysam collector (`_is_depths_source`, `_collect_depths_source_data`)

That makes the module impure and makes `api.py` depend on collection internals that are not actually about pysam collection.

## Design

Keep the existing `alignment_io` package shape and make the boundary cleaner:

- `collect.py` becomes the pure pysam-backed collection engine
- `api.py` owns source classification and the depth-file fallback bridge because that logic belongs to the public contract layer
- no new files are introduced
- no public API names change

## Module Responsibilities After the Cut

### `SpliceGrapher/formats/alignment_io/collect.py`
Keep only pysam-driven collection logic:
- `_collect_pysam_data`
- `_build_splice_junctions`
- `_next_match_anchor`
- `_record_strand`
- depth-array helpers tightly coupled to pysam collection

Remove from this module:
- `_is_depths_source`
- `_collect_depths_source_data`
- `depth_io` imports that exist only for fallback bridging

### `SpliceGrapher/formats/alignment_io/api.py`
Own the public contract decisions:
- source-type branching
- depths-file fallback collection
- calls into `collect._collect_pysam_data(...)` only for actual alignment sources

Add private helpers here if needed:
- `_is_depths_source`
- `_collect_depths_source_data`

## Why This Shape

This is the smallest cleanup that improves the package boundary without repeating the over-splitting mistake from the abandoned `#205` attempt.

It removes the impurity from `collect.py` while avoiding another package explosion.

## Non-Goals

This slice does not:
- change `collect_alignment_data(...)` semantics
- change depth/junction calculations
- split `collect.py` into more files
- rename public alignment API functions

## Testing Strategy

Pin the lean boundary rather than internal bloat:
- `collect.py` should no longer expose `_is_depths_source` or `_collect_depths_source_data`
- `api.py` should own those helpers
- focused fallback tests should patch the new seam in `api.py`
- parity tests for BAM/SAM/CRAM behavior must stay green

## Verification

Required gates:
- `uv run ruff check . --fix`
- `uv run ruff format .`
- `uv run mypy .`
- `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
- `uv build`

## Next Step After This Slice

After this lands, the next cleanup candidate should remain:
- `SpliceGrapher/formats/parsers/gene_model_gff_record_handlers.py`
