# SpliceGraph Hard-Cut Implementation Plan

Date: 2026-03-07
Issue: #173
Branch: `sg-next/173-split-splicegraph-core-parser`

## Step 1. Rewrite tests/imports to the target module layout

Update tests to import from:
- `SpliceGrapher.core.splice_graph`
- `SpliceGrapher.formats.parsers.splice_graph`

At minimum touch:
- `tests/test_splice_graph.py`
- `tests/test_graph_math.py`
- `tests/test_splicing_events.py`
- `tests/test_integration_simple.py`
- `tests/test_core_shared_import_smoke.py`

Run the targeted pytest slice and confirm failures are import-path related.

## Step 2. Create the new core module

Create `SpliceGrapher/core/splice_graph.py` and move:
- graph constants used by runtime/writer/tests
- `NullNode`
- `SpliceGraphNode`
- `SpliceGraph`

Update imports in:
- `SpliceGrapher/core/graph_math.py`
- `SpliceGrapher/core/splicing_events.py`
- `SpliceGrapher/formats/writers/splice_graph.py`

## Step 3. Create the parser module

Create `SpliceGrapher/formats/parsers/splice_graph.py` and move `SpliceGraphParser`.

Keep parser behavior unchanged aside from import/module path adjustments.

## Step 4. Delete the PascalCase module

Delete `SpliceGrapher/SpliceGraph.py` once all runtime and tests import the new modules.

## Step 5. Update guards and references

Update:
- `scripts/ci/check_clean_invariant.py`
- `tests/test_clean_invariant_guard.py`
- `PROVENANCE.md`
- docs references that would become stale inside the repo

## Step 6. Verification

Run:
- `uv run ruff check ...`
- `uv run ruff format --check ...`
- `uv run mypy ...`
- targeted splice-graph pytest slice
- full `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`

## Commit Plan

1. `docs: add SpliceGraph hard-cut design and plan`
2. `test: rewrite SpliceGraph imports to lowercase module targets`
3. `refactor: split SpliceGraph core and parser modules`
