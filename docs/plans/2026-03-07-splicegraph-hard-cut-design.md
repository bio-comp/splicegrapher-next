# SpliceGraph Hard-Cut Design

Date: 2026-03-07
Issue: #173
Branch: `sg-next/173-split-splicegraph-core-parser`

## Goal

Delete `SpliceGrapher/SpliceGraph.py` and replace it with two lowercase modules:

- `SpliceGrapher/core/splice_graph.py`
- `SpliceGrapher/formats/parsers/splice_graph.py`

This is a hard cut. There will be no PascalCase compatibility shim.

## Current Problems

1. `SpliceGrapher/SpliceGraph.py` still violates the repo layout target.
   - PascalCase module name.
   - graph core and parser live in the same file.
2. Internal SGN code still imports the deleted legacy module path.
3. CI guards and smoke tests still encode the legacy path.

## Decisions

### 1. No shim

`SpliceGrapher/SpliceGraph.py` will be removed outright.

Reason:
- SGN has already accepted hard-cut cleanup for `alignment_io.py`.
- Downstream Python 2 compatibility is not a gating concern for this tranche.
- A shim would preserve naming debt and add temporary surface area we already know we want to delete.

### 2. Split core and parser by responsibility

`SpliceGrapher/core/splice_graph.py` will contain:
- graph constants required by runtime code and tests
- `NullNode`
- `SpliceGraphNode`
- `SpliceGraph`

`SpliceGrapher/formats/parsers/splice_graph.py` will contain:
- `SpliceGraphParser`

### 3. Keep the `SpliceGrapher` namespace, not the legacy module path

This is a repo-layout cleanup, not a package-namespace change.
The import namespace stays under `SpliceGrapher`, but callers must use the new lowercase modules.

### 4. Update invariant guards with the new path

Current clean-invariant logic explicitly names `SpliceGrapher/SpliceGraph.py` in enum and overlap protected-path lists. Those protections should move to `SpliceGrapher/core/splice_graph.py` as part of the cut.

## Known Import Blast Radius

Runtime:
- `SpliceGrapher/core/graph_math.py`
- `SpliceGrapher/core/splicing_events.py`
- `SpliceGrapher/formats/writers/splice_graph.py`

Tests:
- `tests/test_graph_math.py`
- `tests/test_splice_graph.py`
- `tests/test_splicing_events.py`
- `tests/test_integration_simple.py`
- `tests/test_core_shared_import_smoke.py`
- `tests/test_clean_invariant_guard.py`

Tooling/docs:
- `scripts/ci/check_clean_invariant.py`
- docs plans and ADR references to the old filename
- `PROVENANCE.md`

## TDD Strategy

1. Rewrite import/smoke tests to the new module paths first.
2. Run the targeted test slice and capture the expected failures.
3. Move code into the new lowercase modules.
4. Delete `SpliceGrapher/SpliceGraph.py`.
5. Update runtime imports, guards, and references.
6. Run the full SGN suite.

## Non-Goals

- No compatibility re-export from `SpliceGrapher/__init__.py`
- No additional graph-behavior changes beyond the file/module split
- No downstream iDiffIR/TAPIS migration work in this branch
