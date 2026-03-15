# Core SpliceGraph Package Cut Design

## Goal

Replace the flat `SpliceGrapher/core/splice_graph.py` monolith with a concrete package layout that separates constants, node behavior, and graph behavior, while updating SGN imports directly to those concrete modules.

## Why This Slice

`SpliceGrapher/core/splice_graph.py` is still a dense runtime monolith. It currently combines:
- graph-level constants and coercion helpers
- the `SpliceGraphNode` data model
- the `SpliceGraph` topology manager
- invariant-sensitive file-path assumptions in tests and CI guards

That shape no longer matches the package-first structure now used across `formats/` and `shared/`. This slice removes that structural pollution without reopening the parser/writer contracts or reintroducing compatibility shims.

## Non-Goals

This slice does not:
- change networkx graph semantics
- change splice-graph parser behavior
- change splice-graph writer output
- add a compatibility shim for the old flat-file path
- preserve broad package-root re-exports for hypothetical downstream callers

## Target Layout

Delete:
- `SpliceGrapher/core/splice_graph.py`

Add:
- `SpliceGrapher/core/splice_graph/__init__.py`
- `SpliceGrapher/core/splice_graph/constants.py`
- `SpliceGrapher/core/splice_graph/node.py`
- `SpliceGrapher/core/splice_graph/graph.py`

## Module Responsibilities

### `constants.py`
Owns shared splice-graph constants and shared type aliases:
- record and attribute constants consumed by parser/writer/tests
- `NodeAttributeValue`
- `_coerce_alt_splicing_event(...)` if it remains shared by the node model

### `node.py`
Owns interval-only and node-level data behavior:
- `NullNode`
- `SpliceGraphNode`
- node attribute, codon, isoform, and alt-splicing helpers

### `graph.py`
Owns graph-level orchestration and topology:
- `SpliceGraph`
- graph-owned node lookup, edge creation, traversal, and validation

### `__init__.py`
Minimal package plumbing only.
It may re-export the handful of symbols SGN intentionally imports from the package root during the transition, but it is not a compatibility layer for removed path names.

## Import Rewrite Policy

Be harsh and explicit:
- parser and writer modules should import constants from `constants.py`
- graph consumers should import `SpliceGraph` from `graph.py`
- node-specific consumers should import `SpliceGraphNode` from `node.py`
- tests should prefer concrete module imports where they are asserting layout or ownership

This avoids recreating another flat aggregator disguised as a package.

## Blast Radius

Known SGN touch points:
- `SpliceGrapher/core/graph_math.py`
- `SpliceGrapher/core/splicing_events.py`
- `SpliceGrapher/formats/parsers/splice_graph.py`
- `SpliceGrapher/formats/writers/splice_graph.py`
- splice-graph tests and integration tests
- invariant guard files that currently hard-code `SpliceGrapher/core/splice_graph.py`
- `PROVENANCE.md` current-source inventory language

## Invariant and Test Impact

The clean-invariant ratchet currently protects the flat file path. This slice must update:
- `scripts/ci/check_clean_invariant.py`
- `tests/test_clean_invariant_guard.py`

New or updated layout checks should verify:
- the flat file is gone
- the package directory exists
- the concrete modules contain the intended owners

## Verification

Required gates after the cut:
- targeted splice-graph/core test slice
- `uv run ruff check . --fix`
- `uv run ruff format .`
- `uv run mypy .`
- `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
- `uv build`

## Recommended Next Step After This Slice

If this cut lands cleanly, the next structural cleanup candidate should shift to `SpliceGrapher/formats/gene_model/model.py`, which remains one of the larger runtime containers after the recent package conversions.
