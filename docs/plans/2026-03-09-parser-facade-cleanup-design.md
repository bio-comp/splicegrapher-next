# Parser Facade Cleanup Design

## Context

`SpliceGrapher.formats.parsers` still has two pieces of parser debt after the gene-model and splice-graph refactors:

1. `SpliceGrapher/formats/parsers/__init__.py` exposes only `gene_model_gff`, which is an incomplete and leaky package facade.
2. `SpliceGrapher/formats/parsers/splice_graph.py` still uses manual iterator bookkeeping and a camelCase loader method:
   - `self._graph_keys`
   - `self.graphId`
   - `__next__`
   - `loadFromFile`

The user explicitly wants a hard cut on `loadFromFile`: remove the old name and fix the fallout in one PR.

## Goal

Tighten the parser package boundary and modernize `SpliceGraphParser` without expanding into a broader parser API rewrite.

Concretely:
- `SpliceGrapher.formats.parsers` should export the parser entrypoints callers actually need
- `SpliceGraphParser` should iterate using standard Python iteration over parsed graph values
- `loadFromFile` should be renamed to `load_from_file`, with any fallout fixed immediately

## Scope

### In scope
- `SpliceGrapher/formats/parsers/__init__.py`
- `SpliceGrapher/formats/parsers/splice_graph.py`
- tests that currently import parser modules through `SpliceGrapher.formats.parsers`
- tests or runtime callers that still reference `loadFromFile`

### Out of scope
- removing eager parsing from `SpliceGraphParser.__init__`
- moving the gene-model parser into a `gene_model/` sub-package
- rewriting `ParseContext` to be functional or immutable
- changing `SpliceGraphParser.graphDict` storage semantics

## Design

### Parser package facade

`SpliceGrapher/formats/parsers/__init__.py` should become an explicit public boundary:

- `from .gene_model_gff import load_gene_model_records`
- `from .splice_graph import SpliceGraphParser`

This reduces namespace leakage and gives callers a stable import surface.

### SpliceGraphParser cleanup

`SpliceGraphParser` should:
- keep eager parsing in `__init__` for now
- rename `loadFromFile` to `load_from_file`
- remove `self._graph_keys` and `self.graphId`
- remove `__next__`
- implement `__iter__` as `return iter(self.graphDict.values())`

This keeps external behavior intact for `for graph in parser` and `list(parser)` while deleting the remaining Python 2-style iterator scaffolding.

## Testing strategy

We need a deliberate break-and-fix cycle.

### Red phase
Write or update tests so they target the desired boundary:
- `from SpliceGrapher.formats.parsers import load_gene_model_records, SpliceGraphParser`
- `for graph in SpliceGraphParser(...)` still works
- `parser.load_from_file()` exists and `loadFromFile` no longer does, if we choose to pin that explicitly

### Green phase
Implement the facade and iterator cleanup, then run:
- targeted parser tests
- targeted gene-model boundary tests
- full SGN `pytest`
- `uv build`

## Risks

### Risk: caller code still reaches into `parsers` for module objects
Current tests do this with `gene_model_gff` and related support modules. Tightening `__init__.py` may require updating those tests to import the real submodules directly instead of relying on `parsers.__init__` leakage.

### Risk: some code still calls `loadFromFile`
That is why this PR must intentionally break the old name and sweep the fallout in one go.

### Risk: changing iteration breaks parser consumers
Low risk, but we still need explicit coverage for `iter(parser)` behavior.

## Recommendation

Implement this as one narrow PR:
1. tighten the parser facade
2. remove `SpliceGraphParser` manual iterator state
3. rename `loadFromFile` to `load_from_file`
4. fix all fallout immediately
