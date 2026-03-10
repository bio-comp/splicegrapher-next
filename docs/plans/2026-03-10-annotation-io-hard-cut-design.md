# #189 Design: annotation_io hard in-place cleanup

## Summary

`SpliceGrapher/formats/annotation_io.py` is still a mixed-responsibility formats module with a weak public boundary: `load_gene_models(path: str, **args: object)`. This slice keeps the module in place but hard-cuts that API to an explicit typed signature, removes `**args`, and tightens the internal helper layout while preserving behavior.

## Scope

In scope:
- keep `SpliceGrapher/formats/annotation_io.py` as a single module
- replace `load_gene_models(path: str, **args: object)` with an explicit typed signature
- preserve supported named parameters:
  - `outdir`
  - `cache_dir`
- make unknown keyword arguments fail immediately via normal Python signature enforcement
- tighten internal helper typing and section boundaries without changing the module’s location

Out of scope:
- splitting `annotation_io.py` into multiple files
- changing the core gffutils-backed loading strategy
- broader `GeneModel` or parser architecture changes

## Why this shape

The file is only 437 lines. It is not a monolith that requires another package split right now. The real pollution is the boundary and loose option handling, not the file path.

A hard in-place cleanup gives the right outcome with less churn:
- explicit public contract
- fewer casts and `object`-typed option flows
- no fake compatibility layer
- no new module topology for reviewers to absorb

## Intended public API

`load_gene_models(`
- `path: str | Path,`
- `*,`
- `outdir: str | Path | None = None,`
- `cache_dir: str | Path | None = None,`
`) -> GeneModel`

That is the only public boundary change intended in this slice.

## Expected fallout

Direct tests already exist in `tests/test_annotation_io.py`, and `tests/test_integration_simple.py` also covers the loader.

The likely failures from the hard cut are:
- tests that assume permissive keyword handling
- typing fallout inside `annotation_io.py` where `cast(...)` is currently compensating for `**args: object`

This should stay a reviewable formats-layer cleanup PR rather than another structural refactor.
