# GeneModel Lean Cleanup Design

## Goal

Remove dead indirection from `SpliceGrapher/formats/gene_model` without adding new internal scaffolding.

## Re-Scope

The original plan for `#205` was to split `SpliceGrapher/formats/gene_model/model.py` into additional internal modules. That turned into unnecessary bloat:
- new tiny files
- mixin machinery
- tests pinning internal layout instead of runtime behavior

`model.py` on `main` is only about 420 lines. That is not the real problem.

The real smells are:
- `SpliceGrapher/formats/gene_model/repository.py` as a pure pass-through shell
- `GeneModel.load_gene_model(...)` as a dead wrapper around the parser entrypoint

## Actual Target Shape

Keep:
- `SpliceGrapher/formats/gene_model/model.py`
- `SpliceGrapher/formats/gene_model/__init__.py`

Delete:
- `SpliceGrapher/formats/gene_model/repository.py`
- `GeneModel.load_gene_model(...)`
- `GeneModelRepository` package export

## Design Decisions

### Keep `model.py` whole
`GeneModel` still owns:
- container state
- mutation helpers
- query methods
- `from_gff`
- `make_sorted_model`
- writer delegation

That keeps the runtime shape simple and avoids over-modularizing a file that is not yet too large.

### Inline parser and writer calls
Replace the repository shell with direct runtime calls:
- `from_gff` calls `load_gene_model_records(...)` directly
- `write_gff` calls `write_gene_model_gff(...)` directly
- `write_gtf` calls `write_gene_model_gtf(...)` directly

### Remove dead wrapper surface
`GeneModel.load_gene_model(...)` is not used by SGN runtime. Only tests were calling it.

That means it should be deleted instead of preserved.

## Non-Goals

This slice does not:
- split `model.py` further
- change parser behavior
- change writer behavior
- rename the surviving public `GeneModel` query methods

## Testing Strategy

Pin the lean contract:
- `GeneModelRepository` is gone from the package
- `GeneModel.load_gene_model(...)` is gone from the class
- parser-loading tests call `load_gene_model_records(...)` directly
- writer tests patch the direct writer call sites in `model.py`

## Verification

Required gates:
- `uv run ruff check . --fix`
- `uv run ruff format .`
- `uv run mypy .`
- `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
- `uv build`

## Next Step After This Slice

After this lands, the next cleanup candidates remain:
- `SpliceGrapher/formats/alignment_io/collect.py`
- `SpliceGrapher/formats/parsers/gene_model_gff_record_handlers.py`
