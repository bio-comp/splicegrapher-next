# GeneModel Internal Split Design

## Goal

Reduce `SpliceGrapher/formats/gene_model/model.py` from a single container/query monolith into smaller internal modules, while keeping `SpliceGrapher.formats.gene_model.GeneModel` as the stable public type.

## Why This Slice

`SpliceGrapher/formats/gene_model/model.py` is still one of the largest remaining runtime modules in SGN. It currently mixes three distinct concerns:
- container state and mutation (`add_gene`, `add_chromosome`)
- query and lookup behavior (`get_gene`, `get_genes_in_range`, `isoform_dict`, donor/acceptor aggregation)
- loading/index/serialization delegation (`from_gff`, `load_gene_model`, `make_sorted_model`, `write_gff`, `write_gtf`)

The package boundary for `SpliceGrapher.formats.gene_model` already exists, so this slice can remove structural pollution without another namespace migration.

## Non-Goals

This slice does not:
- change GeneModel semantics
- change parser or writer behavior
- rename public GeneModel methods
- move domain entities out of `SpliceGrapher.formats.models`
- add compatibility aliases beyond the current `SpliceGrapher.formats.gene_model` package surface

## Target Shape

Keep:
- `SpliceGrapher/formats/gene_model/model.py`
- `SpliceGrapher/formats/gene_model/__init__.py`
- `SpliceGrapher/formats/gene_model/repository.py`

Add:
- `SpliceGrapher/formats/gene_model/queries.py`
- `SpliceGrapher/formats/gene_model/loading.py`

## Responsibility Split

### `model.py`
Owns the `GeneModel` dataclass and core state only:
- field definitions
- minimal state mutation helpers (`add_gene`, `add_chromosome`, parent lookup helper)
- thin delegators to query/loading helpers

### `queries.py`
Owns read/query logic over an existing `GeneModel` state:
- all-gene and per-chromosome lookups
- range and parent lookups
- donor/acceptor aggregation
- isoform dictionary generation
- annotation parsing helpers if they remain query-adjacent

### `loading.py`
Owns model population and index maintenance helpers:
- `from_gff` support path validation
- `load_gene_model` delegation wrapper
- `make_sorted_model`
- any small helpers tied directly to initialization/loading flow

### `repository.py`
Remains the boundary for loading and writing. This slice should not move repository behavior back into `model.py`.

## Import Policy

Be conservative at the package boundary and explicit internally:
- `SpliceGrapher.formats.gene_model.GeneModel` stays where it is
- `model.py` may import internal helper functions from `queries.py` and `loading.py`
- tests that assert internal layout should import the concrete internal modules directly
- external SGN code should continue to prefer `SpliceGrapher.formats.gene_model`

## Expected Blast Radius

Primary touch points:
- `SpliceGrapher/formats/gene_model/model.py`
- `SpliceGrapher/formats/gene_model/__init__.py`
- new internal helper modules
- `tests/test_gene_model.py`
- `tests/test_gene_model_package_layout.py`
- any tests that inspect method signatures or monkeypatch `GeneModel` methods

## Testing Strategy

This should be TDD and structural:
- first pin the internal package layout and current public `GeneModel` boundary in tests
- then extract one concern at a time behind the existing public class
- run targeted `gene_model` tests first
- then run full SGN quality gates

## Verification

Required gates after the split:
- targeted `tests/test_gene_model.py`
- targeted `tests/test_gene_model_package_layout.py`
- `uv run ruff check . --fix`
- `uv run ruff format .`
- `uv run mypy .`
- `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
- `uv build`

## Next Step After This Slice

If this lands cleanly, the next cleanup candidate should shift to either:
- `SpliceGrapher/formats/alignment_io/collect.py`, or
- `SpliceGrapher/formats/parsers/gene_model_gff_record_handlers.py`

Both are still large internal modules, but `GeneModel` is the cleaner structural win right now.
