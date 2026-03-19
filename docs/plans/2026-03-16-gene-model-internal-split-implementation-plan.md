# GeneModel Lean Cleanup Implementation Plan

Goal: remove the `GeneModelRepository` shim and the dead `GeneModel.load_gene_model(...)` wrapper without splitting `model.py` further.

## Task 1: Pin the lean contract in tests

Files:
- `tests/test_gene_model_package_layout.py`
- `tests/test_gene_model.py`

Changes:
- assert `GeneModelRepository` is no longer exported
- assert `GeneModel.load_gene_model(...)` no longer exists
- rewrite parser-loading tests to use `load_gene_model_records(...)` directly
- keep writer seam tests, but patch the actual direct call sites in `model.py`

Verification:
- `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_gene_model_package_layout.py tests/test_gene_model.py`

## Task 2: Delete the repository shell

Files:
- delete `SpliceGrapher/formats/gene_model/repository.py`
- modify `SpliceGrapher/formats/gene_model/__init__.py`
- modify `SpliceGrapher/formats/gene_model/model.py`

Changes:
- remove `GeneModelRepository` import/export
- inline direct parser call in `from_gff`
- inline direct writer calls in `write_gff` and `write_gtf`
- delete `GeneModel.load_gene_model(...)`

Verification:
- `uv run ruff check SpliceGrapher/formats/gene_model tests/test_gene_model.py tests/test_gene_model_package_layout.py`
- `uv run mypy SpliceGrapher/formats/gene_model tests/test_gene_model.py tests/test_gene_model_package_layout.py`
- `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_gene_model_package_layout.py tests/test_gene_model.py`

## Task 3: Run full SGN verification

Run:
- `uv run ruff check . --fix`
- `uv run ruff format .`
- `uv run mypy .`
- `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
- `uv build`

Expected outcome:
- full repo stays green
- diff is net-negative and removes a dead layer rather than adding new internal modules
