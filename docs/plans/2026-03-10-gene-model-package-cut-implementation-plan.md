# GeneModel Package Cut Implementation Plan

1. Create a failing layout/boundary test for the `formats.gene_model` package shape.
2. Convert `SpliceGrapher/formats/gene_model.py` into a package directory.
3. Extract `GeneModelRepository` into `repository.py`.
4. Extract `GeneModel` into `model.py` and keep internal imports coherent.
5. Rebuild `__init__.py` as the import-stable facade and re-export boundary.
6. Run focused checks for gene-model modules and tests.
7. Run full SGN gates:
   - `uv run ruff check . --fix`
   - `uv run ruff format .`
   - `uv run mypy .`
   - `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
   - `uv build`
