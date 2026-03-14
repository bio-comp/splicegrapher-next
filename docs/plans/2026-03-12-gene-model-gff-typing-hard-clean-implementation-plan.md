# Gene Model GFF Helper Typing Hard-Clean Implementation Plan

1. Add a failing boundary test that asserts the parser helper signatures no longer use `object` for runtime inputs.
2. Tighten helper signatures in `gene_model_gff_records.py` and `gene_model_gff_resolution.py`.
3. Run focused parser checks and tests.
4. Run full SGN gates:
   - `uv run ruff check . --fix`
   - `uv run ruff format .`
   - `uv run mypy .`
   - `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
   - `uv build`
