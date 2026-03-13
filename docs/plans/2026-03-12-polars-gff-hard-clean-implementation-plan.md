# Polars GFF Hard-Clean Implementation Plan

1. Add a failing boundary test that asserts `polars_gff.py` no longer relies on `object` in runtime helper signatures.
2. Introduce typed protocols for flattened feature extraction.
3. Replace `object`-typed helpers in `polars_gff.py` with protocol-based helpers.
4. Run focused checks for `polars_gff` and its tests.
5. Run full SGN gates:
   - `uv run ruff check . --fix`
   - `uv run ruff format .`
   - `uv run mypy .`
   - `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
   - `uv build`
