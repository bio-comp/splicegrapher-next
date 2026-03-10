# #187 Implementation plan: junction/domain hard cut

1. Rewrite the junction and depth tests to the snake_case surface first.
2. Rename the `SpliceJunction` field and methods in `SpliceGrapher/formats/junction.py`.
3. Update the `depth_io.py` protocol and junction filtering calls.
4. Update the `alignment_io.py` junction filtering calls.
5. Run focused verification:
   - `tests/test_junction.py`
   - `tests/test_depth_io.py`
   - alignment/depth boundary tests if fallout reaches them
6. Run full repo verification:
   - `uv run ruff check . --fix`
   - `uv run ruff format .`
   - `uv run mypy .`
   - `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
   - `uv build`
7. Keep local tracker files unstaged and open the PR with migration-impact notes.
