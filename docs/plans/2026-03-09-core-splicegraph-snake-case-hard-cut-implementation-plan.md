# #185 Implementation plan: Core SpliceGraph snake_case hard cut

1. Rewrite the splice-graph boundary tests and integration tests to the snake_case surface.
2. Rename the node-level API in `SpliceGrapher/core/splice_graph.py` and update direct callers.
3. Rename the graph-level API in `SpliceGrapher/core/splice_graph.py` and update direct callers.
4. Update parser/writer/runtime modules:
   - `SpliceGrapher/formats/parsers/splice_graph.py`
   - `SpliceGrapher/formats/writers/splice_graph.py`
   - `SpliceGrapher/core/splicing_events.py`
   - `SpliceGrapher/core/graph_math.py`
5. Run focused verification:
   - splice-graph tests
   - graph-math tests
   - splicing-events tests
   - integration roundtrip tests
6. Run full repo verification:
   - `uv run ruff check . --fix`
   - `uv run ruff format .`
   - `uv run mypy .`
   - `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
   - `uv build`
7. Commit the docs and implementation separately if the slice stays reviewable.
