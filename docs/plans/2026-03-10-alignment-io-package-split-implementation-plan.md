# Alignment I/O Package Split Implementation Plan

1. Add/adjust layout tests to pin the new `alignment_io` package structure.
2. Extract shared type aliases and protocols into `types.py`.
3. Extract source/path opening and alignment streaming helpers into `sources.py`.
4. Extract pysam-backed collection logic into `collect.py`.
5. Extract gene-depth calculation helpers into `depths.py`.
6. Replace the flat module with `alignment_io/__init__.py` and `api.py` re-exports.
7. Run focused alignment/depth verification, then full repo gates.
8. Commit the structural split as one reviewable refactor tranche.
