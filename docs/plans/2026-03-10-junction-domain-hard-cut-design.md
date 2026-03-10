# #187 Design: junction/domain hard cut

## Summary

The remaining live legacy API residue at the SGN junction/depth boundary is concentrated in `SpliceGrapher/formats/junction.py` and propagated into `depth_io.py`, `alignment_io.py`, and the junction/depth tests. This slice hard-cuts that boundary to snake_case and cleans the adjacent consumer code in the same PR.

## Scope

In scope:
- hard-cut `SpliceJunction.sjCode` to `sj_code`
- hard-cut `SpliceJunction.minAnchor()` to `min_anchor()`
- hard-cut `SpliceJunction.toString()` to `to_string()`
- update direct SGN callers in `depth_io.py`, `alignment_io.py`, and tests
- clean the related protocol/fixture expectations while fixing fallout

Out of scope:
- broader `depth_io.py` architectural changes
- unrelated parser/package reorganization
- compatibility aliases or wrapper properties

## Boundary

The design keeps `SpliceGrapher/formats/junction.py` and `parse_junction_record()` in place. The hard cut is limited to the junction model surface and the narrow band of adjacent depth/alignment consumers that currently depend on the legacy names.

## Caller blast radius

Confirmed runtime/test fallout is limited to:
- `SpliceGrapher/formats/junction.py`
- `SpliceGrapher/formats/depth_io.py`
- `SpliceGrapher/formats/alignment_io.py`
- `tests/test_junction.py`
- `tests/test_depth_io.py`

That keeps the slice reviewable while still removing the full live residue pocket instead of preserving compatibility.

## Hard-cut policy

- no `sjCode` field after this PR
- no `minAnchor()` method after this PR
- no `toString()` method after this PR
- no compatibility aliases
- callers and tests move in the same PR
