# Polars Benchmark Relocation Design

## Goal

Move `SpliceGrapher/formats/polars_gff_benchmark.py` into a repo-local top-level
`benchmarks/` package so benchmark helpers stay importable for local
development, test runs, and CI without shipping inside the runtime wheel.

## Constraints

- `benchmarks/` must remain repo-local and must not be included in the built
  wheel.
- The relocation is structural only. Benchmark behavior should not change.
- No compatibility shim should remain under `SpliceGrapher/formats/`.
- This slice does not include `SpliceGraph.py` shim work or any benchmark API
  redesign.

## Current State

- The wheel target in `pyproject.toml` includes only `SpliceGrapher`, so a new
  top-level `benchmarks/` package will remain out of the shipped distribution by
  default.
- Runner and test code currently import benchmark helpers from
  `SpliceGrapher.formats.polars_gff_benchmark`.
- A top-level `scripts/benchmarks/` area already exists for CLI runners, so the
  benchmark module should move beside the runtime package rather than inside the
  `scripts/` tree.

## Chosen Approach

Create a top-level `benchmarks/` package and move the benchmark module there:

- `benchmarks/__init__.py`
- `benchmarks/polars_gff_benchmark.py`

Then update runner and test imports to consume:

- `benchmarks.polars_gff_benchmark`

The old runtime module will be deleted outright. No wrapper module or re-export
will remain under `SpliceGrapher/formats/`.

## Rejected Alternatives

### Keep it under `scripts/benchmarks/`

Rejected because the module is imported by tests and runners. Putting it under
`scripts/` would turn an importable helper module into a script-path problem and
invite `sys.path` hacks.

### Leave a shim under `SpliceGrapher/formats/`

Rejected because it preserves runtime pollution and weakens the exact boundary
this slice is supposed to establish.

### Split the benchmark module while relocating it

Rejected for this slice because it would mix structural cleanup with internal
benchmark redesign and produce a larger, noisier PR.

## Files Affected

- Create: `benchmarks/__init__.py`
- Create: `benchmarks/polars_gff_benchmark.py`
- Modify: `scripts/benchmarks/run_gff_loader_benchmarks.py`
- Modify: `tests/test_polars_gff_benchmark.py`
- Modify: `tests/test_polars_gff_benchmark_runner.py`
- Delete: `SpliceGrapher/formats/polars_gff_benchmark.py`

## Verification Strategy

- Run the benchmark-focused pytest slice:
  - `tests/test_polars_gff_benchmark.py`
  - `tests/test_polars_gff_benchmark_runner.py`
- Add or tighten one assertion that the wheel target still points only at
  `SpliceGrapher`.
- Run build verification so the packaging boundary is checked from the actual
  build config, not inferred.
