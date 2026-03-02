# Interval Overlap Benchmark Gate

This benchmark tracks parity and performance for interval-overlap helper paths.
It is designed as a lightweight CI guardrail and a reproducible local check.

## Default policy

- parity must hold between legacy overlap scan and indexed overlap helper results
- indexed slowdown ratio (`indexed_mean / legacy_mean`) must be `<= 1.35`

The ratio threshold intentionally allows some runtime variance across CI runners
while still blocking meaningful regressions.

## Reproducible command

```bash
uv run python scripts/benchmarks/run_interval_overlap_benchmarks.py \
  --interval-count 2000 \
  --query-count 2000 \
  --iterations 3 \
  --assert-threshold \
  --max-slowdown-ratio 1.35 \
  --json-out docs/testing/interval_overlap_benchmark_results.json \
  --markdown-out docs/testing/interval_overlap_benchmark_results.md
```

## CI gate command

The quality-gates workflow runs the same benchmark policy with CI-friendly
runtime parameters and fails on threshold breaches:

```bash
uv run python scripts/benchmarks/run_interval_overlap_benchmarks.py \
  --interval-count 1200 \
  --query-count 1200 \
  --iterations 2 \
  --assert-threshold \
  --max-slowdown-ratio 1.35 \
  --json-out /tmp/interval_overlap_benchmark_results.json \
  --markdown-out /tmp/interval_overlap_benchmark_results.md
```

## Artifacts

- JSON result payload: metrics + threshold fields
- Markdown report: human-readable status summary
