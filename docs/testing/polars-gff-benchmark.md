# Polars GFF Single-Cycle Go/No-Go Evaluation

Issue: #94 (follow-up from #62)

## Scope

This benchmark implements SGN's timeboxed go/no-go policy for optional backend paths.

Compared paths:

- baseline object path: `SpliceGrapher.formats.gene_model.GeneModel`
- baseline row path: `SpliceGrapher.formats.polars_gff.load_gff_rows`
- candidate path: `SpliceGrapher.formats.polars_gff.load_gff_to_polars`

Measured workload classes:

- ingest microbenchmark
- end-to-end analytics benchmark with shared exon grouping, transcript/exon aggregation,
  and donor/acceptor-style coordinate derivations

## Decision Rule

Polars is adoptable only if all are true:

- >=20% runtime win on at least 2 of 3 real datasets
- <=10% peak memory regression on winning datasets
- no correctness drift against baseline analytics signatures

If not met, recommendation is `defer` and Polars remains optional analysis helper only.

## Methodology

- One warm-up run per loader/workload before measurement.
- 3 timed iterations per loader/workload.
- Metrics captured:
  - mean and max wall-clock seconds
  - peak traced Python heap MiB (`tracemalloc`)
  - analytics signature for end-to-end correctness checks
- Synthetic defaults:
  - `small`: 100 genes
  - `medium`: 750 genes
  - `large`: 3000 genes

## Real Dataset Set (This Cycle)

From local idiffir SpliceGrapher corpus:

- `at5`: `iDiffIR/SpliceGrapher/examples/AT5G03770_model.gff`
- `at1`: `iDiffIR/SpliceGrapher/examples/AT1G13440_model.gff`
- `at2`: `iDiffIR/SpliceGrapher/tutorial/AT2G04700.gff`

## Reproducible Command

From repository root:

```bash
uv run --with polars python scripts/benchmarks/run_gff_loader_benchmarks.py \
  --iterations 3 \
  --synthetic-work-dir .benchmarks/gff \
  --real-dataset at5=/Users/mikeh/repos/idiffir/iDiffIR/SpliceGrapher/examples/AT5G03770_model.gff \
  --real-dataset at1=/Users/mikeh/repos/idiffir/iDiffIR/SpliceGrapher/examples/AT1G13440_model.gff \
  --real-dataset at2=/Users/mikeh/repos/idiffir/iDiffIR/SpliceGrapher/tutorial/AT2G04700.gff
```

Generated artifacts:

- `docs/testing/polars_gff_benchmark_results.json`
- `docs/testing/polars_gff_benchmark_results.md`

## Result (Single Cycle)

Current decision: **`defer (threshold_not_met)`**

- evaluated real datasets: 3
- compared real datasets: 3
- winning datasets: 0

Observed notes:

- Polars did not hit the >=20% end-to-end runtime win threshold on required real datasets.
- For `at2`, `GeneModel` reports `No gene models found`; row and Polars paths still benchmark successfully and are compared for go/no-go decision.
- No default path migration is recommended.

## Recommendation

Keep Polars path as optional helper only.

- Keep `SpliceGrapher.formats.polars_gff` for analysis/experimentation.
- Do not switch SGN default GFF loader/model path based on this cycle.
- Re-open adoption only if a future cycle meets runtime/memory/correctness gates.
