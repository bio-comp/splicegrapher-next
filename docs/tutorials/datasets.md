# Tutorial Dataset Curation (Plant + Human)

This document defines the canonical modern datasets used by
`splicegrapher-next` tutorial notebooks and e2e/integration tracks.

Canonical selection is pinned in `docs/tutorials/dataset_manifest.toml`.
Species entries use short stable keys (`ath`, `hsa`) and store scientific names
and taxonomy IDs in explicit fields.

## Why These Datasets

Goals for canonical tutorial/e2e data:

- modern, publication-backed provenance
- clear relevance to transcript and alternative splicing analysis
- tractable subset size for local and CI execution
- stable accession metadata for reproducible acquisition

Selected canonical studies:

- Plant: Arabidopsis AtRTD3 study (`PRJNA755474`)
- Human: DTU benchmark context (`GSE172421`, short-read subset)

## Tiering Model

Two-tier model:

- Notebook tier:
  - source data fetched from canonical public studies
  - suitable for local hands-on tutorials
- CI/e2e tier:
  - tiny derived fixtures generated from canonical studies
  - tracked under `tests/fixtures/` with checksums

CI should not depend on full public data downloads.

## Acquisition Methods

Both methods are supported by policy:

- Primary: Python/HTTP metadata and file acquisition (implemented with `httpx`)
- Alternative: SRA Toolkit (`prefetch`, `fasterq-dump`) for users who prefer it

### Python/HTTP Method (Default)

Planned script contract:

```bash
uv run python scripts/data/fetch_tutorial_data.py \
  --manifest docs/tutorials/dataset_manifest.toml \
  --species plant \
  --method python-httpx
```

### SRA Toolkit Method (Alternative)

Example command shape:

```bash
prefetch SRR17001820 SRR17001821
fasterq-dump SRR17001820 SRR17001821 --threads 8 --outdir data/raw/plant
```

Human example:

```bash
prefetch SRR14684074 SRR14684079
fasterq-dump SRR14684074 SRR14684079 --threads 8 --outdir data/raw/human
```

## Reproducibility Requirements

- Keep accession IDs pinned in `dataset_manifest.toml`.
- Keep study references and citation metadata updated when dataset choices
  change.
- Add and verify checksums for derived fixture bundles.
- If acquisition logic changes, document it in the related issue/PR.
- Any use of sync HTTP or `requests` requires explicit issue/PR justification
  under repository HTTP policy.

## Curation Update Rules

When rotating canonical runs:

1. Open/update a scoped issue.
2. Update `dataset_manifest.toml` with explicit accession changes.
3. Rebuild fixture tier and refresh checksums.
4. Verify notebook smoke execution still passes.
5. Call out migration/compatibility impact in PR notes.

## References

- AtRTD3 (Arabidopsis) data availability and publication:
  - <https://doi.org/10.1186/s13059-022-02711-0>
- DTU benchmark publication:
  - <https://doi.org/10.1038/s41592-023-02026-3>
- GEO study context for human benchmark:
  - <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172421>
