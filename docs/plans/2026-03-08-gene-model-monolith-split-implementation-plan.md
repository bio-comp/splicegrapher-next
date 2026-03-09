# Gene-Model Monolith Split Implementation Plan

Date: 2026-03-08
Issue: #175
Branch: `sg-next/175-genemodel-monolith-plan`

## Phase 0. Inventory and baseline

- Confirm direct import surface for `formats.models` and `parsers.gene_model_gff`
- Keep `gene_model.py` as the orchestration façade
- Keep runtime behavior unchanged during the first split

## Phase 1. Split `SpliceGrapher/formats/models.py`

### Step 1. Create `model_index.py`
Move:
- `IntervalQuery`
- `ChromosomeGeneIndex`
- `Chromosome`
- `gene_type_filter`
- `default_gene_filter`
- `feature_cmp`
- `feature_sort_key`
- `gene_sort_key`
- `gtf_feature_sort_key`
- `feature_overlaps`
- `feature_contains`
- supporting interval protocols

### Step 2. Create `model_domain.py`
Move:
- constants and aliases used by the domain layer
- `cds_factory`
- `Locus`
- `BaseFeature`
- `TranscriptRegion`
- `Exon`, `CDS`, `FpUtr`, `TpUtr`
- `Transcript`
- `Gene`
- `PseudoGene`

### Step 3. Reduce `models.py` to an aggregator
- re-export the moved symbols from `model_domain.py` and `model_index.py`
- keep existing import paths working

### Verification
Run:
- `tests/test_gene_model.py`
- `tests/test_gene_model_serializers.py`
- `tests/test_gene_model_gff_parser.py`
- `tests/test_integration_simple.py`
- targeted `ruff` and `mypy`

## Phase 2. Split `SpliceGrapher/formats/parsers/gene_model_gff.py`

### Step 1. Create `gene_model_gff_context.py`
Move:
- `GeneModelLike`
- `ParseStats`
- `ParsedRecord`
- `ParseContext`

### Step 2. Create `gene_model_gff_records.py`
Move:
- normalization helpers
- annotation parsing helpers
- pending-child queue helpers
- record-line parsing helpers

### Step 3. Create `gene_model_gff_handlers.py`
Move:
- parent resolution helpers
- gene/mrna/exon/transcript-region/misc/chromosome handlers

### Step 4. Reduce `gene_model_gff.py` to orchestration
Keep:
- `load_gene_model_records`
- top-level record-dispatch flow
- imports from the new parser helper modules

### Verification
Run:
- `tests/test_gene_model_gff_parser.py`
- `tests/test_gene_model.py`
- `tests/test_gene_model_serializers.py`
- targeted `ruff` and `mypy`

## Commit Plan

1. `docs: add gene-model monolith split design and plan`
2. `refactor(formats): split gene-model domain and index modules`
3. `refactor(formats): split gene-model gff parser internals`

## Guardrails

- Do not stage local `AGENTS.md` or `LOCAL_TODO.md`
- Do not widen public API surface during the split
- Prefer aggregator/facade stabilization over another hard-cut import break unless tests show the façade itself is the problem
