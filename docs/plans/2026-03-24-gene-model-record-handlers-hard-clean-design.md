# Gene Model Record Handlers Hard Clean Design

## Goal

Reduce branch density and duplicated dispatch logic in
`SpliceGrapher/formats/parsers/gene_model_gff_record_handlers.py` while
preserving parser behavior and public API boundaries.

## Existing Smells

- Repeated transcript-parent-gene lookup logic across exon/CDS-region handlers.
- Repeated model-link registration (`mrna_forms`/`mrna_gene`) plus pending-drain flow.
- Mixed dispatch and bookkeeping concerns inside large handler functions.

## Proposed Shape

Introduce internal helper seams inside the handlers module:

1. `_parent_gene_for_transcript(...)`
- Centralize transcript -> parent gene lookup with fail-fast error handling.

2. `_register_transcript_links(...)`
- Centralize transcript registration to model maps and pending child draining.

3. `_resolve_region_target(...)`
- Centralize CDS/UTR parent transcript resolution including pending queue behavior.

Handlers retain current signatures and routing through `RECORD_HANDLERS`.

## Behavioral Constraints

- No changes to public parser entrypoints (`load_gene_model_records`).
- No new compatibility wrappers or alias APIs.
- Preserve existing fail/queue semantics and verbose messaging intent.

## Verification Strategy

- New boundary tests for helper seams and queue behavior.
- Existing parser and gene model tests as regression net.
- Full repo gates before PR.
