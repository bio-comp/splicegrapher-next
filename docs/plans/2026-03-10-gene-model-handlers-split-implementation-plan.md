# Gene-Model Parser Handlers Split Implementation Plan

1. Add/adjust parser layout tests to pin the new internal handler modules.
2. Extract pending-child queue helpers from `gene_model_gff_handlers.py` into `gene_model_gff_pending.py`.
3. Extract parent / isoform / strand resolution helpers into `gene_model_gff_resolution.py`.
4. Extract record-type handler functions into `gene_model_gff_record_handlers.py`.
5. Replace imports in `gene_model_gff.py` with the new handler module boundary.
6. Remove `gene_model_gff_handlers.py`.
7. Run focused parser/gene-model verification, then full repo gates.
8. Commit the structural split as one reviewable refactor tranche.
