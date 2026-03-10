# Gene-Model Domain Split Implementation Plan

1. Add a module-layout test that pins the new internal package shape.
2. Extract shared constants/factories from `domain.py` into `constants.py`.
3. Extract `Locus` into `locus.py`.
4. Extract `BaseFeature`, `TranscriptRegion`, `Exon`, `CDS`, `FpUtr`, and `TpUtr` into `features.py`.
5. Extract `Transcript` into `transcript.py`.
6. Extract `Gene` and `PseudoGene` into `gene.py`.
7. Replace `domain.py` imports in `SpliceGrapher/formats/models/__init__.py` with explicit re-exports from the new modules.
8. Update internal imports to use package-relative boundaries and resolve circular typing with `TYPE_CHECKING` where needed.
9. Run focused gene-model tests, then full repo verification.
10. Commit the structural split as one reviewable refactor tranche.
