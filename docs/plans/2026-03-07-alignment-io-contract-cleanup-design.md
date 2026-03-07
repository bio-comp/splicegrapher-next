# alignment_io.py Hard-Break Rewrite Design

## Goal

Replace the legacy camelCase `alignment_io.py` public API with a snake_case-only alignment I/O surface, then repair the SGN call sites and tests against that new boundary.

## Scope

This tranche deliberately breaks the old API.

In scope:
- delete camelCase public entrypoints from `SpliceGrapher/formats/alignment_io.py`
- replace them with explicit snake_case functions
- keep the corrected tuple-shape contract in the new `collect_alignment_data(...)` API
- update SGN call sites and tests to the new names
- tighten touched typing and keep the module internally compartmentalized by concern

Out of scope:
- compatibility aliases in `alignment_io.py`
- benchmark relocation
- `SpliceGraph.py` lowercase/shim migration
- downstream iDiffIR/TAPIS repair
- broad parser/package extraction outside the touched alignment surface

## Replacement API

Delete:
- `getSamReadData`
- `getSamDepths`
- `getSamJunctions`
- `getSamAlignments`
- `getSamHeaders`
- `getSamHeaderInfo`
- `getSamSequences`
- `pysamReadDepths`
- `pysamStrand`
- `isBamFile`
- `isCramFile`
- `makeChromosomeSet`

Replace with:
- `collect_alignment_data`
- `read_alignment_depths`
- `read_alignment_junctions`
- `read_alignment_spans`
- `read_alignment_headers`
- `read_alignment_chromosome_info`
- `read_alignment_sequences`
- `calculate_gene_depths`

Private helpers will be snake_case and private by default, for example:
- `_is_bam`
- `_is_cram`
- `_make_chromosome_set`
- `_record_strand`

## Contract Decisions

### 1. No compatibility wrappers

This branch does not preserve the old names. The repo should fail loudly anywhere that still imports or calls them. Those failures are the migration map.

### 2. Keep the modern tuple contract

`collect_alignment_data(..., include_alignments=True)` returns a 3-tuple on both branches:
- depths map
- junction map
- alignment map

Depths-file fallback returns an empty alignment map instead of collapsing to a 2-tuple.

### 3. Keep the module split by concern

The file remains internally grouped into:
1. type aliases/protocols/constants
2. source normalization and opening helpers
3. pysam-backed collection helpers
4. depths-file fallback bridge
5. public snake_case API

This keeps the rewrite reviewable and prepares future file splits without combining that work into `#165`.

## Testing Strategy

1. Rewrite alignment I/O tests to the new snake_case API first.
2. Run the focused alignment slice and observe the import/call failures.
3. Rename the public functions in `alignment_io.py` and delete the camelCase wrappers.
4. Sweep SGN tests/callers to the new API.
5. Verify the touched alignment surface with Ruff, MyPy, and targeted pytest.

## Expected Blast Radius

The initial failures should come from:
- `tests/test_alignment_io_parity.py`
- `tests/test_alignment_io_process_utils_boundary.py`
- `tests/test_splicegrapher_alignment_io.py`
- any SGN callers/imports of `getSamReadData`, `getSamDepths`, `getSamJunctions`, `getSamAlignments`, or header helpers

That blast radius is intentional.

## Follow-up Work (Not This Tranche)

1. Move `polars_gff_benchmark.py` out of the runtime package.
2. Migrate `SpliceGraph.py` to a lowercase compatibility shim model.
3. Split `alignment_io.py` into smaller modules after the new API has settled.
