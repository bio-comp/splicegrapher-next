# Modern Assay Manifesto

`splicegrapher-next` is a generally usable Python 3 modernization/continuation
of SpliceGrapher, with compatibility priorities for iDiffIR and TAPIS during
active migration.

This manifesto defines how modern assay support is added without breaking
legacy workflows or overpromising unsupported biology.

## Positioning Contract

1. `splicegrapher-next` is not iDiffIR-only or TAPIS-only.
2. iDiffIR and TAPIS compatibility remain explicit migration priorities.
3. Public identity and language must match:
   - `README.md`
   - `pyproject.toml`
   - `NOTICE`
   - `PROVENANCE.md`

## Architecture Direction

1. Compatibility core first:
   - preserve `SpliceGrapher` import compatibility
   - keep behavior parity before deep rewrites
2. Evidence abstraction second:
   - normalize assay evidence into typed interval/evidence contracts
   - keep parser boundaries fail-fast and typed
3. Visualization remains first-class:
   - `SpliceGrapher/view/*` and `SpliceGrapher/plot/*` stay on the critical path

## Modern Assay Scope

Supported incrementally via adapters and aggregation layers:

1. Bulk short-read RNA-seq
2. Grouped/pseudobulk single-cell summaries
3. Bulk long-read evidence paths

Future/optional expansions:

1. Per-cell and protocol-heavy workflows through plugin-style adapters
2. AnnData/interoperability integrations that do not destabilize core APIs

## Roadmap Links

Current execution and hardening tracks:

1. Identity/docs alignment (this issue): `#17`
2. Overlap architecture umbrella: `#103`
3. Cleanup/type tracks:
   - `#100` snake_case modernization path
   - `#101` remove `typing.Any` from shared utility surfaces
   - `#102` broaden `GeneModel` type coverage

Milestones:

1. SGN-M2 Stabilize+Adopt:
   - <https://github.com/bio-comp/splicegrapher-next/milestone/2>
2. SGN-M3 Clean-Invariant:
   - <https://github.com/bio-comp/splicegrapher-next/milestone/3>

## Non-Goals

1. No big-bang graph-engine rewrite.
2. No silent runtime/API behavior changes.
3. No dependency expansion without issue-scoped justification.

## Engineering Rules

1. Keep issue -> branch -> PR isolation.
2. Land parity tests with structural refactors.
3. Keep clean-invariant ratchet monotonic.
4. Document migration impact explicitly in PR summaries.
