# SGN Enum Harmonization Design

**Date:** 2026-03-01
**Status:** Approved
**Decision Mode:** Aggressive staged rollout with fail-fast semantics and willingness to accept breaking changes when they materially simplify architecture.

## Context

`splicegrapher-next` currently relies on widespread ad-hoc string constants in core control-flow logic (notably `SpliceGrapher/SpliceGraph.py` and `SpliceGrapher/formats/GeneModel.py`). This creates duplicate domains, weak validation, and brittle behavior around parser and serialization boundaries.

The project direction is to harmonize finite domains via enums and use them throughout SGN, prioritizing architectural clarity and rapid surfacing of integration defects.

## Goals

1. Replace magic-string control flow with canonical enums across SGN core modules.
2. Centralize finite-domain semantics (record types, strands, node dispositions, edge types, attribute keys).
3. Fail hard and fast in internal logic when invalid enum values are encountered.
4. Preserve external file format behavior where required (`GFF/GTF` text output remains contract-stable unless intentionally changed and documented).
5. Make migration progress explicit and testable through tranche-scoped issues and parity guards.

## Non-Goals

1. Full legacy API compatibility for all internal constants.
2. One-shot perfect migration of every module in a single commit.
3. Silent coercion of unknown values in internal core paths.

## Approved Architecture

### 1) Canonical enum module

Create a shared enum module at:
- `SpliceGrapher/core/enums.py`

Initial enum set:
- `Strand`
- `RecordType`
- `NodeDisposition`
- `EdgeType`
- `AttrKey`

Use `str`-backed enums (`class X(str, Enum)`) so boundary serialization remains straightforward (`.value`) and integrations expecting strings remain manageable at adapters.

### 2) Boundary coercion + validation

Add narrow conversion helpers for parser and input boundaries (e.g. `coerce_record_type`, `coerce_strand`).

Rules:
- Internal control flow uses enums, not raw strings.
- Unknown/invalid values raise explicit `ValueError` in internal logic.
- Boundary adapters may normalize legacy input where explicitly intended.

### 3) Serialization contract handling

At output boundaries (`gffString`, `gtfString`, related formatter paths), use enum `.value` explicitly to emit text contracts.

### 4) Fast-fail posture

Migration defaults to strict failures for invalid states in touched code paths. This is intentional to surface breakpoints quickly and reveal downstream assumptions.

## Alternatives Considered

1. Policy-only tracking: rejected as insufficient; does not reduce current ambiguity.
2. Hybrid with broad compatibility shims: rejected for now because it masks integration defects.
3. Global codemod in one pass: accepted as guiding posture, implemented in tranches for reviewability.

## Tranche Plan

### E0: Foundation

- Add `SpliceGrapher/core/enums.py`.
- Add coercion/validation helpers.
- Add baseline enum unit tests.

### E1: SpliceGraph-first

- Migrate `SpliceGrapher/SpliceGraph.py` control-flow domains to enums.
- Remove duplicated string constants that become obsolete.
- Preserve serialized output via `.value`.

### E2: GeneModel + formats

- Migrate `SpliceGrapher/formats/GeneModel.py` domains to enums.
- Continue into adjacent format modules touching record/attribute/strand handling.

### E3: shared + scripts/CLI

- Migrate remaining shared modules and script/CLI control-flow domains to enums.
- Add/extend boundary adapters where external inputs still arrive as strings.

### E4: Hardening

- Remove temporary migration shims/adapters that are no longer needed.
- Enforce no-magic-string policy in touched modules.

## Verification Strategy

Per tranche:
- Add targeted tests before behavior changes where risk is nontrivial.
- Validate parser coercion success/failure behavior.
- Verify serialization parity where contracts are stable.

Global gates:
- `uv run ruff check .`
- `uv run ruff format --check .`
- `uv run pytest -q`
- `uv run python scripts/ci/check_clean_invariant.py`

Guardrail extension:
- Add banned-pattern checks for magic-string control flow in touched modules (allowlisted for explicit boundary adapters only).

## Risk and Rollback

Risks:
- Intentional breakage of undocumented internal assumptions.
- Downstream callers relying on legacy constants.

Mitigation:
- Keep tranches and commits small and reviewable.
- Include explicit migration notes in PRs.
- Revert by tranche if instability is detected; keep enum foundation where safe.

## Tracking Updates Required

1. Add explicit enum harmonization policy to `AGENTS.md`.
2. Add enum epic and tranche checklist to `LOCAL_TODO.md`.
3. Open umbrella GitHub issue for enum harmonization and link tranche issues.

## Exit Criteria

1. Core control-flow domains in `SpliceGraph.py` and `GeneModel.py` are enum-backed.
2. New touched modules do not introduce ad-hoc string domain constants.
3. CI and clean-invariant gates pass with expanded anti-magic-string guardrails.
4. Migration impact and remaining debt are explicitly tracked.
