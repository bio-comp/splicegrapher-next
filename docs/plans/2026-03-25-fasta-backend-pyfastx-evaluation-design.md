# FASTA Backend Pyfastx Evaluation Design

## Goal

Evaluate `pyfastx` as a FASTA backend candidate without changing SGN runtime behavior in this slice.

## Constraints

- Behavior-neutral for production FASTA APIs in `SpliceGrapher.formats.fasta`.
- No dependency changes in `pyproject.toml`.
- Produce a recommendation with concrete parity/performance/index evidence.

## Approach

1. Add a repo-local evaluation probe in `benchmarks/`:
- compare SGN current behavior against `pyfaidx` and optional `pyfastx`
- measure iteration + random-access timing, plus index-build overhead
- check contract compatibility (notably open file handle support)

2. Persist one benchmark artifact under `docs/testing/` from a reproducible command:
- `uv run --with pyfastx python benchmarks/fasta_backend_probe.py ...`

3. Record decision in ADR:
- adopt / optional / reject with compatibility implications for SGN/iDiffIR/TAPIS

## Recommendation Target

Use evidence to decide whether pyfastx should be default, optional, or rejected.
Default adoption requires preserving SGN FASTA contracts and migration safety.
