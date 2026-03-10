# Shared Process Hard Cut Design

## Summary

Hard-clean the shared helper layer by deleting `SpliceGrapher/shared/process_utils.py`, replacing it with a modern `SpliceGrapher/shared/process.py` module, and removing the remaining Python 2 iterator residue from `SpliceGrapher/shared/progress.py`.

## Goals

- Remove the legacy shared process helper namespace pollution.
- Replace camelCase shared process APIs with snake_case names and typed interfaces.
- Update all SGN callers and tests in the same PR, with no compatibility aliases.
- Remove `next = __next__` and other touched legacy residue from `shared/progress.py`.
- Keep runtime behavior stable where tests already pin it.

## Non-Goals

- No broader reorganization of the `shared/` package beyond this module rename.
- No compatibility wrappers for deleted names like `runCommand`, `logMessage`, or `getAttribute`.
- No redesign of progress behavior beyond the necessary cleanup and caller/test updates.

## Current Problems

`SpliceGrapher/shared/process_utils.py` still carries explicit legacy debt:

- camelCase API names (`getAttribute`, `logMessage`, `runCommand`, `writeStartupMessage`)
- `**args: object` typing smell
- direct `sys.stderr.write(...)` side effects in touched runtime code

`SpliceGrapher/shared/progress.py` is cleaner, but still contains a Python 2 compatibility alias:

- `next = __next__`

These modules sit in the shared layer, so leaving them untouched keeps legacy style alive in code that other parts of SGN still import.

## Proposed Module Layout

```text
SpliceGrapher/shared/
├── process.py          # new modern process helper boundary
└── progress.py         # cleaned in place
```

Delete:

- `SpliceGrapher/shared/process_utils.py`

## Public API After the Cut

`SpliceGrapher.shared.process`
- `get_attribute`
- `id_factory`
- `log_message`
- `run_command`
- `run_logged_command`
- `write_startup_message`

`SpliceGrapher.shared.progress`
- `ProgressIndicator`
- `RandomListIterator`
- no `next = __next__`

Deleted names:
- `getAttribute`
- `idFactory`
- `logMessage`
- `runCommand`
- `writeStartupMessage`

## Design Details

### Process Module

`run_command` remains the low-level subprocess wrapper.

`run_logged_command` replaces the old high-level `runCommand` behavior:
- emits structured start/failure logs
- mirrors messages to an optional `logstream`
- defaults stdout/stderr to `DEVNULL` when not explicitly provided
- raises on non-zero return when configured

`get_attribute` should stop using variadic keyword capture and instead accept a typed mapping input:

```python
def get_attribute(mapping: Mapping[str, object], key: str, default: AttrT) -> AttrT:
    ...
```

`log_message` and `write_startup_message` should stop writing directly to `stderr` and instead use the existing logger plus optional stream mirroring where needed.

### Progress Module

Keep the file in place and modernize only the touched residue:
- remove `next = __next__`
- keep the typed iterator implementation intact
- preserve `ProgressIndicator` behavior because it is a live runtime dependency for parsers and writers

## Testing Strategy

Tests need a hard-cut update in the same PR:

- `tests/test_process_utils.py` becomes `tests/test_process.py`
- smoke/invariant tests should import `SpliceGrapher.shared.process`
- delete or rewrite any assertions that exist only to pin deleted legacy names
- keep targeted progress tests to ensure `ProgressIndicator` behavior remains stable after cleanup

Verification gates for the slice:
- targeted shared-helper tests
- full `ruff`, `pytest`, `mypy SpliceGrapher tests`, and `uv build`

## Migration Impact

SGN internal callers and tests must move to the new module path and snake_case names in the same PR.

Downstream repos using the deleted legacy shared helper names will break and must be updated explicitly. This is acceptable because the goal of this slice is to remove shared-layer legacy pollution rather than preserve compatibility wrappers.
