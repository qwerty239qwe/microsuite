# Phase 2 Review

Reviewer: spawned review agent for the differential abundance cleanup.

Accepted findings and fixes:

- Added a subprocess test proving the ANCOM-BC wrapper invokes
  `scripts/r/ancombc.R` with the expected command arguments.
- Added a compatibility test for the legacy `microsuite diffab ancombc` command.
- Renamed backend-loop variables in the method catalog from `method` to
  `backend`.
- Added CLI documentation noting that `diffab ancombc` remains as a deprecated
  compatibility command for `diff_abundance --backend ancombc`.

Verification after fixes:

```text
uv run --extra dev pytest
uv run --extra dev ruff check .
uv run --extra dev ruff format --check .
uv run --extra dev ty check
```
