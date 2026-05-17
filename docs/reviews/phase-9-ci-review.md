# Phase 9 Review: GitHub Actions CI

Reviewer: Carver

Scope:

- `.github/workflows/ci.yml`
- `tests/test_ci_workflow.py`
- `README.md`

Findings addressed:

- The CI contract test checked commands but did not lock the Python version
  matrix or uv cache behavior. It now asserts Python 3.11/3.12, uv cache
  enablement, and `uv.lock` cache dependency configuration.

Verification:

```bash
uv run --extra dev pytest tests\test_ci_workflow.py tests\test_three_api_docs.py
uv run --extra dev pytest
uv run --extra dev ruff check .
uv run --extra dev ruff format --check .
uv run --extra dev ty check
uv build
```

Residual risk:

- CI is Ubuntu-only and does not yet run real external microbiome tools, Docker
  builds, or Nextflow workflows. Those should become separate jobs when the
  container/workflow layer is ready for execution.
