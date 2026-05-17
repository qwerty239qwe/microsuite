# Phase 7 Review: microsuite Rename

Reviewer: Herschel

Scope:

- `pyproject.toml`
- `src/microsuite`
- `tests`
- `README.md`
- `docs`
- `containers/microsuite`

Findings addressed or accepted:

- Stale generated run provenance under `runs/` still used the old toolbox name.
  The `runs/` directory is ignored and is not part of the committed source.
- CLI compatibility alias exists for `microbiome`, but Python import
  compatibility for `microbiome_suite` does not. Accepted as a clean pre-0.1
  rename; `microsuite` is the package import going forward.

Verification:

```bash
uv run --extra dev pytest
uv run --extra dev ruff check .
uv run --extra dev ruff format --check .
uv run --extra dev ty check
uv run --extra dev microsuite --help
uv run --extra dev microbiome --help
uv build
```

Residual risk:

- Docker image builds were not run. Container files were statically checked by
  tests and review.
