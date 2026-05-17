# Phase 1 Review

Reviewer: spawned review agent for the Phase 1 API documentation and SDK facade.

Accepted findings and fixes:

- `docs/toolbox.md` described the CLI as the only public surface. Updated it to
  state that the CLI is one public surface alongside Nextflow and the Python SDK.
- `docs/api-nextflow.md` showed a concrete `nextflow run` command before the
  Nextflow files exist. Added wording that the entry point is planned until
  Phase 4.
- `tests/test_three_api_docs.py` was too weak to catch boundary-language drift.
  Added targeted assertions for the CLI boundary, planned Nextflow status, and
  SDK `.h5ad` limitation.
- `docs/api-python.md` omitted the current `.h5ad` limitation for `read_table`
  and `write_table`. Added the constraint explicitly.

Verification after fixes:

```text
uv run --extra dev pytest
uv run --extra dev ruff check .
uv run --extra dev ruff format --check .
uv run --extra dev ty check
```
