# Phase 5 Review: Reporting Layer

Reviewer: Hypatia

Scope:

- `src/microsuite/methods/report.py`
- `src/microsuite/cli/method_cmd.py`
- `tests/test_report_method.py`
- `README.md`
- `docs/api-cli.md`
- `docs/toolbox.md`
- `docs/three-api-roadmap.md`

Findings addressed:

- Roadmap acceptance overclaimed report depth. It now describes the implemented
  provenance summary: `run.json` plus optional `outputs.json`.
- Tests did not prove `outputs.json` was read. The synthetic run now records a
  distinct output only in `outputs.json` and asserts it appears in the HTML.
- Malformed JSON surfaced as `json.JSONDecodeError`. Report inputs now raise
  `MicrobiomeSuiteError` for invalid JSON and non-object JSON.

Verification:

```bash
uv run --extra dev pytest tests\test_report_method.py
uv run --extra dev ruff check src\microsuite\methods\report.py tests\test_report_method.py
uv run --extra dev ty check
```

Residual risk:

- The native report is intentionally a metadata/provenance summary. TSV
  summaries, figure embedding, and formal Nextflow run-folder validation remain
  deferred.
