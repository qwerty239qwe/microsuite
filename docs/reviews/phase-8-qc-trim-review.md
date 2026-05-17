# Phase 8 Review: QC and Trim Methods

Reviewer: Kuhn

Scope:

- `src/microsuite/methods/qc.py`
- `src/microsuite/methods/trim.py`
- `src/microsuite/cli/method_cmd.py`
- `tests/test_qc_trim_methods.py`
- `README.md`
- `docs/api-cli.md`
- `docs/toolbox.md`
- `docs/three-api-roadmap.md`

Findings addressed:

- `--force` for `qiime2-demux` allowed an existing `.qzv` through local
  validation but would still leave a target that QIIME 2 could reject. Existing
  output files are now unlinked when `force=True`.
- `multiqc --force` was not propagated to the external command. It is now added
  when `force=True`.
- CLI coverage only exercised `trim`. Added CLI-level `qc --backend fastqc`
  missing-tool coverage.

Verification:

```bash
uv run --extra dev pytest tests\test_qc_trim_methods.py
uv run --extra dev pytest
uv run --extra dev ruff check .
uv run --extra dev ruff format --check .
uv run --extra dev ty check
```

Residual risk:

- External tool behavior is statically tested with monkeypatched subprocesses.
  Real FastQC, MultiQC, QIIME 2, and fastp smoke tests should be added once
  containers are ready.
