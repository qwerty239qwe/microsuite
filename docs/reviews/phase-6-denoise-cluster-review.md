# Phase 6 Review: Denoise and Cluster Methods

Reviewer: Euclid

Scope:

- `src/microsuite/methods/denoise.py`
- `src/microsuite/methods/cluster.py`
- `src/microsuite/cli/method_cmd.py`
- `tests/test_denoise_cluster_methods.py`
- `README.md`
- `docs/api-cli.md`
- `docs/toolbox.md`

Findings addressed:

- Paired-end DADA2 command construction was untested. Added coverage for the
  `denoise-paired` branch and paired trimming/truncation flags.
- Deblur positive command construction was untested. Added coverage for the
  `qiime deblur denoise-16S` command.
- Cluster CLI path was not exercised. Added a CLI failure test for missing QIIME
  on `microsuite cluster --backend vsearch`.
- README grouped planned `dada2-r` with implemented denoise backends. The table
  now distinguishes implemented QIIME 2 backends from planned `dada2-r`.

Verification:

```bash
uv run --extra dev pytest tests\test_denoise_cluster_methods.py
uv run --extra dev pytest
uv run --extra dev ruff check .
uv run --extra dev ruff format --check .
uv run --extra dev ty check
```

Residual risk:

- These wrappers are statically tested with monkeypatched subprocess calls. A
  real QIIME 2/container smoke test is still needed before treating raw-read
  workflows as production-ready.
