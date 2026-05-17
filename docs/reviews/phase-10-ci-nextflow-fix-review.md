# Phase 10 Review: CI Dependency Fix and Nextflow Smoke

Reviewer: Descartes

Scope:

- `.github/workflows/ci.yml`
- `tests/test_ci_workflow.py`
- `.gitignore`
- Local Docker-based Nextflow smoke command

Findings addressed:

- Local Nextflow smoke runs can leave published outputs under `results/`.
  Added `results/` to `.gitignore`.

Verification:

```bash
uv sync --all-extras --locked
uv run --extra dev pytest tests\test_ci_workflow.py tests\test_nextflow_skeleton.py
uv run --extra dev pytest
uv run --extra dev ruff check .
uv run --extra dev ruff format --check .
uv run --extra dev ty check
uv build
docker run --rm -v ${PWD}:/workspace -w /workspace nextflow/nextflow:25.10.0 \
  nextflow run workflows/nextflow/main.nf \
  -profile local \
  --manifest tests/fixtures/moving_pictures_small/table.tsv \
  --metadata tests/fixtures/moving_pictures_small/metadata.tsv \
  --classifier tests/fixtures/moving_pictures_small/taxonomy.tsv \
  --outdir /workspace/runs/nextflow-smoke \
  -work-dir /tmp/microsuite-nextflow-work
```

Residual risk:

- The Nextflow smoke job exercises placeholder modules only. Real external-tool,
  Docker, and Singularity profile smoke tests should be added when those modules
  stop being placeholders.
