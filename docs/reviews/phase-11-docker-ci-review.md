# Phase 11 Review: Docker CI Builds

Reviewer: Bohr

Scope:

- `.github/workflows/ci.yml`
- `.dockerignore`
- `containers/kraken2/Dockerfile`
- `docs/containers.md`
- `tests/test_ci_workflow.py`
- `tests/test_container_skeletons.py`

Findings addressed:

- `.dockerignore` used unanchored `data`, `runs`, `results`, and `work`
  patterns. These are now root-anchored so `src/microsuite/data` remains in the
  Docker build context.
- Container docs implied Dockerfiles were only checked statically. They now
  state that CI builds the lighter `microsuite` and `kraken2` images by default,
  while heavy QIIME 2 and R images are manual-gated.
- CI tests did not lock heavy-image gating. They now assert two light and two
  heavy image matrix entries and smoke-test command presence.
- CI built images but did not smoke-test them. The default Docker images now run
  `microsuite --help` and `kraken2 --version` after build.

Verification:

```bash
docker build -f containers/microsuite/Dockerfile -t microsuite/microsuite:ci .
docker build -f containers/kraken2/Dockerfile -t microsuite/kraken2:ci .
docker run --rm microsuite/microsuite:ci --help
docker run --rm --entrypoint kraken2 microsuite/kraken2:ci --version
uv run --extra dev pytest
uv run --extra dev ruff check .
uv run --extra dev ruff format --check .
uv run --extra dev ty check
uv build
```

Residual risk:

- Heavy `qiime2-amplicon` and `r-diffab` images are not built by default because
  they are large and network-sensitive. They remain available through manual CI
  dispatch with `build-heavy-containers=true`.
