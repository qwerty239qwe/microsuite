# External Integration Tests

Default CI tests wrapper validation without requiring large third-party
toolchains or databases. Real external-tool integration tests are opt-in:

```bash
MICROSUITE_RUN_EXTERNAL_INTEGRATION=1 uv run pytest tests/integration
```

These tests call the Python API and execute real tools on `PATH`. Tests skip
individually when a tool is not installed.

GitHub Actions also runs a lightweight external integration job. It installs
`fastp`, `vsearch`, `mafft`, and `fasttree` with apt, installs `cutadapt` and
`multiqc` into the uv environment, and runs this same test module with
`MICROSUITE_RUN_EXTERNAL_INTEGRATION=1`.

No-database tests currently cover:

- `fastp`
- `cutadapt`
- `multiqc`
- `vsearch`
- `mafft` plus `FastTree`

Database-dependent tests require extra environment variables:

```bash
MICROSUITE_RUN_EXTERNAL_INTEGRATION=1 \
MICROSUITE_KRAKEN2_DB=/path/to/kraken2-db \
uv run pytest tests/integration/test_external_tools.py::test_kraken2_python_api_real_database
```

```bash
MICROSUITE_RUN_EXTERNAL_INTEGRATION=1 \
MICROSUITE_BRACKEN_DB=/path/to/bracken-db \
MICROSUITE_BRACKEN_REPORT=/path/to/kraken-report.tsv \
uv run pytest tests/integration/test_external_tools.py::test_bracken_python_api_real_database
```

```bash
MICROSUITE_RUN_EXTERNAL_INTEGRATION=1 \
MICROSUITE_METAPHLAN_DB=/path/to/metaphlan-db \
uv run pytest tests/integration/test_external_tools.py::test_metaphlan_python_api_real_database
```

Optional Bracken overrides:

- `MICROSUITE_BRACKEN_LEVEL`, default `S`
- `MICROSUITE_BRACKEN_READ_LENGTH`, default `150`
