# MicroSuite User-Experience Hardening — Multi-Agent Implementation Plan

> **For agentic workers:** implement only the assigned work package. Do not edit
> another package's reserved files. Every package must finish with tests, a concise
> handoff note, and one or more focused commits. The integration agent owns shared
> CLI registration, capability aggregation, release metadata, and cross-package
> validation.

**Goal:** eliminate the recurring late failures, version mismatches, FASTQ-layout
guessing, classifier download problems, incomplete run metadata, and invalid
repeated-measures inference observed while running the oral microbiome workflows.

**Target release:** `0.2.0`

**Architecture:** MicroSuite owns reusable method execution, validated data and
artifact contracts, runtime diagnostics, and stage metadata. Project repositories
own biological study decisions, metadata recoding, project DAG selection, and
presentation-specific reporting. New contracts are JSON-first and versioned; TSV
is retained as a human-readable projection where useful.

**Tech stack:** Python 3.11/3.12, Typer, dataclasses, JSON Schema, pandas/numpy,
Docker or Podman, R/vegan for formula-based ecological statistics, pytest,
ruff, ty.

---

## 1. Why This Work Exists

Recent changes already fixed or added:

- DADA2 output ownership, diagnostic plots, sample-name preservation, retention
  warnings, overlap checks, parameter sweeps, and 454 alignment controls.
- TSV table normalization/export and metadata-aware visualizations.
- Separate differential-abundance images and ANCOM-BC 2.4 compatibility.
- `stage-result.v1`, `resolved-config.v1`, metrics with units, and initial stage
  wrappers for denoising, taxonomy, and diversity.
- Native one-factor PERMDISP.

Do not reimplement those features. This plan addresses the remaining gaps:

1. A workflow can request a capability such as PERMDISP while an analysis PC
   silently runs an older checkout.
2. Expected `MicrobiomeSuiteError` failures currently produce full tracebacks.
3. Runtime, output permission, image, memory, and QIIME compatibility problems
   are discovered only after expensive stages begin.
4. Mixed FASTQ naming, incomplete pairs, technical splits, and multi-region
   accessions require project-specific Python.
5. Pretrained SILVA classifier acquisition and validation still happens outside
   the core reference-database service.
6. Most commands do not emit stage-result envelopes, so Microboard still needs
   filename guessing.
7. Native beta significance is one-factor with unrestricted permutations and is
   unsuitable for paired or longitudinal confirmatory inference.
8. Differential-abundance design failures are reported too late and without a
   reusable machine-readable design report.

---

## 2. Scope Boundaries

### MicroSuite owns

- Stable version/capability reporting and environment diagnostics.
- FASTQ discovery, validation, technical-run staging, and primer evidence.
- Reference artifact acquisition, caching, integrity checks, and compatibility
  preflight.
- Container invocation, UID/GID handling, image provenance, and runtime checks.
- Versioned metadata schemas and writers.
- Stage-level resolved parameters, artifacts, versions, metrics, and errors.
- Reusable statistical engines and design validation.

### Project workflow repositories continue to own

- Which samples, marker regions, primers, and classifiers are scientifically
  appropriate for a study.
- Study-specific metadata recoding and cohort exclusions.
- The end-to-end branch DAG and stage selection.
- FOMC algorithm policy until its interface is stabilized as a standalone method.
- Project-specific plot composition and narrative reports.
- Cloud-sync behavior. MicroSuite may detect unreadable/placeheld files, but it
  must not attempt to control Dropbox or another sync client.

---

## 3. Multi-Agent Rules

### Branches and worktrees

Create one branch/worktree per package from the same integration base:

| Package | Suggested branch |
|---|---|
| UX foundation | `feat/ux-foundation` |
| Reads manifest | `feat/reads-manifest` |
| Classifier artifacts | `feat/classifier-artifacts` |
| Metadata contracts | `feat/workflow-metadata-contracts` |
| Ecological statistics | `feat/restricted-beta-tests` |
| Diffab validation | `feat/diffab-design-validation` |
| Integration/release | `feat/ux-hardening-integration` |

Agents must not pull or merge another feature branch into their branch. The
integration agent cherry-picks or merges completed packages in dependency order.

### Reserved shared files

Only the named owner may edit these until the integration wave:

| File or area | Owner |
|---|---|
| `src/microsuite/cli/app.py` | UX foundation, then integration |
| `src/microsuite/__init__.py` | UX foundation, then integration |
| `pyproject.toml` | UX foundation, then integration |
| `src/microsuite/runtime/container.py` | UX foundation |
| `src/microsuite/metadata/models.py` | Metadata contracts |
| `src/microsuite/metadata/schemas.py` | Metadata contracts |
| `src/microsuite/metadata/validate.py` | Metadata contracts |
| `src/microsuite/metadata/_schema/` | Metadata contracts |
| `docs/methods.md` | Integration |
| `README.md` | Integration |
| capability registry entries | Integration after package handoffs |

Feature agents may create CLI modules but must not register them in `app.py`.
Their handoff must list the import and `app.add_typer(...)` line required.

### Required handoff from every agent

Create a short handoff in the PR/commit message or agent response containing:

1. Commits to integrate.
2. Interfaces added or changed.
3. Shared-file registrations still required.
4. Tests run and results.
5. Known limitations or deferred work.
6. Capability identifiers to add to the registry.

---

## 4. Locked Cross-Package Contracts

These contracts are decided for this implementation. Agents must not invent
incompatible alternatives without stopping for architecture review.

### 4.1 CLI error behavior

- A known `MicrobiomeSuiteError` prints exactly one `Error: ...` message and exits
  with code 1 without a traceback.
- Typer usage errors retain Typer's standard message and exit code 2.
- Unexpected exceptions retain a traceback.
- Global `--debug` causes known errors to include their traceback.

### 4.2 Version and capabilities

Commands:

```text
microsuite version [--json]
microsuite capabilities [--json]
microsuite doctor [--json] [--require CAPABILITY ...]
```

`version --json` returns at least:

```json
{
  "name": "microsuite",
  "version": "0.2.0",
  "source": "installed|editable",
  "commit": "<sha-or-null>",
  "python": "3.x.y"
}
```

`capabilities --json` uses schema `microsuite-capabilities.v1`:

```json
{
  "schema_version": "microsuite-capabilities.v1",
  "producer": {"name": "microsuite", "version": "0.2.0"},
  "capabilities": {
    "diversity.beta_significance.permdisp.native": {"available": true, "api": 1},
    "reads.manifest": {"available": true, "api": 1}
  }
}
```

Capability keys are stable API identifiers, not Git commit hashes. Project
wrappers must test capability keys rather than parsing help output or requiring a
specific commit.

### 4.3 Doctor result

Schema `microsuite-doctor.v1` contains checks with:

```text
id, status(pass|warn|fail|skip), message, details, remediation
```

Human output is a compact table. JSON output is stable and contains no ANSI.
Exit code is 0 when there are no failures, 1 when any required check fails, and 2
only for CLI usage errors.

### 4.4 Reads manifest

Canonical format is `reads-manifest.v1.json`. It contains one biological sample
record with one or more technical read sets:

```json
{
  "schema_version": "reads-manifest.v1",
  "samples": [{
    "sample_id": "S01",
    "layout": "PE",
    "read_sets": [
      {"run_id": "lane1", "read1": "/abs/S01_L1_R1.fastq.gz", "read2": "/abs/S01_L1_R2.fastq.gz"},
      {"run_id": "lane2", "read1": "/abs/S01_L2_R1.fastq.gz", "read2": "/abs/S01_L2_R2.fastq.gz"}
    ],
    "metadata": {}
  }],
  "warnings": []
}
```

Rules:

- Paths are absolute in the canonical JSON.
- One sample cannot mix SE and PE read sets.
- R2 without R1 is an error by default and a warning only under an explicit
  permissive flag.
- Technical sets are never silently discarded.
- Staging concatenates corresponding gzip streams atomically into deterministic
  sample-level files and writes `source_manifest.tsv`.
- TSV sample sheets remain accepted as an input adapter, but are converted to the
  canonical JSON model before validation or execution.

### 4.5 Classifier artifact identity

Classifier cache identity includes:

```text
database, database_release, training_framework, marker_region, primer_set,
distribution(uniform|weighted), habitat, provider, artifact_format
```

The downloaded bytes are not installed into the cache until archive validation
and checksum validation succeed. Cache records include SHA-256 and source URL.

### 4.6 Metadata artifacts

MicroSuite publishes and validates:

```text
stage-result.v1
resolved-config.v1
workflow.v1
run-manifest.v1
reads-manifest.v1
microsuite-capabilities.v1
microsuite-doctor.v1
```

Project runners may write workflow/run instances through MicroSuite utilities,
but MicroSuite must not infer a workflow DAG by scanning result folders.

### 4.7 Confirmatory beta significance

- Existing native one-factor behavior remains backward-compatible.
- Formula and restricted permutations use a new `vegan` backend.
- CLI accepts `--formula`, optional `--strata`, `--permutations`, and `--seed`.
- Output always records backend, formula, strata, permutation scheme, sample
  count, group/design degrees of freedom, statistic, and p-value.
- PERMANOVA and PERMDISP are separate methods but share aligned samples and the
  same declared permutation design.

---

## 5. Dependency Waves

```text
Wave 0: UX foundation
    |
    +-- Wave 1A: reads manifest
    +-- Wave 1B: classifier artifacts
    +-- Wave 1C: workflow metadata contracts
    +-- Wave 1D: restricted beta tests
    +-- Wave 1E: diffab design validation
              |
              +-- Wave 2: integration, capability registry, docs, release
              |
              +-- Wave 3: oral workflow migration (separate repository)
```

Wave 1 packages may run concurrently after Wave 0 contracts are merged. Metadata
stage coverage within a feature package should use the existing `stage_execution`
API; additions to the core metadata models wait for Wave 1C.

---

## 6. Work Package UX: CLI Foundation, Version, Capabilities, Doctor

**Owner:** UX foundation agent

**Dependencies:** none

**Reserved files:** `cli/app.py`, `__init__.py`, `pyproject.toml`,
`runtime/container.py`

**Create:**

- `src/microsuite/cli/system_cmd.py`
- `src/microsuite/system/__init__.py`
- `src/microsuite/system/version.py`
- `src/microsuite/system/capabilities.py`
- `src/microsuite/system/doctor.py`
- `tests/test_cli_error_handling.py`
- `tests/test_system_version.py`
- `tests/test_system_capabilities.py`
- `tests/test_system_doctor.py`

### Tasks

- [ ] Add failing tests proving a known `MicrobiomeSuiteError` exits 1 without
  `Traceback` or `click.exceptions.Exit` in output.
- [ ] Implement clean error handling and a global `--debug` switch.
- [ ] Replace the duplicated hard-coded `0.1.0` version with one authoritative
  version source. Prefer `importlib.metadata.version("microsuite")`, with a safe
  source-tree fallback for tests.
- [ ] Implement `version` human and JSON output. For editable source checkouts,
  include a best-effort Git commit without making Git a runtime requirement.
- [ ] Implement a static capability registry with schema versioning. Initial
  registry describes capabilities present on the Wave 0 base only.
- [ ] Implement `doctor` checks for Python version, output/cache writability,
  Docker/Podman presence, Docker socket access, host UID/GID, and optional
  executable discovery.
- [ ] Add reusable container probes to `runtime/container.py`: engine info,
  image existence, image inspect, and a non-destructive bind-mount write test.
- [ ] Ensure every probe has a timeout and returns structured data instead of
  printing internally.
- [ ] Add `--require` capability validation with actionable remediation.
- [ ] Register `version`, `capabilities`, and `doctor` in `cli/app.py`.

### Acceptance tests

```bash
uv run pytest tests/test_cli_error_handling.py tests/test_system_*.py -v
uv run microsuite version --json
uv run microsuite capabilities --json
uv run microsuite doctor --json
```

Tests mock Docker and Git. They must not require a running daemon.

### Capability handoff

```text
system.version api=1
system.capabilities api=1
system.doctor api=1
```

---

## 7. Work Package READS: Manifest, Validation, Staging, Primer Evidence

**Owner:** reads agent

**Dependencies:** Wave 0

**Do not edit:** `cli/app.py`, metadata schema registry, runtime container files

**Create:**

- `src/microsuite/reads/__init__.py`
- `src/microsuite/reads/models.py`
- `src/microsuite/reads/discovery.py`
- `src/microsuite/reads/validate.py`
- `src/microsuite/reads/stage.py`
- `src/microsuite/reads/primers.py`
- `src/microsuite/cli/reads_cmd.py`
- `tests/test_reads_discovery.py`
- `tests/test_reads_manifest.py`
- `tests/test_reads_stage.py`
- `tests/test_reads_primers.py`
- `tests/fixtures/reads/`

### CLI contract

```text
microsuite reads manifest INPUT_DIR --output manifest.json
  [--sample-sheet samples.tsv] [--recursive] [--allow-orphans]

microsuite reads validate manifest.json [--json]

microsuite reads stage manifest.json --output-dir prepared/
  [--method concat|copy|symlink] [--force]

microsuite reads inspect-primers manifest.json --primer NAME=SEQUENCE ...
  --output primer-evidence.tsv [--reads-per-set 1000]
```

### Tasks

- [ ] Implement frozen dataclasses and dependency-free JSON serialization for
  `reads-manifest.v1`.
- [ ] Detect `.fastq`, `.fq`, and gzip forms using established suffix families:
  `R1/R2`, `_1/_2`, `read1/read2`, and `forward/reverse`, including optional lane
  and `_001` components.
- [ ] Keep unmatched singleton files as SE. Never infer PE from alternating file
  order.
- [ ] Report duplicate mates, R2-only samples, mixed layout, empty files, invalid
  gzip streams, and duplicate sample IDs with targeted messages.
- [ ] Import TSV sample sheets with `sample_id` or `Run`, `read1`, optional
  `read2`, and optional `source_runs`. Resolve relative paths relative to the
  sample sheet, not process CWD.
- [ ] Stage multiple technical read sets atomically. `concat` copies gzip members
  byte-for-byte; do not decompress/recompress. Validate the resulting stream.
- [ ] Preserve source paths, sizes, and run IDs in `source_manifest.tsv` and a
  JSON provenance file.
- [ ] Make staging idempotent by comparing a source fingerprint manifest; a same-
  sized but changed source must not be treated as current.
- [ ] Implement IUPAC-aware primer matching on bounded prefixes/suffixes. Report
  evidence and orientation; do not automatically drop or reclassify samples.
- [ ] Add fixtures based on the observed classes: normal PE, singleton SE,
  incomplete mate, mixed naming, two technical lanes, and two marker-primer
  signatures. Use synthetic tiny FASTQs, not project data.
- [ ] Create `reads_cmd.app`, but leave registration to the integration agent.

### Acceptance tests

```bash
uv run pytest tests/test_reads_*.py -v
uv run microsuite reads manifest tests/fixtures/reads/mixed --output /tmp/reads.json
uv run microsuite reads validate /tmp/reads.json --json
```

### Capability handoff

```text
reads.manifest api=1
reads.validate api=1
reads.stage api=1
reads.primer_evidence api=1
```

---

## 8. Work Package REFDB: Pretrained Classifier Artifacts

**Owner:** classifier/refDB agent

**Dependencies:** Wave 0

**May modify:** `refdb/`, `qiime2/artifact.py`, `cli/refdb_cmd.py`,
`methods/tax_classify.py`

**Do not edit:** `runtime/container.py`, `cli/app.py`, core metadata models

**Create or extend:**

- `src/microsuite/refdb/classifiers.py`
- `src/microsuite/refdb/providers/biodbs.py`
- `src/microsuite/refdb/registry.py`
- `src/microsuite/qiime2/compatibility.py`
- `tests/test_refdb_classifiers.py`
- `tests/test_qiime2_compatibility.py`
- `tests/fixtures/qza/`

### CLI contract

```text
microsuite refdb classifier list silva --release 138.2 [--json]

microsuite refdb classifier fetch silva \
  --release 138.2 --training-framework 2025.7 \
  --region V4V5-515f-926r --distribution weighted \
  --habitat human-oral [--locator PATH] [--force]

microsuite refdb classifier inspect CLASSIFIER.qza [--json]

microsuite refdb classifier validate CLASSIFIER.qza \
  [--runtime local|docker] [--image IMAGE] [--engine docker|podman]
```

### Tasks

- [ ] Add `ClassifierSpec` and `ClassifierArtifact` without weakening existing
  raw-reference `RefDbSpec` behavior.
- [ ] Extend the registry so classifier cache keys include every identity field
  in section 4.5 and records include source URL, byte size, SHA-256, QIIME type,
  framework version, and archive version.
- [ ] Implement provider listing separately from direct download. An empty listing
  is not proof that a known locator is absent.
- [ ] Download to a unique temporary file in the target filesystem, reject HTTP
  error/HTML bodies, inspect the QIIME ZIP, verify checksum when supplied, then
  atomically install.
- [ ] Preserve an invalid existing artifact as `.invalid.<UTC timestamp>` only
  under explicit repair/fetch behavior; validation alone is read-only.
- [ ] Improve `inspect_artifact` error normalization to catch `BadZipFile`,
  malformed metadata, and missing payloads as concise `MicrobiomeSuiteError`s.
- [ ] Implement authoritative runtime compatibility validation by running
  `qiime tools validate` in the configured local environment or container. Use
  Wave 0 container probes and UID/GID mapping.
- [ ] Make `tax_classify --backend qiime2` run cheap local archive inspection and
  optional runtime compatibility preflight before classification.
- [ ] Emit reference identity/checksum through the existing taxonomy stage record.
- [ ] Test invalid HTML, truncated ZIP, valid minimal QZA, unsupported archive,
  checksum mismatch, atomic download failure, cache hit, and force repair.
- [ ] Do not perform real network access in unit tests. Put live provider probes
  behind a separately marked integration test.

### Acceptance tests

```bash
uv run pytest tests/test_refdb_classifiers.py tests/test_qiime2_compatibility.py \
  tests/test_refdb_tax_integration.py -v
```

### Capability handoff

```text
refdb.classifier_artifacts api=1
qiime2.classifier_preflight api=1
taxonomy.reference_provenance api=1
```

---

## 9. Work Package META: Workflow/Run Schemas And Aggregation Utilities

**Owner:** metadata-contract agent

**Dependencies:** Wave 0

**Reserved files:** all core metadata model/schema/validator files

**Create or modify:**

- `src/microsuite/metadata/models.py`
- `src/microsuite/metadata/schemas.py`
- `src/microsuite/metadata/validate.py`
- `src/microsuite/metadata/workflow.py`
- `src/microsuite/metadata/manifest.py`
- `src/microsuite/metadata/_schema/workflow.v1.schema.json`
- `src/microsuite/metadata/_schema/run-manifest.v1.schema.json`
- `tests/test_workflow_metadata.py`
- `tests/test_run_manifest_metadata.py`
- `tests/fixtures/workflow/`
- `tests/fixtures/run_manifest/`

### Public Python API

```python
write_workflow(run_dir, workflow)
write_run_manifest(run_dir, manifest)
aggregate_stage_results(run_dir, *, expected_stages)
validate_workflow(payload)
validate_run_manifest(payload)
workflow_hash(payload)
```

### Tasks

- [ ] Translate the approved Microboard integration requirements into published
  `workflow.v1` and `run-manifest.v1` JSON Schemas.
- [ ] Represent DAG nodes, dependencies, branches, method selections, and
  intermediate-file definitions explicitly. Do not infer them from directories.
- [ ] Store immutable workflow version and SHA-256 over canonical JSON excluding
  only the hash field itself.
- [ ] Represent dataset/run identity, progress, stage attempts, provenance,
  artifacts, and summary metrics in `run-manifest.v1`.
- [ ] Aggregate `stage-results/*.json` by declared stage and attempt. Ignore
  diagnostics/non-JSON files. Preserve failed attempts and identify the latest
  attempt without deleting history.
- [ ] Define progress deterministically from expected DAG nodes and latest attempt
  statuses.
- [ ] Reuse atomic writes and secret redaction. Do not redact biological sample
  IDs unless explicitly configured by the caller.
- [ ] Publish schema fixtures and parity-test them with both the dependency-free
  validator and `jsonschema` when installed.
- [ ] Add APIs to `metadata/__init__.py`.
- [ ] Keep existing stage and resolved-config schemas backward-compatible.
- [ ] Document how Microboard presentation overlays remain separate from immutable
  execution records.

### Acceptance tests

```bash
uv run pytest tests/test_metadata_*.py tests/test_workflow_metadata.py \
  tests/test_run_manifest_metadata.py -v
```

### Capability handoff

```text
metadata.workflow.v1 api=1
metadata.run_manifest.v1 api=1
metadata.stage_aggregation api=1
```

---

## 10. Work Package BETA: Formula-Based Restricted PERMANOVA/PERMDISP

**Owner:** ecological-statistics agent

**Dependencies:** Wave 0

**May modify:** `diversity/`, `cli/diversity_cmd.py`, new R resources and
ecology container files

**Do not edit:** `cli/app.py`, core metadata schemas, shared diffab R files

**Create or modify:**

- `src/microsuite/diversity/vegan.py`
- `src/microsuite/diversity/r/beta_significance.R`
- `src/microsuite/cli/diversity_cmd.py`
- `containers/r-ecology/Dockerfile`
- `.github/workflows/` only if an isolated ecology smoke job can be added without
  editing another active agent's workflow file
- `tests/test_beta_significance_vegan.py`
- `tests/fixtures/diversity/`

### CLI contract

Existing calls remain valid:

```text
microsuite diversity beta-significance distances.tsv \
  --metadata metadata.tsv --column group --method permanova
```

New confirmatory calls:

```text
microsuite diversity beta-significance distances.tsv \
  --metadata metadata.tsv --backend vegan \
  --formula "site_group + phase_code" --strata subject_code \
  --method permanova --runtime docker --image IMAGE
```

### Tasks

- [ ] Preserve the native backend output and defaults exactly.
- [ ] Add `backend=native|vegan`; require `vegan` for formulas, covariates, or
  strata.
- [ ] Use R `vegan::adonis2` for PERMANOVA and `betadisper` plus restricted
  `permutest` for PERMDISP. Do not hand-roll multifactor sums of squares.
- [ ] Align distance matrix and metadata once, reject duplicate sample IDs, and
  report dropped metadata-only samples.
- [ ] Validate formula columns, missing values, group counts, residual degrees of
  freedom, strata sizes, and permutation feasibility before the expensive test.
- [ ] Make permutation design deterministic from seed. Record unrestricted,
  blocked, or no-permutation status explicitly.
- [ ] Return tidy TSV with one row per model term for PERMANOVA. PERMDISP returns
  the tested grouping factor and must reject a formula with no unambiguous factor
  unless `--column` is supplied.
- [ ] Add local and Docker execution through a dedicated `r-ecology` image. Run as
  host UID/GID and record image digest sidecar.
- [ ] Wrap vegan execution in `stage_execution` and emit formula, strata,
  permutation count, statistic, p-value, and sample count metrics.
- [ ] Test paired synthetic data where blocked and unrestricted p-values differ,
  singleton strata warnings, confounded formulas, missing metadata, deterministic
  seeds, and container command construction.
- [ ] Add a real-input build smoke fixture small enough for CI.

### Acceptance tests

```bash
uv run pytest tests/test_ecological_statistics.py \
  tests/test_beta_significance_vegan.py -v
docker build --progress=plain -t microsuite-r-ecology-test containers/r-ecology
```

### Capability handoff

```text
diversity.beta_significance.formula.vegan api=1
diversity.beta_significance.strata.vegan api=1
diversity.permdisp.vegan api=1
```

---

## 11. Work Package DIFFAB: Design Validation And Stable Results Contract

**Owner:** differential-abundance agent

**Dependencies:** Wave 0

**May modify:** `diffab/`, `cli/diffab_cmd.py`, backend-specific R scripts/tests

**Do not edit:** core metadata files, `runtime/container.py`, `cli/app.py`

**Create or modify:**

- `src/microsuite/diffab/design.py`
- `src/microsuite/cli/diffab_cmd.py`
- `src/microsuite/diffab/r/validate_design.R`
- `tests/test_diffab_design.py`
- `tests/test_diffab_result_contract.py`

### CLI contract

```text
microsuite diffab validate-design TABLE.h5ad \
  --fix-formula "site_group + phase_code" \
  [--rand-formula "(1|subject_code)"] [--reference col=level ...] \
  --output design-report.json
```

### Tasks

- [ ] Implement a read-only design preflight shared by ANCOM-BC2 and MaAsLin2.
- [ ] Validate referenced columns and levels, all-missing/constant columns,
  duplicate sample IDs, non-finite counts, zero library samples, and random-effect
  grouping cardinality.
- [ ] Build the definitive R model matrix using `model.matrix`; report numerical
  rank, column count, aliased columns, and likely confounded source variables.
- [ ] Produce `diffab-design-report.v1` JSON with `valid`, errors, warnings,
  formula, random formula, references, dimensions, and factor levels.
- [ ] Make ANCOM-BC2 and MaAsLin2 call the same preflight before model execution.
  A failed design must not create a partial result TSV.
- [ ] Define a normalized result schema across supported ANCOM-BC versions:
  feature, term/contrast, estimate, standard error, statistic, p-value, q-value,
  significance, method, plus backend-specific columns.
- [ ] Preserve raw backend output alongside normalized output rather than dropping
  version-specific fields.
- [ ] Ensure missing expected q-value families fail with a targeted contract error
  that lists actual columns.
- [ ] Wrap preflight and model execution in one differential-abundance stage
  envelope, including container image/digest and design report provenance.
- [ ] Test rank deficiency, invalid references, repeated measures, one-level
  factors, ANCOM-BC 2.4 column families, and older fixture column families.

### Acceptance tests

```bash
uv run pytest tests/test_diffab_design.py tests/test_diffab_result_contract.py \
  tests/test_ancombc*.py tests/test_diffab_runner.py -v
```

### Capability handoff

```text
diffab.design_validation api=1
diffab.normalized_results api=1
```

---

## 12. Work Package INTEGRATION: Registration, Metadata Coverage, Release

**Owner:** integration agent

**Dependencies:** all accepted Wave 1 packages

This agent resolves shared-file changes and is the only agent that performs the
final release-wide refactor.

### Tasks

- [ ] Merge/cherry-pick Wave 0 first, then each Wave 1 package independently.
  Run that package's focused tests after each integration.
- [ ] Register `reads_cmd.app` and any new nested refDB commands in `cli/app.py`.
- [ ] Add all handed-off capability identifiers to the capability registry.
- [ ] Add schemas for `reads-manifest.v1`, capability, doctor, and diffab design
  to the central schema registry if their package kept local validation during
  parallel development.
- [ ] Resolve stage-result wrapping consistently across:
  - reads staging
  - trim and QC
  - import/export and normalization
  - abundance and assignment QC
  - refDB acquisition
  - differential abundance
  - visualization
- [ ] For every wrapped stage, declare effective parameters after defaults and
  overrides, typed inputs/outputs, software/container versions, metrics with
  units, and required provenance files.
- [ ] Do not write duplicate envelopes when a high-level stage calls a lower-level
  helper. The public biological stage is the envelope boundary.
- [ ] Add `microsuite doctor --classifier ... --image ...` composition so project
  runners can preflight QIIME archive compatibility, image availability, memory,
  disk, and output writability in one command.
- [ ] Update `docs/methods.md`, `docs/containers.md`, `docs/api-python.md`, and
  `README.md` with final CLI examples and capability identifiers.
- [ ] Update the historical complaint document or add a release note marking
  resolved items. Do not leave obsolete items described as open.
- [ ] Add `CHANGELOG.md` entry for 0.2.0 and finalize package version handling.
- [ ] Add a migration guide from help/commit probing and project-specific FASTQ
  detection to capabilities and reads manifests.

### Release-wide acceptance gates

```bash
uv run ruff check .
uv run ruff format --check .
uv run ty check
uv run pytest -q
python -m build
```

Container gates:

```bash
docker build --progress=plain -t microsuite-r-dada2-test containers/r-dada2
docker build --progress=plain -t microsuite-r-ecology-test containers/r-ecology
docker build --progress=plain -t microsuite-r-diffab-ancombc-test containers/r-diffab-ancombc
```

CLI smoke gates:

```bash
uv run microsuite version --json
uv run microsuite capabilities --json
uv run microsuite doctor --json
uv run microsuite reads --help
uv run microsuite refdb classifier --help
uv run microsuite diversity beta-significance --help
uv run microsuite diffab validate-design --help
```

Backward-compatibility gates:

- Existing native diversity commands produce the same columns and defaults.
- Existing `refdb fetch` remains valid.
- Existing direct FASTQ `trim`, `qc`, and `denoise` calls remain valid.
- Existing `stage-result.v1` and `resolved-config.v1` fixtures remain valid.
- Existing local and Docker differential-abundance calls remain valid.

---

## 13. Separate Repository Migration After 0.2.0

This work occurs in `oral_microbiome_gingivitis` only after the integrated
MicroSuite release is available.

### Tasks

- [ ] Add a single preflight at the start of `run_amplicon_end_to_end.sh` using
  `microsuite doctor --require ...` and JSON capability checks.
- [ ] Replace project FASTQ detection and SRP009299 technical concatenation with
  `microsuite reads manifest/validate/stage` while keeping committed study sample
  selections and metadata mappings.
- [ ] Replace direct SILVA classifier download logic with
  `microsuite refdb classifier fetch/validate`.
- [ ] Replace script-owned stage metadata where MicroSuite now emits the required
  envelope; retain the project runner's workflow/run aggregation.
- [ ] Use restricted vegan PERMANOVA/PERMDISP for paired/longitudinal datasets.
  Keep native unrestricted outputs only as explicitly labelled exploratory tests.
- [ ] Run dry runs for every canonical project and one real small workflow smoke.
- [ ] Do not delete compatibility fallback code until one complete run imports
  correctly into Microboard using metadata files rather than folder scanning.

---

## 14. End-to-End Acceptance Scenarios

The integrated release is not complete until these synthetic or fixture-backed
scenarios pass.

### Scenario A: stale capability

A caller requires `diversity.beta_significance.permdisp.native` from a build that
does not advertise it. `doctor --require` fails before any analysis, names the
missing capability, prints the installed version/commit, and suggests upgrading.

### Scenario B: expected CLI error

`beta-significance --method invalid` exits 1 with one concise error and no
traceback. The same call with `--debug` includes a traceback.

### Scenario C: mixed reads

A directory contains `S1_R1/R2`, `S2_1/2`, singleton `S3`, two lanes for `S4`,
and orphan `S5_R2`. Strict manifest creation reports the orphan and creates no
silent partial pair. Permissive mode records a warning. Staging concatenates S4
lanes and preserves all source provenance.

### Scenario D: corrupt classifier

A `.qza` path contains HTML. Inspection identifies a non-QIIME archive before
QIIME starts. Fetch downloads atomically, validates ZIP/type/checksum, and installs
one cache artifact. A classifier created by an incompatible QIIME framework fails
the runtime validation probe before classification.

### Scenario E: cloud/output permissions

Doctor receives an unwritable output directory or a readable metadata file with
an unreadable/empty FASTQ target. It reports the exact path and remediation before
the workflow creates a run bundle.

### Scenario F: paired ecological test

A paired synthetic dataset is tested with `--strata subject`. Output identifies
blocked permutations and differs from the unrestricted exploratory result. The
same model's PERMDISP result uses the declared design and reports its method.

### Scenario G: rank-deficient differential abundance

Metadata contains perfectly confounded time and phase columns. `validate-design`
fails before model fitting, names the aliased model columns, and writes a valid
design report. No result TSV is created.

### Scenario H: Microboard contract

A small workflow produces stage results, immutable workflow snapshot/hash,
resolved configuration, and run manifest. Schema validation passes. Microboard can
identify stages, attempts, artifacts, metrics, methods, and failures without
searching folders.

---

## 15. Definition Of Done

- [ ] All package acceptance tests and release-wide gates pass.
- [ ] No known user-facing `MicrobiomeSuiteError` prints a traceback by default.
- [ ] A machine-readable capability check can prevent old-checkout failures.
- [ ] FASTQ manifests represent SE, PE, technical splits, and incomplete pairs
  without filename guessing in project scripts.
- [ ] Pretrained QIIME classifiers are acquired atomically and validated before
  classification.
- [ ] Every core analysis stage emits a valid stage-result envelope with resolved
  parameters and meaningful provenance.
- [ ] Workflow/run metadata can be consumed without directory scanning.
- [ ] Restricted formula-based PERMANOVA and PERMDISP are available for repeated
  designs.
- [ ] Differential-abundance design failures are reported before model fitting.
- [ ] Existing 0.1 command lines remain functional unless explicitly documented
  in the migration guide.
- [ ] The oral workflow repository has a follow-up migration issue/plan, but its
  project-specific biological decisions remain outside MicroSuite.
