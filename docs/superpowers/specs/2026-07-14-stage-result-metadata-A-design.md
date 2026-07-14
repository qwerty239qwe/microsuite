# Structured run metadata for Microboard — Sub-project A (schema + envelope foundation) — Design

- **Date:** 2026-07-14
- **Status:** Approved (design, with adjustments), pending implementation plan
- **Origin:** Microboard integration. The neighbouring `microboard` dashboard
  currently imports `microsuite-results.json` (schema `microsuite-results.v1`)
  and otherwise *guesses* stage metrics by scanning result folders
  (`dada2_qc_summary.json`, `end_to_end_runs/*/manifest.tsv`). We want microsuite
  to emit a first-class, versioned, per-stage metadata envelope so Microboard
  stops filename-guessing and can reliably compare benchmark runs (e.g. DADA2 vs
  Deblur). Agreed decomposition: **A** (this — schema registry, validator,
  `StageResult` writer, `run_command` integration), **B** (software/refDB
  versions + `resolved_config.json`), **C** (metrics with units).

## Scope of A

Deliver the metadata *foundation*: a versioned schema registry, a dependency-free
validator, typed `StageResult` / `Artifact` / `ProvenanceFile` models, an atomic
per-stage-attempt writer, stable workflow/stage identifiers, and `run_command`
integration that emits a `stage-result.v1.json` on **success, failure, and
timeout** — from explicit artifact declarations, never from folder/argv guessing.
The `metrics`, `software`, and `reference_db` slots ship **empty** in A; B and C
populate them.

### Out of scope for A
- Populating `metrics` (C), `software` / `reference_db` / `tool_versions` (B),
  and writing `resolved_config.json` (B).
- Any change to `microsuite-results.json`. It **remains backward-compatible and
  unchanged during A** (B/C will later add `tool_versions` / `metrics` to it).
- The workflow topology (nodes, edges, branches, selected methods). That belongs
  in a separate `workflow-run.v1.json` (a later sub-project). A `stage-result`
  describes one *execution*, never the whole graph.
- Wiring every microsuite stage to declare typed artifacts. A ships the
  mechanism + generic `run_command` integration + **one worked reference stage
  (`denoise`)**; converting the remaining stages is incremental follow-up.
- Microboard-side ingestion changes (documented here for context only).

## Verified current state

- `runtime/results.py` writes/append-updates `microsuite-results.json`
  (`microsuite-results.v1`): `executions[]` (`task, backend, command, started_at,
  inputs, outputs, params, duration_sec, exit_code`) + derived `artifacts[]`.
  Written by `run_command` **only on exit 0**, only when `run_dir` is set.
- `runtime/runner.py` `run_command` already records `command.txt`, `events.jsonl`,
  `run.json`, `stdout/stderr.log` per `run_dir`; `CommandLog` carries
  `task, backend, inputs, outputs, params` (all `Any`/dict today).
- No schema registry, no validator, no per-stage envelope, no secret redaction,
  no failure-path metadata.
- No `jsonschema`/`pydantic` dependency (deps are stdlib + scientific stack) →
  the validator is dependency-free by necessity and design.
- Microboard (`microboard/src/domain.js`) ingestion the envelope must eventually
  slot into: `microsuite-results.v1` → sets `metrics: {}`; metrics only from
  folder scan. Target ingestion priority (Microboard side, not A): supported
  `stage-result.*.json` → legacy `microsuite-results.json` → folder discovery.

## Design

### File layout & lifecycle (adjustment #1)

Never overwrite. One file **per stage attempt**, under a subdirectory of the
run directory:

```
run_dir/
  stage-results/
    denoise--dada2-r--attempt-1.json
    denoise--dada2-r--attempt-2.json      # a retry
    taxonomy--silva--attempt-1.json
  microsuite-results.json                 # legacy aggregate, unchanged
```

- Filename: `<stage>--<backend>--attempt-<N>.json`, each component slugified
  (`[^a-z0-9]+` → `-`, lower-cased; `backend` omitted → `none`).
- `attempt N` = `1 + count(existing stage-results/<stage>--<backend>--attempt-*.json)`.
  This is the retry/benchmark counter and is also written into the payload.
- If `run_dir` is `None`, no envelope is written (mirrors current behaviour).

### Identifiers (adjustments #1, #2)

Every envelope carries execution identity; workflow identity is optional and
injected by the orchestrator via environment (no new CLI flags):

| Field | Source | Notes |
|---|---|---|
| `run_id` | env `MICROSUITE_RUN_ID`, else `run_dir.name` | stable per workflow execution |
| `stage_run_id` | generated `stage-run-<uuid4hex12>` | unique per stage invocation/attempt |
| `attempt` | computed from existing files | int ≥ 1 |
| `workflow_id` | env `MICROSUITE_WORKFLOW_ID` | optional; null when absent |
| `workflow_run_id` | env `MICROSUITE_WORKFLOW_RUN_ID` | optional; distinguishes two DADA2 runs in different benchmarks |
| `dataset_id` | env `MICROSUITE_DATASET_ID` | optional |

A `metadata/context.py` `workflow_context()` reads these once (with optional
explicit overrides for tests).

### Models (`metadata/models.py`) (adjustment #3)

Frozen dataclasses; the writer serialises them (omitting `None`/absent optionals
except where the schema requires an explicit `null`).

```python
@dataclass(frozen=True)
class ArtifactCount:
    value: int
    unit: str            # what is counted: "samples", "features", "reads", ...

@dataclass(frozen=True)
class Artifact:
    label: str
    path: str            # relative to run_dir when under it; else absolute
    format: str | None = None      # "tsv", "fastq.gz", "directory", ...
    kind: str | None = None        # "feature_table", "representative_sequences", ...
    count: ArtifactCount | None = None   # declared by the stage; writer does NOT count
    bytes: int | None = None       # filled by writer via stat(); None if dir/missing
    external: bool = False         # true for inputs outside run_dir
    exists: bool | None = None     # filled by writer for outputs

@dataclass(frozen=True)
class ProvenanceFile:
    kind: str            # "dada2_manifest", "ancombc_provenance", "container", ...
    path: str            # relative to run_dir
    exists: bool | None = None      # filled by writer
```

`CommandLog` (in `runtime/runner.py`) gains **new optional** typed fields
alongside the untouched legacy `inputs/outputs/params` (which still feed
`microsuite-results.json`):

```python
stage: str | None = None                       # conceptual stage; defaults to task
artifact_inputs: tuple[Artifact, ...] = ()
artifact_outputs: tuple[Artifact, ...] = ()    # declared/expected outputs
provenance_files: tuple[ProvenanceFile, ...] = ()
```

The writer **enriches declared `artifact_outputs`** post-execution: resolve path,
set `exists`, set `bytes = stat().st_size` for regular files (`None` for
directories / missing). It never invents artifacts and never derives `count`
(counts are declared by the stage when known, else `None`).

### Envelope schema `stage-result.v1` (adjustments #2, #4, #5)

```json
{
  "schema_version": "stage-result.v1",
  "run_id": "results_erp120510_end_to_end",
  "stage_run_id": "stage-run-9f3c1a2b7d40",
  "workflow_id": "oral-standard",
  "workflow_run_id": "oral-standard--erp120510--001",
  "dataset_id": "ERP120510",
  "attempt": 1,
  "stage": "denoise",
  "task": "denoise",
  "backend": "dada2-r",
  "status": "completed",
  "exit_code": 0,
  "error": null,
  "timing": {
    "started_at": "2026-07-14T09:12:03Z",
    "finished_at": "2026-07-14T09:12:15Z",
    "duration_sec": 12.34
  },
  "command": ["Rscript", "denoise.R", "--token", "***"],
  "params": {"trunc_len_f": 240, "auth_token": "***"},
  "inputs":  [{"label": "reads", "path": "../input", "format": "directory",
               "external": true}],
  "outputs": [{"label": "feature table", "path": "feature-table.tsv",
               "format": "tsv", "kind": "feature_table",
               "count": {"value": 1842, "unit": "features"},
               "bytes": 928104, "exists": true}],
  "provenance_files": [{"kind": "dada2_manifest",
                        "path": "dada2_denoise_manifest.json", "exists": true}],
  "metrics": {},
  "software": {},
  "reference_db": null,
  "producer": {"name": "microsuite", "version": "0.x.y"}
}
```

Field rules:
- `status` enum: `running | completed | failed | timed_out | cancelled`. A emits
  `completed` (exit 0), `failed` (non-zero exit), `timed_out` (TimeoutExpired).
  `running`/`cancelled` are reserved for future emitters.
- `exit_code`: nullable — `null` for `timed_out`/`cancelled` (no meaningful code).
- `error`: `null` on success; else `{"type": str, "message": str}` with a **safe,
  redacted, truncated** message (e.g. `NonZeroExit` + last stderr line, run
  through the same redactor as `command`/`params`).
- Paths: **relative to `run_dir`** when the target is under it; otherwise
  absolute with `external: true` (inputs) so Microboard can tell them apart.
- `bytes`: nullable (missing outputs, directories).
- Timestamps: UTC RFC 3339 with `Z`.
- `command` and `params` are redacted by **one shared mechanism** (below).
- Unknown fields are permitted (forward-compat, #7); required/core fields are
  strictly validated.

### Failure & timeout semantics (adjustment #5)

`run_command` writes the envelope from a `finally`-style path that **preserves the
original exception**. On `failed`/`timed_out` the envelope still records: timing,
`params`, redacted `command`, declared `inputs`, declared `outputs` enriched with
whichever partial outputs exist (`exists`/`bytes`), `exit_code` when available,
structured `error`, and the `provenance_files` that were actually written
(`exists: true`). The envelope write must never mask or replace the pipeline's
real failure — any error inside the writer is swallowed (best-effort) so the
original `MicrobiomeSuiteError` propagates unchanged.

### Secret redaction (`metadata/redact.py`) (adjustment #4)

One mechanism used for both `params` and `command`:
- `SENSITIVE_KEY_RE = (?i)(token|secret|password|passwd|api[-_]?key|credential|auth)`.
- `redact_params(mapping)` → deep-copy; any key matching the regex → `"***"`
  (recurse into nested dicts/lists).
- `redact_command(argv)` → mask the value after a sensitive flag
  (`--token X` → `--token ***`), inline `--api-key=X` → `--api-key=***`, and any
  bare token equal to a value redacted from `params` (captured secret values).
- `redact_text(s, secrets)` → replace captured secret substrings in free text
  (used for `error.message`).

### Writer (`metadata/stage_result.py`)

`write_stage_result(run_dir, log, *, status, started_at, finished_at, exit_code,
error=None, command)`:
1. Build identity (`run_id`, `stage_run_id`, workflow ctx, `attempt` from disk).
2. Build the payload from `CommandLog` typed fields (relativise paths, enrich
   outputs with `exists`/`bytes`, enrich provenance existence), redact `params`
   and `command`, attach `producer` version, empty `metrics`/`software`, null
   `reference_db`.
3. `validate_stage_result(payload)` — on errors, log a warning and still write
   (never block the pipeline on a schema nit), so bugs surface without data loss.
4. **Atomic write**: `tmp = target.with_name(target.name + ".tmp.<pid>")`;
   `tmp.write_text(json.dumps(..., indent=2, sort_keys=True))`; `os.replace(tmp,
   target)`. Microboard never sees a partial file.

### Schema registry & validator (`metadata/schemas.py`, `metadata/validate.py`) (adjustments #6, #7)

A **MicroSuite-specific** schema format (a small subset inspired by JSON Schema —
*not* a JSON Schema validator; documented as such in the module docstring):

```python
STAGE_RESULT_V1 = {
    "required": ["schema_version", "run_id", "stage_run_id", "attempt",
                 "stage", "status", "producer"],
    "allow_unknown": True,          # forward-compat: B/C may add fields
    "fields": {
        "schema_version": {"type": "str", "const": "stage-result.v1"},
        "status": {"type": "str",
                   "enum": ["running","completed","failed","timed_out","cancelled"]},
        "exit_code": {"type": "int", "nullable": True},
        "attempt": {"type": "int"},
        "error": {"type": "object", "nullable": True,
                  "fields": {"type": {"type": "str"}, "message": {"type": "str"}}},
        "timing": {"type": "object", "fields": {
            "started_at": {"type": "str"}, "finished_at": {"type": "str"},
            "duration_sec": {"type": "number"}}},
        "inputs":  {"type": "array", "items": {"type": "object"}},
        "outputs": {"type": "array", "items": {"type": "object"}},
        "provenance_files": {"type": "array", "items": {"type": "object"}},
        "producer": {"type": "object",
                     "fields": {"name": {"type": "str"}, "version": {"type": "str"}}},
        # metrics/software/reference_db present but unconstrained in A
    },
}
SCHEMAS = {"stage-result.v1": STAGE_RESULT_V1}
```

`validate(payload, schema_name) -> list[str]` walks the schema: checks required
keys present, `type` (`str|int|number|bool|object|array`), `const`, `enum`,
`nullable` (allows explicit `None`), nested `fields`, and `array` `items`; unknown
keys allowed when `allow_unknown`. Returns human-readable error strings (empty ⇒
valid). `validate_stage_result(payload)` is the convenience wrapper.

### `run_command` integration (`runtime/runner.py`)

Guarded by `run_dir is not None` (unchanged). In the success branch (after the
existing `write_results_manifest`) and in the failure and timeout paths, call
`write_stage_result(...)` with the mapped `status`/`exit_code`/`error`, wrapped so
a writer exception cannot suppress the pipeline's real error. `microsuite-results.json`
writing is unchanged. Wire `denoise` (`methods/denoise.py`) as the reference:
populate `CommandLog.stage`, `artifact_inputs/outputs`, and a `provenance_files`
reference to `dada2_denoise_manifest.json`.

## Testing (all offline / unit)

`tests/test_metadata_validate.py`, `tests/test_metadata_stage_result.py`,
`tests/test_metadata_redact.py`, plus additions to `tests/test_runtime_runner.py`:

1. **Validator**: valid payload → `[]`; missing required → named error; wrong
   type → error; bad `status` enum → error; wrong `schema_version` const → error;
   nullable `exit_code: None` accepted; **unknown top-level field accepted**
   (forward-compat); malformed (non-dict) input → error, not crash.
2. **Redaction**: `redact_params` masks `auth_token`/`api_key`/nested; leaves
   benign keys; `redact_command` masks `--token X`, `--api-key=X`, and bare
   values equal to captured secrets; `redact_text` masks secrets in a message.
3. **Writer success**: envelope at `stage-results/<stage>--<backend>--attempt-1.json`;
   validates clean; `outputs` enriched (`exists`, `bytes` from a real temp file);
   directory output → `bytes: null`; relative path for in-`run_dir` output,
   `external: true` for an outside input; `metrics/software` empty,
   `reference_db` null.
4. **Writer failure/timeout**: `status: failed` with `exit_code`, structured
   `error`; `status: timed_out` with `exit_code: null`; partial existing output
   still recorded; original exception preserved (the pipeline still raises).
5. **Retry**: two writes for the same stage+backend → `attempt-1.json` and
   `attempt-2.json`, distinct `stage_run_id`, `attempt` 1 then 2; neither
   overwrites the other.
6. **Atomic replacement**: no `*.tmp*` left behind after a successful write; the
   target is only ever the fully-rendered JSON (assert parseable + no temp files).
7. **Identifiers**: env `MICROSUITE_WORKFLOW_RUN_ID`/`MICROSUITE_DATASET_ID` land
   in the payload; absent → `null`; `run_id` falls back to `run_dir.name`.
8. **`run_command`**: success writes both `microsuite-results.json` (unchanged)
   and a `stage-result` envelope; a non-zero command writes a `failed` envelope
   **and** still raises `MicrobiomeSuiteError`; `run_dir=None` writes nothing.

## Success criteria

1. `microsuite/metadata/` exists with schema registry, dependency-free validator,
   `StageResult`/`Artifact`/`ProvenanceFile` models, redactor, and an atomic
   writer producing unique per-attempt files under `run_dir/stage-results/`.
2. `run_command` emits a validating `stage-result.v1.json` on success, failure,
   and timeout, from explicit `CommandLog` declarations (no folder/argv guessing),
   with `run_id` + `stage_run_id` + optional workflow identity, redacted
   `params`/`command`, and empty `metrics`/`software`/`reference_db` slots.
3. `microsuite-results.json` is byte-for-byte unchanged by A.
4. `denoise` declares typed artifacts + a provenance reference as the worked
   example. Full suite + `ty check` + ruff + format green.

## Open questions / follow-ups (not blocking A)
- **B**: `software`/`reference_db`/`tool_versions` capture + `resolved_config.json`
  (full resolved config incl defaults, secrets excluded via this redactor).
- **C**: `metrics` with units (`retainedReads`, `meanShannon`, `taxaCount`, …)
  into the envelope and the aggregate; lets Microboard stop folder-guessing.
- Separate `workflow-run.v1.json` (topology: nodes/edges/branches/method
  selection) — a later sub-project; Microboard reads it for the DAG.
- Convert the remaining stages (taxonomy, diversity, diffab, …) to declare typed
  artifacts + provenance references; fold sub-project I's `<backend>_container.json`
  in as a `ProvenanceFile` once both land.
