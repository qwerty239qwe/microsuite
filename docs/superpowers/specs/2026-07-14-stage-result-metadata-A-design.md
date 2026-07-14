# Structured run metadata for Microboard — Sub-project A (schema + envelope foundation) — Design

- **Date:** 2026-07-14
- **Status:** Approved (design, revised after architecture review), pending
  implementation plan
- **Origin:** Microboard integration. The neighbouring `microboard` dashboard
  imports `microsuite-results.json` (schema `microsuite-results.v1`) and otherwise
  *guesses* stage metrics by scanning result folders. We want microsuite to emit a
  first-class, versioned, **per-stage** metadata envelope so Microboard stops
  filename-guessing and can reliably compare benchmark runs. Decomposition:
  **A** (this — schema registry, validator, models, stage-execution boundary +
  writer), **B** (software/refDB versions + `resolved_config.json`), **C**
  (metrics with units).

## Key revision (architecture review)

The stage boundary is **not** the subprocess boundary. `run_command` runs one
subprocess; a biological stage (e.g. DADA2) then runs *output validation, QC
summary, and provenance writes in Python afterwards* (see
`methods/denoise.py:729` `_run()` and the later manifest write near
`denoise.py:745`). Writing the envelope inside `run_command` would (a) report
`completed` before Python post-processing fails, (b) record the provenance file as
missing (it is written after), (c) miss QC-derived outputs, (d) turn a
multi-subprocess stage into several files that look like retries, and (e) produce
no `failed` envelope for a Python post-processing failure.

**Therefore A introduces a stage-level boundary — a `stage_execution(...)` context
manager — that finalizes exactly one envelope after the whole stage
succeeds or fails.** `run_command` stays generic: it keeps writing its
command-level artifacts (`command.txt`, `events.jsonl`, `run.json`, logs, the
unchanged `microsuite-results.json`) and, *when a stage is active*, contributes
subprocess details (exit code, duration, timeout flag) to that stage. `run_command`
does **not** write stage-result envelopes.

## Scope of A

- `microsuite/metadata/` package: versioned **schema registry**, a dependency-free
  **MicroSuite-specific** validator (recursive + status invariants), typed models
  (`Artifact`, `ArtifactCount`, `ProvenanceFile`, `StageError`), a mutable
  `StageRecord` accumulator, a shared **secret redactor**, and an **atomic**
  writer producing concurrency-safe, unique per-attempt files.
- A `stage_execution(...)` context manager that owns envelope finalization on
  **success, failure, and timeout** (from `finally`, preserving the original
  exception), from **explicit** artifact declarations — never folder/argv guessing.
- `run_command` integration: feed subprocess details to the active stage via a
  `ContextVar`; no behavioural change when no stage is active.
- **One worked reference stage: `denoise` (dada2-r)** — wrapped so its subprocess,
  output validation, QC, and `dada2_denoise_manifest.json` provenance all land in
  a single finalized envelope.
- `metrics`, `software`, `reference_db` slots ship **empty/null** in A.

### Out of scope for A
- Populating `metrics` (C), `software`/`reference_db`/`tool_versions` (B), and
  writing `resolved_config.json` (B).
- Any change to `microsuite-results.json` — it **remains backward-compatible and
  unchanged during A**.
- Workflow topology (nodes/edges/branches/method selection) — a later
  `workflow-run.v1.json`. A `stage-result` describes one execution.
- Converting stages other than `denoise` to declare typed artifacts — incremental
  follow-up.
- Redacting the *other* run-directory files (see Limitation below).
- Microboard-side ingestion changes.

## Verified current state

- `runtime/results.py` writes `microsuite-results.json` (`executions[]`,
  `artifacts[]`, `producer`) from `run_command`, **only on exit 0**, only when
  `run_dir` set.
- `runtime/runner.py` `run_command` substitutes an empty `CommandLog()` when none
  is passed (`runner.py:52`); some callers/tests invoke it with no log
  (`tests/test_runtime_runner.py:120`). So `stage`/`task` cannot be assumed present
  at the `run_command` layer — another reason the envelope boundary is the stage,
  not `run_command`.
- No schema registry, validator, envelope, redaction, or failure-path metadata.
- No `jsonschema`/`pydantic` dependency → dependency-free validator by design.

## Design

### Stage boundary (`metadata/stage.py`) — finding #1, #3

```python
_ACTIVE: ContextVar[StageRecord | None] = ContextVar("active_stage", default=None)

@contextmanager
def stage_execution(
    run_dir: Path | None, *, stage: str, task: str | None = None,
    backend: str | None = None, params: Mapping[str, Any] | None = None,
    inputs: Iterable[Artifact] = (), outputs: Iterable[Artifact] = (),
    provenance_files: Iterable[ProvenanceFile] = (),
) -> Iterator[StageRecord]:
    record = StageRecord(run_dir, stage=stage, task=task, backend=backend,
                         params=dict(params or {}), inputs=list(inputs),
                         outputs=list(outputs), provenance=list(provenance_files))
    token = _ACTIVE.set(record)
    try:
        yield record            # stage adds outputs/provenance as it produces them
    except BaseException as exc:
        record.mark_failure(exc)      # status failed|timed_out; exit_code; StageError
        _publish(record, on_failure=True)   # best-effort; original exception re-raised
        raise
    else:
        record.mark_success()         # status completed
        _publish(record, on_failure=False)  # strict: invalid/write failure raises
    finally:
        _ACTIVE.reset(token)
```

`StageRecord` (mutable) exposes `add_output(Artifact)`, `add_provenance(ProvenanceFile)`,
`add_input(Artifact)`, `note_subprocess(command, exit_code, duration_sec, *, timed_out)`,
and `to_payload() -> dict`. `stage` is **required** here (schema-required, and the
stage always knows it) — resolving finding #3 without touching generic `run_command`.

`run_command` change: after a subprocess completes (or times out) and `run_dir`
is set, if `_ACTIVE.get()` is not `None`, call
`record.note_subprocess(command, exit_code, duration_sec, timed_out=...)`. Nothing
else changes; with no active stage, `run_command` behaves exactly as today. The
last successful subprocess supplies the envelope's `command`/`exit_code`; on
timeout the record is marked before `run_command` re-raises.

### Status & invariants (finding #4)

`mark_success()` → `completed`. `mark_failure(exc)`: if any noted subprocess timed
out → `timed_out` (`exit_code = null`); else `failed` (`exit_code =` the failing
subprocess code if one failed, else `null` for a Python post-processing failure).
`StageError = {"type": exc.__class__.__name__, "message": <safe, redacted, truncated>}`.

Status invariants enforced by the validator:
- `completed` → `error is null` **and** `exit_code in (0, null)` (null when the
  stage ran no subprocess — a deliberate relaxation of "must be 0" so pure-Python
  stages validate).
- `failed` → `error` is a non-null object.
- `timed_out` → `exit_code is null` and `error` non-null.

### File layout & concurrency-safe naming (findings #1, #2)

One file per stage attempt, never overwritten, under a run-dir subdirectory:

```
run_dir/
  stage-results/
    denoise--dada2-r--attempt-1--stage-run-9f3c1a2b7d40.json
    taxonomy--silva--attempt-1--stage-run-1b77e0c4a219.json
  microsuite-results.json                 # legacy aggregate, unchanged
```

- Name: `<stage>--<backend>--attempt-<N>--<stage_run_id>.json`; components
  slugified (`[^a-z0-9]+`→`-`, lower-cased; `backend` `None` → `none`).
- **Uniqueness comes from `stage_run_id`** (a UUID), so two concurrent processes
  cannot collide even if both compute `attempt-1`. `attempt` (`1 + count(existing
  matching files)`) is **informational/best-effort**, not a uniqueness guarantee.
- Invalid-on-publish payloads go to a **diagnostic** name
  (`<same>.invalid.json`), never the normal path (finding #6).
- `run_dir is None` → no envelope (unchanged behaviour).

### Identifiers (findings #1, #2)

| Field | Source | Notes |
|---|---|---|
| `run_id` | env `MICROSUITE_RUN_ID`, else `run_dir.name` | stable per workflow execution |
| `stage_run_id` | generated `stage-run-<uuid4hex12>` | unique per attempt; also in filename |
| `attempt` | `1 + count(existing files)` | informational; `>= 1` |
| `workflow_id` | env `MICROSUITE_WORKFLOW_ID` | optional; `null` when absent |
| `workflow_run_id` | env `MICROSUITE_WORKFLOW_RUN_ID` | optional; separates two DADA2 runs in different benchmarks |
| `dataset_id` | env `MICROSUITE_DATASET_ID` | optional |

`metadata/context.py:workflow_context(overrides=None)` reads env once (overridable
for tests).

### Models (`metadata/models.py`) — findings #3, #4, #5

Frozen dataclasses. **Declared artifact paths MUST be absolute** (the caller passes
absolute `Path`s); the writer serialises them **relative to `run_dir` when the
target is under it**, else keeps them absolute and sets `external: true`. `external`
applies uniformly to inputs, outputs, and provenance (DADA2's provenance is written
beside `output_stats`, which may sit outside `run_dir`).

```python
@dataclass(frozen=True)
class ArtifactCount:
    value: int          # >= 0
    unit: str           # what is counted: "samples","features","reads",...

@dataclass(frozen=True)
class Artifact:
    label: str
    path: str | Path              # absolute when declared; serialised rel/abs by writer
    format: str | None = None
    kind: str | None = None
    count: ArtifactCount | None = None   # declared by the stage; writer never counts
    external: bool = False               # set by writer if outside run_dir
    # writer fills on serialise: bytes (int>=0|null), exists (bool)

@dataclass(frozen=True)
class ProvenanceFile:
    kind: str
    path: str | Path              # absolute when declared
    # writer fills: external, exists

@dataclass(frozen=True)
class StageError:
    type: str
    message: str        # safe, redacted, truncated
```

Writer **enrichment** (from the absolute path, before relativising): `exists` and,
for regular files, `bytes = stat().st_size` (`null` for directories/missing). It
never invents artifacts and never derives `count`.

`CommandLog` is **unchanged** (legacy dict `inputs/outputs/params` still feed
`microsuite-results.json`). Typed artifacts live on the stage, not `CommandLog`.

### Envelope schema `stage-result.v1`

```json
{
  "schema_version": "stage-result.v1",
  "run_id": "results_erp120510_end_to_end",
  "stage_run_id": "stage-run-9f3c1a2b7d40",
  "workflow_id": "oral-standard",
  "workflow_run_id": "oral-standard--erp120510--001",
  "dataset_id": "ERP120510",
  "attempt": 1,
  "stage": "denoise", "task": "denoise", "backend": "dada2-r",
  "status": "completed", "exit_code": 0, "error": null,
  "timing": {"started_at":"2026-07-14T09:12:03Z","finished_at":"2026-07-14T09:12:15Z","duration_sec":12.34},
  "command": ["Rscript","denoise.R","--token","***"],
  "params": {"trunc_len_f":240,"auth_token":"***"},
  "inputs":  [{"label":"reads","path":"/data/input","format":"directory","external":true,"exists":true}],
  "outputs": [{"label":"feature table","path":"feature-table.tsv","format":"tsv","kind":"feature_table","count":{"value":1842,"unit":"features"},"bytes":928104,"external":false,"exists":true}],
  "provenance_files": [{"kind":"dada2_manifest","path":"dada2_denoise_manifest.json","external":false,"exists":true}],
  "metrics": {}, "software": {}, "reference_db": null,
  "producer": {"name":"microsuite","version":"0.x.y"}
}
```

Timestamps are UTC RFC 3339 with `Z`. `command`/`params` share one redactor.

### Failure & timeout semantics (finding #5)

On `failed`/`timed_out`, the envelope still records timing, `params`, redacted
`command`, declared `inputs`, declared `outputs` (each enriched with `exists`/`bytes`
so partial outputs are visible), `exit_code` when available, structured `error`,
and the `provenance_files` that were actually written (`exists: true`). Envelope
publication runs from the `except` path and **re-raises the original exception**;
any error inside publication is swallowed there (best-effort) so the pipeline's real
failure propagates unchanged.

### Publish contract (finding #6)

`_publish(record, *, on_failure)`:
1. `payload = record.to_payload()` — enrich, relativise, JSON-safe-coerce, redact.
2. `errors = validate_stage_result(payload)`.
3. **On the success path (`on_failure=False`):** if `errors`, write the payload to
   the **diagnostic** path and raise `MicrobiomeSuiteError` — never publish an
   invalid document at a normal path, and never report success with broken
   metadata. If valid, atomic-write to the normal path; a write failure raises.
4. **On the failure path (`on_failure=True`):** if `errors`, write to the
   diagnostic path and warn; else atomic-write to the normal path. Any exception
   here is swallowed (warn only) to preserve the original pipeline exception.
- **Atomic write:** `tmp = target.with_name(target.name + f".tmp.{os.getpid()}")`;
  render JSON; `os.replace(tmp, target)`. Microboard never reads a partial file.

### Secret redaction (`metadata/redact.py`) — finding #7

One mechanism for `params` and `command`:
- `SENSITIVE_KEY_RE = (?i)(token|secret|password|passwd|api[-_]?key|credential|auth)`.
- `redact_params(mapping)` → deep copy; sensitive keys → `"***"` (recurse
  dicts/lists). Also **JSON-safe-coerce** every value: `Path`→str, `Enum`→`.value`,
  `tuple`/`set`→list, other non-JSON types→`str(...)`.
- Capture the set of redacted secret *values*; ignore empty/whitespace ones.
- `redact_command(argv, secrets)` → mask value after a sensitive flag
  (`--token X`→`***`), inline `--api-key=X`→`--api-key=***`, and any bare arg equal
  to a captured secret.
- `redact_text(s, secrets)` → replace captured secrets in free text, **longest-first**,
  skipping secrets shorter than 4 chars (avoid corrupting unrelated text). Used for
  `error.message`.

**Limitation (documented, finding #7):** A redacts only the stage-result envelope.
The existing `command.txt`, `events.jsonl`, and `run.json` are still written from
raw values, so a run directory is **not** wholesale safe to share after A. Redacting
those is a follow-up.

### Schema registry & validator (`metadata/schemas.py`, `metadata/validate.py`) — findings #4, #7

A **MicroSuite-specific** schema format — a small subset inspired by JSON Schema,
**not** a JSON Schema validator (stated in the module docstring). Supports:
`required`, `allow_unknown`, per-field `type` (`str|int|number|bool|object|array`),
`const`, `enum`, `nullable`, `min` (numeric), nested `fields`, `array` `items`, and
a schema-level `invariants` hook (callables returning error strings) for the
status/exit_code/error cross-field rules.

`stage-result.v1` **required** (every consistently emitted top-level field):
`schema_version, run_id, stage_run_id, attempt, stage, task, backend, status,
exit_code, error, timing, command, params, inputs, outputs, provenance_files,
metrics, software, reference_db, producer`. Nullable values: `backend, exit_code,
error, command, reference_db` (and the optional `workflow_id/workflow_run_id/
dataset_id` when present). Recursive item schemas validate `Artifact`
(`label:str`, `path:str` required; `format/kind` nullable str; `count` →
`ArtifactCount{value:int min0, unit:str}`; `bytes` int min0 nullable; `external`
bool; `exists` bool), `ProvenanceFile{kind:str, path:str, external:bool, exists:bool}`,
`StageError{type:str, message:str}`, `Timing{started_at:str, finished_at:str,
duration_sec:number min0}`, `Producer{name:str, version:str}`. Numeric
constraints: `attempt>=1`, `bytes>=0`, `duration_sec>=0`, `count.value>=0`.
`allow_unknown: true` at top level and on nested objects → B/C add fields without
breaking older Microboard (finding: forward-compat). `metrics`/`software` are
`object` (may be empty); unconstrained in A beyond type.

`validate(payload, schema_name) -> list[str]` returns human-readable errors (`[]`
⇒ valid); non-dict input → a single error, never a crash.
`validate_stage_result(payload)` is the wrapper (applies the invariants hook).

## Testing (all offline / unit)

New: `tests/test_metadata_validate.py`, `test_metadata_redact.py`,
`test_metadata_stage.py`; plus `tests/test_runtime_runner.py` additions; plus a
`denoise` integration test.

1. **Validator:** valid → `[]`; each missing required → named error; wrong type;
   bad `status` enum; wrong `schema_version` const; nullable `exit_code:null` ok;
   `attempt:0`/`bytes:-1`/`duration_sec:-1` rejected; `{"outputs":[{}]}` rejected
   (Artifact requires `label`,`path`); **unknown top-level & nested field
   accepted** (forward-compat); invariants: `completed`+`error!=null` rejected,
   `failed`+`error==null` rejected, `timed_out`+`exit_code!=null` rejected; non-dict
   input → error not crash.
2. **Redaction:** `redact_params` masks `auth_token`/`api_key`/nested, JSON-safe-
   coerces `Path`/`Enum`/`tuple`; `redact_command` masks `--token X`,
   `--api-key=X`, bare captured secret; empty secret ignored; `redact_text`
   longest-first and skips <4-char secrets.
3. **Writer success:** envelope at
   `stage-results/denoise--dada2-r--attempt-1--stage-run-*.json`; validates clean;
   output enriched (`exists`,`bytes` from a real temp file); directory output →
   `bytes:null`; in-`run_dir` path serialised relative, outside path
   `external:true` + absolute; `metrics/software` empty, `reference_db` null;
   **no `*.tmp*` left behind**.
4. **Failure/timeout:** `status:failed` with `exit_code` + structured `error`;
   `status:timed_out` with `exit_code:null`; a partial existing output still
   recorded; the original exception still propagates.
5. **FileNotFoundError / launch error** inside the stage → a `failed` envelope is
   still published, exception preserved.
6. **Multiple subprocesses in one stage** → exactly **one** envelope (not several
   "retries").
7. **Retry:** two `stage_execution` runs for the same stage+backend → two files
   with distinct `stage_run_id`, `attempt` 1 then 2; neither overwrites.
8. **Concurrent allocation:** two records that both compute `attempt-1` write to
   distinct filenames (different `stage_run_id`) — no overwrite.
9. **Missing task/stage:** `stage` required by `stage_execution`; `run_command`
   with no log / no active stage writes **no** envelope (existing test still
   passes) and still writes `microsuite-results.json`.
10. **Success-path publish strictness:** an intentionally invalid payload on the
    success path raises and lands only at `*.invalid.json`; a writer failure after
    a successful stage raises `MicrobiomeSuiteError`.
11. **Relative paths with custom `cwd`:** an output declared absolute but produced
    under a `run_command(cwd=...)` still serialises relative to `run_dir` when
    under it, else `external:true`.
12. **DADA2 reference:** the `denoise` envelope is finalized **after** output
    validation, QC, and `dada2_denoise_manifest.json` are written — the provenance
    file is `exists:true`, QC outputs are present, and a forced Python
    post-processing failure yields `status:failed`.
13. **Documented log limitation:** a test asserts raw `command.txt`/`run.json`
    still contain unredacted values (documents that A does not sanitise them).

## Success criteria

1. `microsuite/metadata/` provides the schema registry, recursive+invariant
   dependency-free validator, typed models, redactor, and atomic writer with
   concurrency-safe unique filenames.
2. `stage_execution(...)` finalizes exactly one validating `stage-result.v1.json`
   per stage on success, failure, and timeout, from explicit declarations, with
   dual identifiers + optional workflow identity, redacted `command`/`params`, and
   empty `metrics`/`software`/`reference_db`. Invalid documents never reach a
   normal path; success is never reported with unwritten/broken metadata.
3. `run_command` stays generic (contributes subprocess details to the active stage
   only) and `microsuite-results.json` is byte-for-byte unchanged.
4. `denoise` is wrapped as the worked stage, with provenance/QC inside the
   envelope boundary. Full suite + `ty check` + ruff + format green.

## Open questions / follow-ups (not blocking A)
- **B**: `software`/`reference_db`/`tool_versions` + `resolved_config.json` (defaults
  in; secrets excluded via this redactor).
- **C**: `metrics` with units into the envelope + aggregate.
- Redact `command.txt`/`events.jsonl`/`run.json` so a run dir is shareable.
- Separate `workflow-run.v1.json` (topology).
- Convert remaining stages to typed declarations; fold sub-project I's
  `<backend>_container.json` in as a `ProvenanceFile`.
