# Metadata sub-project B — software/refDB versions + resolved_config — Design

- **Date:** 2026-07-14
- **Status:** Approved (design locked via decision questions), implementing
- **Origin:** Fills the empty `software`/`reference_db` slots left by sub-project
  **A** and adds the `resolved_config.json` writer the end-to-end runner needs.
  See [[microboard-metadata-integration]]. Follows A (`688453d`).

## Decisions (locked)
1. **Software versions are declaratively lifted, never live-probed.** A stage sets
   `software` from data it already produced. No `<tool> --version` calls.
2. **`resolved_config.json` is a writer util + schema for the orchestrator.**
   No single microsuite command owns the whole config; the end-to-end runner
   (oral repo) aggregates and calls it. Per-stage resolved params already live in
   `stage-result.params`.
3. **refDB: API + schema now; wire `tax_classify` as a stage later.** B adds
   `set_reference_db()` + the documented shape; denoise's software lift is the
   worked example.

## Scope of B
- `StageRecord.set_software(mapping)` and `set_reference_db(info)` — populate the
  A envelope's `software` / `reference_db` slots (default stays `{}` / `null`, so
  A's fixtures/tests are unchanged).
- `denoise` (dada2-r): after the R params file is written, lift
  `{"microsuite": {version}, "dada2": {version}, "R": {version}}` into `software`
  (best-effort: absent/unparseable params → leave `software` empty).
- `metadata/config.py`: `write_resolved_config(run_dir, config, *, name=…)` —
  redact secrets (reuse `redact_params`), wrap as `resolved-config.v1`, validate,
  atomic-write.
- Schema: add `resolved-config.v1` to the internal registry + a published
  `_schema/resolved-config.v1.schema.json` (draft 2020-12) + valid/invalid
  fixtures + parity test. `stage-result.v1` is **unchanged** (its `software`/
  `reference_db` were already permissive objects — forward-compatible).

### Out of scope for B
- Wrapping `tax_classify`/other stages (reference_db population end-to-end) — a
  follow-up once more stages are stage-wrapped.
- Metrics with units — sub-project **C**.
- Live version probing of any kind.

## Design

### Envelope slots (no schema change)
`stage-result.v1` keeps `software` (`object`, permissive) and `reference_db`
(`object|null`, permissive). B documents the shapes it writes:
- `software`: `{ "<tool>": {"version": str, …}, … }` (e.g. `dada2`, `R`,
  `microsuite`).
- `reference_db`: `{ "name": str, "version": str, "build_target": str,
  "checksum": str, "provider": str }` or `null`.

### StageRecord additions (`metadata/stage.py`)
- `self._software: dict[str, Any] = {}`, `self._reference_db: dict | None = None`.
- `set_software(mapping)` → `self._software.update(mapping)` (JSON-safe-coerced).
- `set_reference_db(info)` → store (JSON-safe-coerced) or `None`.
- `to_payload` emits `"software": self._software` and
  `"reference_db": self._reference_db` (instead of the A constants). Empty stays
  `{}` / `null`.

### denoise software lift (`methods/denoise.py`)
In `_declare_dada2_stage_outputs` (already runs after `_emit_dada2_manifest`),
best-effort read `dada2_r_params.json` (`dada2_manifest.read_r_params`) and, when
present, `stage.set_software({"microsuite": {"version": _MICROSUITE_VERSION},
"dada2": {"version": r_params.get("dada2_version")}, "R": {"version":
r_params.get("r_version")}})`, dropping keys whose version is `None`. A missing
params file → no software (stays `{}`); never fails the run.

### resolved-config (`metadata/config.py`, `resolved-config.v1`)
```json
{
  "schema_version": "resolved-config.v1",
  "generated_at": "2026-07-14T09:12:03Z",
  "producer": {"name": "microsuite", "version": "0.1.0"},
  "config": { … redacted resolved configuration, defaults included … }
}
```
`write_resolved_config(run_dir, config, *, name="resolved_config.json") -> Path`:
`masked, _ = redact_params(config)` (secrets → `***`, JSON-safe-coerced); build
the envelope; `validate(payload, "resolved-config.v1")` → raise
`MicrobiomeSuiteError` on invalid; atomic temp+rename (reuse the stage writer's
`_atomic_write`). Internal schema: required `schema_version` (const),
`generated_at` (rfc3339), `producer` (`{name,version}`), `config` (object,
allow_unknown); `allow_unknown: true` at top. Published JSON Schema mirrors it.

## Testing (offline)
- `set_software`/`set_reference_db` populate the envelope; empty defaults keep A's
  behaviour; a populated envelope still `validate_stage_result == []`.
- denoise integration: the completed envelope's `software` has `dada2`/`R`/
  `microsuite` versions lifted from a faked `dada2_r_params.json`; a missing params
  file → `software == {}`.
- `write_resolved_config`: masks `auth_token`/nested secrets, includes defaults,
  writes a valid `resolved-config.v1`, atomic (no `*.tmp*`), invalid config-envelope
  path raises.
- schema-contract: `resolved-config.v1` valid/invalid fixtures pass the Python
  validator and the published JSON Schema (jsonschema parity, skip if absent).

## Success criteria
1. Stages can populate `software`/`reference_db`; denoise lifts dada2/R/microsuite
   versions as the worked example; A's empty-slot behaviour is preserved.
2. `write_resolved_config` emits a redacted, validating `resolved-config.v1`
   snapshot (defaults in, secrets out) for the orchestrator.
3. Full suite + `ty` + `ruff` + `format` green; `microsuite-results.json` and the
   `stage-result.v1` structural schema unchanged.
