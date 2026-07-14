# Metadata sub-project C — metrics with units — Design

- **Date:** 2026-07-14
- **Status:** Approved (design locked via decision questions), implementing
- **Origin:** Fills the `metrics` slot left by sub-project **A** so Microboard
  reads real, unit-tagged metrics from the stage-result envelope instead of
  folder-guessing. Follows A (`688453d`) + B (`c64982f`). See
  [[microboard-metadata-integration]].

## Decisions (locked)
1. **Metrics live in the stage-result envelope only.** `microsuite-results.json`
   stays byte-stable (A/B froze it); the envelope is Microboard's new primary path.
2. **Wrap denoise + diversity + taxonomy stage boundaries.** denoise emits real
   metrics from its native QC summary now. `diversity_calc` / `tax_classify` output
   QIIME 2 `.qza` artifacts, so their `meanShannon`/`taxaCount`/`assignedReads`
   need a `.qza` extractor — **deferred to a focused follow-up**. C wraps their
   boundaries (envelope + subprocess capture + declared output) with `metrics: {}`.

### Out of scope for C
- `.qza` metric extraction for diversity/taxonomy (follow-up).
- Metrics for the exotic taxonomy backends (emu/bracken/metaphlan/kraken2).
- Any change to `microsuite-results.json` or the `stage-result.v1` structural
  schema (`metrics` is already a permissive object — forward-compatible).

## Design

### Metrics slot (no schema change)
`stage-result.v1` keeps `metrics` as a permissive `object`. C documents + writes
the shape `{ "<name>": {"value": number, "unit": str}, … }`.

### StageRecord additions (`metadata/stage.py`)
- `self._metrics: dict[str, Any] = {}`.
- `add_metric(name, value, unit)` → `self._metrics[name] = {"value": value,
  "unit": unit}` (value JSON-safe-coerced; NaN/inf rejected → skip, best-effort).
- `set_metrics(mapping)` → merge a `{name: {"value","unit"}}` mapping.
- `to_payload` emits `"metrics": self._metrics` (was `{}` in A). Empty stays `{}`.

### denoise metrics (`methods/denoise.py`)
In `_declare_dada2_stage_outputs` (best-effort, never fails the run): compute
`summary = dada2_qc.summarize_dada2_stats(output_stats)` and emit:
- `input_reads` = `overall.input` (`unit: "reads"`),
- `filtered_fraction` = `overall.filtered_frac` (`unit: "fraction"`),
- `merged_fraction` = `overall.merged_frac` when paired (`unit: "fraction"`),
- `nonchimeric_fraction` = `overall.nonchim_frac` (`unit: "fraction"`),
- `nonchimeric_reads` = `round(overall.input * overall.nonchim_frac)`
  (`unit: "reads"`).
A missing/unparseable stats file → no metrics (stays `{}`).

### diversity + taxonomy stage boundaries
Wrap each dispatcher body in `stage_execution` (run_dir = `output.parent`):
- `methods/diversity_calc.py` `diversity_calc(...)`: `stage="diversity"`,
  `backend=<backend>`, params = a small resolved dict (metric, threads), one
  declared output (`output`, `kind="alpha_diversity"`, **`required=False`** —
  these lack denoise's `validate` flag and tests stub the subprocess without
  producing the `.qza`; existence is still recorded). `metrics: {}`.
- `methods/tax_classify.py` `tax_classify(...)`: `stage="taxonomy"`,
  `backend=<backend>`, one declared output (`output`, `kind="taxonomy_table"`,
  `required=False`). `metrics: {}`. (refDB capture via B's `set_reference_db`
  stays a follow-up — `resolve_classifier` returns only a path.)
The subprocess is auto-captured by `run_command` via the `ContextVar`. Errors
inside → failed envelope + re-raise (same as denoise).

## Testing (offline)
- `add_metric`/`set_metrics` populate `metrics`; empty default unchanged; a
  populated envelope still `validate_stage_result == []`; NaN value skipped.
- denoise integration: the completed envelope's `metrics` has `input_reads`
  (unit reads), `nonchimeric_fraction` (unit fraction), `nonchimeric_reads`
  from a faked stats file; missing stats → `metrics == {}`.
- diversity + taxonomy: a stubbed-`subprocess.run` run emits one envelope under
  `output.parent/stage-results/` with the declared output, the captured
  subprocess (`completed`), `metrics == {}`; a non-zero subprocess → `failed`
  envelope + raise.

## Success criteria
1. Stages can populate `metrics` with `{value, unit}`; denoise emits real
   retention metrics; A's empty-slot behaviour is preserved.
2. diversity + taxonomy emit stage-result envelopes (boundary wrapped, subprocess
   captured, output declared); their unit metrics are a documented follow-up.
3. Full suite + `ty` + `ruff` + `format` green; `microsuite-results.json` and the
   `stage-result.v1` structural schema unchanged.
