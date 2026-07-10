# DADA2 parameter provenance manifest (Round-2 A) — Design

- **Date:** 2026-07-10
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** Codex's ERP120510 friction list, complaint **#7** (DADA2 parameter
  provenance too weak — run dirs should include a resolved machine-readable
  manifest of the *effective* params after defaults + CLI/config overrides, plus
  tool versions). Sub-project **A** of the DADA2 round-2 roadmap (see
  [[dada2-improvement-roadmap]]); **B** (QC warnings) merged, **C** (sweep) follows.

## Scope

A writes, on every successful `dada2-r` run, a `dada2_denoise_manifest.json`
beside the outputs recording the **effective** DADA2 parameters (after R-side
defaults resolve) plus `dada2`/R versions and the wrapper-level run facts
(runtime, image, microsuite version, mode, threads, paths). The R backend is the
authoritative source for the resolved dada2 params + versions; Python merges in
the wrapper facts. Works identically for `--runtime local` and `docker`.

### Out of scope for A
- Parameter-sensitivity sweep (#10) — **C**.
- Changing any DADA2 default — A only *records* effective values.
- Provenance for non-dada2-r backends, or for failed runs (manifest is written
  on success only, matching the existing results-manifest convention).
- Replacing the existing per-run `run.json`/`events.jsonl` (left untouched; the
  new manifest is the resolved-provenance artifact #7 asks for).

## Verified context

- `run_command(run_dir=...)` already writes `run.json` from `CommandLog.params`,
  but `_dada2_log_params(**params)` (`denoise.py:820`) **drops every `None`**, so
  only explicitly-set params are recorded — never the resolved defaults. And
  `run.json` only appears when a `--run-dir` is passed.
- The R script (`methods/r/dada2_denoise.R`) resolves each default inline via
  `value_after(flag, default)` — e.g. `maxEE` default `"2"`/`"2"`, `truncLen`
  `"0"`, `minOverlap` `"12"`, `maxMismatch` `"0"`, `truncQ` `"2"`, pooling
  `"independent"`, chimera `"consensus"`, `n_reads_learn` `"1000000"`. Only R
  knows these effective values and the `dada2`/R versions.
- `runtime`/`image` are Python-only concepts the R side never sees, so the
  manifest must be assembled on the Python side from both sources.
- B established the "write a sidecar beside `output_stats`" pattern
  (`dada2_qc_summary.{json,tsv}`) and the pure-helper-module pattern
  (`methods/dada2_qc.py`). A mirrors both.

## Design

### Component 1 — R emits resolved params (`dada2_denoise.R`)

Hoist the scattered `value_after(flag, default)` calls for the tunable params
into named `resolved_*` variables near the top (single source of truth), and use
those variables in `filterAndTrim`/`dada`/`mergePairs` (both paired and single
branches). This removes the current duplication of default literals across the
two branches and makes the manifest report exactly what was used.

Add a `--params-out PATH` flag. After a successful run, write a **flat JSON**
object to that path via a small dependency-free writer (no `jsonlite`
assumption — hand-roll from `paste`, quoting strings and emitting numbers/bools
bare/`null`). Keys (effective values actually used):

```
mode ("paired"|"single"), trim_left_f, trim_left_r, trunc_len_f, trunc_len_r,
trim_left, trunc_len, max_ee_f, max_ee_r, max_ee, trunc_q, max_n, rm_phix,
pooling_method, chimera_method, min_fold_parent_over_abundance, allow_one_off,
n_reads_learn, min_overlap, max_merge_mismatch, trim_overhang,
dada2_version, r_version
```

Single-mode-only keys (`trim_left`, `trunc_len`, `max_ee`) and paired-only keys
(`*_f`/`*_r`, `min_overlap`, `max_merge_mismatch`, `trim_overhang`) are emitted
for the branch that ran; the other branch's keys are `null`. Versions from
`as.character(packageVersion("dada2"))` and `R.version.string`. If `--params-out`
is absent, R skips the write (keeps the script runnable standalone).

### Component 2 — `methods/dada2_manifest.py` (pure, offline-testable)

Mirrors `dada2_qc.py`:

```python
def read_r_params(path: Path) -> dict:
    """Parse the flat JSON the R backend wrote; raise MicrobiomeSuiteError on
    missing/invalid JSON."""

def build_manifest(r_params: dict, wrapper: dict) -> dict:
    """Merge R's effective dada2 params + versions with wrapper-level run facts
    into the final manifest structure: {"tool": {...versions...},
    "dada2_params": {...effective...}, "run": {...wrapper facts...}}."""

def write_manifest(manifest: dict, out_dir: Path) -> Path:
    """Write dada2_denoise_manifest.json into out_dir; return its path."""
```

`wrapper` carries: `microsuite_version`, `backend` (`"dada2-r"`), `runtime`,
`image`, `mode`, `paired`, `threads` (resolved), `input_dir`, `output_table`,
`output_rep_seqs`, `output_stats`, `output_plot_dir`, `created_at` (UTC ISO),
`command` (the argv string). `build_manifest` splits R's flat dict into
`tool` (dada2_version, r_version) and `dada2_params` (everything else), and drops
`null`-valued dada2 keys so the manifest shows only the branch's real params.

### Component 3 — wiring in `denoise_dada2_r`

Pass `--params-out <output_stats.parent>/dada2_r_params.json` to the R command
(the output dir is rw-mounted under docker, so R writes it in-container and
Python reads it on the host). After a successful `_run` (in **both** the docker
and local branches, ungated by `validate`):

```python
r_params_path = output_stats.parent / "dada2_r_params.json"
try:
    r_params = dada2_manifest.read_r_params(r_params_path)
    manifest = dada2_manifest.build_manifest(r_params, wrapper_facts)
    dada2_manifest.write_manifest(manifest, output_stats.parent)
except MicrobiomeSuiteError as exc:
    warnings.warn(f"Could not write DADA2 provenance manifest: {exc}", stacklevel=2)
finally:
    r_params_path.unlink(missing_ok=True)
```

Provenance is best-effort: a missing/corrupt `dada2_r_params.json` emits
`warnings.warn` and does **not** fail an otherwise-successful denoise. The
intermediate R file is removed after merge, leaving one canonical manifest.

No new CLI flag is needed for the user — the manifest is always written. (The
`--params-out` flag is an internal wrapper→R detail.)

## Testing (offline)

- `read_r_params`: valid flat JSON → dict; missing file and malformed JSON →
  `MicrobiomeSuiteError`.
- `build_manifest`: given an R-params dict (paired) + wrapper dict → manifest has
  `tool.dada2_version`, `tool.r_version`, `dada2_params` with the effective
  values and no `null` keys, and `run` with runtime/image/mode/paths.
- `write_manifest`: writes `dada2_denoise_manifest.json`; JSON round-trips.
- `denoise_dada2_r` wiring (monkeypatched subprocess writing a fake
  `dada2_r_params.json` next to the stats path): a successful run leaves
  `dada2_denoise_manifest.json` with merged effective params + versions, and the
  intermediate `dada2_r_params.json` is removed; a run where the stub writes **no**
  params file still succeeds but emits a `warnings.warn` (assert via
  `pytest.warns`) and writes no manifest.

## Success criteria

1. A successful `dada2-r` run writes `dada2_denoise_manifest.json` beside the
   stats table with: effective dada2 params (post-default), `dada2`/R versions,
   and wrapper facts (runtime, image, microsuite version, mode, threads, paths,
   timestamp, command).
2. The effective params in the manifest equal what the R backend actually used
   (single source of truth via the hoisted `resolved_*` vars) — e.g.
   `min_overlap` shows `12` when unset, not absent.
3. Manifest generation is best-effort: a missing/unparseable R params file warns
   and does not fail the run.
4. Works identically for `--runtime local` and `docker`; the existing
   `run.json`/`events.jsonl` are unchanged.
5. The full offline suite stays green.

## Open questions / follow-ups (not blocking A)

- An opt-in real-dada2 e2e assertion that a live run's manifest carries a real
  `dada2_version` could extend `tests/integration/test_dada2_naming_contract_live.py`;
  left as a follow-up.
- Config-file provenance (if microsuite later grows a config file) would add a
  `config` block to `run`; out of scope now.
