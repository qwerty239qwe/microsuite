# Post-run output validation for `denoise` (P2) — Design

- **Date:** 2026-07-08
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** Codex's ERP120510 friction list, complaints #8 (weak output
  validation), #9 (Dropbox zero-block placeholder files completing "successfully"),
  and the testable half of #10 (wrapper/backend naming contract). **P2** of the
  four-part DADA2 roadmap (see [[dada2-improvement-roadmap]]); P1 (docker runtime)
  is merged.

## Scope

P2 adds post-run output validation to `microsuite denoise`: a reusable
integrity helper wired through `denoise._run` (all denoise backends), plus a
`dada2-r`-specific semantic check that the ASV table's sample columns match the
input sample set. Fail-by-default with a `--no-validate` escape hatch.

### Out of scope for P2
- Validation for methods other than `denoise` (the integrity helper is reusable
  but only wired into denoise here).
- DADA2 docker docs (#5/#6) — **P3**.
- The naming-contract end-to-end test (#10) — **P4** builds on the
  `_expected_sample_ids` helper introduced here.

## Verified context

- `denoise._run(command, failure_message, *, run_dir, timeout, backend,
  inputs=None, outputs=None, params=None)` wraps `run_command` with a
  `CommandLog`. `run_command` checks the exit code and records a results
  manifest but performs **no** output-file validation. The `outputs` dict
  (name→path string) is exactly what to validate.
- Every denoise backend path calls `_run(...)` with its `outputs` dict, so
  validating inside `_run` covers all backends.
- Post-#8, the dada2-r ASV table columns are clean sample IDs (`names(dada_out)
  <- samples`); P2 guards against regressions and cross-checks against inputs.

## Design

### Component 1 — reusable integrity helper `src/microsuite/runtime/validation.py`

```python
def validate_output_file(path: Path, *, allow_empty: bool = False) -> None:
    """Raise MicrobiomeSuiteError if the output is missing, empty, or (when the
    name ends '.gz') not a readable gzip."""
```
Behavior:
- missing → `MicrobiomeSuiteError(f"Expected output was not created: {path}")`.
- exists but 0 bytes (and not `allow_empty`) → `MicrobiomeSuiteError(f"Output is
  empty: {path}. This can mean an incomplete run or an unsynced cloud-storage
  placeholder.")` (catches #9).
- name ends `.gz`: open with `gzip.open(path, 'rb')` and read one 64 KiB chunk;
  on `gzip.BadGzipFile`/`EOFError`/`OSError` → `MicrobiomeSuiteError(f"Output is
  not a valid gzip file: {path}.")` (catches #8; truncated placeholder gz).

```python
def validate_outputs(outputs: dict[str, str], *, allow_empty: bool = False) -> None:
    for path_str in outputs.values():
        validate_output_file(Path(path_str), allow_empty=allow_empty)
```

### Component 2 — wire into `denoise._run`

Add `validate: bool = True` to `denoise._run`. After `run_command(...)` returns
(exit 0), if `validate`, call `validate_outputs(outputs or {})`. `denoise()`
gains `validate: bool = True` and passes it to every `_run` call. Result: any
denoise backend whose declared output is missing/empty/bad-gzip now fails loudly
instead of returning success.

### Component 3 — `dada2-r` semantic check (the #10 testable half)

New helpers in `denoise.py`:

```python
def _expected_sample_ids(input_dir: Path) -> set[str]:
    """Derive the expected sample IDs from the FASTQs in input_dir using the same
    PE/SE grouping the R backend uses (R1/R2, read1/read2, _1/_2, ..._001; single
    files are their own sample)."""

def _validate_dada2_asv_samples(output_table: Path, input_dir: Path) -> None:
    """Read the ASV table header (tab-separated); the columns are sample IDs.
    Raise MicrobiomeSuiteError if: no columns; any column is empty/duplicated;
    any column contains a FASTQ artifact ('.fastq', '.fq', '.filtered'); or the
    column set != _expected_sample_ids(input_dir) (name the unexpected/missing
    samples)."""
```

`_expected_sample_ids` uses the same regex family as the R script's
`read_suffix_pattern`/`stem_without_fastq_suffix` and the codex manifest builder
(patterns: `(?P<sample>.+?)[._-]R[12]([._-]001)?$`, `...[._-]read[12]...`,
`...[._-][12]...`; unmatched files are single-end samples named by their stem).
Deliberate independent re-derivation — the point is to catch a wrapper/backend
disagreement, not to share code with the R side. P4 formalizes and end-to-end
tests this contract.

In the `dada2-r` path of `denoise_dada2_r`, after `_run(...)` returns and when
`validate`, call `_validate_dada2_asv_samples(output_table, input_dir)`. This
runs identically for `--runtime local` and `--runtime docker` (P1) — the ASV
table is on the host either way. Thread `validate` into `denoise_dada2_r`.

### Component 4 — CLI

`method_features_cmd.py` `denoise` command gains:
```python
        no_validate: Annotated[
            bool, typer.Option("--no-validate", help="Skip post-run output validation.")
        ] = False,
```
and passes `validate=not no_validate` to `denoise(...)`. Default: validation on.

## Testing

Offline, in `tests/`:

1. **`validate_output_file`** (`tests/test_runtime_validation.py`): missing path
   → raises "not created"; empty file → raises "empty" (assert the
   placeholder/cloud wording); a `.gz` written from real gzip bytes → passes; a
   `.gz` whose bytes are not gzip → raises "not a valid gzip"; a non-empty `.tsv`
   → passes. `validate_outputs` raises on the first bad entry.
2. **denoise wiring** (`tests/test_denoise_cluster_methods.py`): monkeypatch
   `subprocess.run` to exit 0 but create **no** output files → `denoise(...)`
   raises `MicrobiomeSuiteError` (validation); the same call with
   `validate=False` (and the CLI `--no-validate`) does **not** raise.
3. **dada2-r ASV check**: `_expected_sample_ids` groups PE (`s_1/s_2`) and SE
   files correctly; `_validate_dada2_asv_samples` passes for a header whose
   columns equal the input samples, raises for a column containing
   `.filtered.fastq.gz`, and raises for a sample-set mismatch.

Existing denoise tests are updated where they monkeypatch a successful run but
don't create outputs — those must pass `validate=False` or write stub outputs, so
the new default validation doesn't break them. No real tool execution added.

## Success criteria

1. A denoise run whose declared output is missing / 0-byte / invalid-gzip fails
   with an actionable `MicrobiomeSuiteError` (default); `--no-validate` restores
   the prior no-check behavior.
2. A `dada2-r` run whose ASV columns contain FASTQ artifacts, or don't match the
   input sample set, fails with a clear message.
3. Validation is runtime-agnostic (local and docker) and lives in a reusable
   `runtime/validation.py` helper plus denoise-local semantic helpers.
4. The full offline suite stays green (existing tests adjusted for the new
   default).

## Open questions / follow-ups (not blocking P2)

- Whether the integrity helper should later be wired into other method wrappers
  (trim/qc/cluster) — reusable now, opt-in per method.
- `_expected_sample_ids` is P2's independent derivation; P4 turns the
  wrapper⇄backend naming agreement into an end-to-end test.
