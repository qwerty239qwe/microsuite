# DADA2 QC warnings & summaries (Round-2 B) — Design

- **Date:** 2026-07-09
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** Codex's ERP120510 round-2 complaints #13 (no retention warning),
  #15 (QC plots need machine-readable summaries), #16 (no overlap check for
  long amplicons); absorbs #12 (default maxEE too strict — surfaced via a
  warning rather than silently changing defaults). Sub-project **B** of the
  round-2 roadmap (see [[dada2-improvement-roadmap]]); A (#11 provenance) and C
  (#14 sweep) follow.

## Scope

B adds DADA2 run-quality feedback to the `dada2-r` backend: a post-run retention
summary (TSV+JSON) with warnings on low retention, a pre-run overlap check for
paired long amplicons, and a `--strict-qc` flag that turns the warnings fatal.
Pure Python (no R-script change); works for both `--runtime local` and `docker`.

### Out of scope for B
- Parameter provenance / resolved-config manifest (#11) — **A**.
- Parameter-sensitivity sweep (#14) — **C**.
- Changing DADA2's default `maxEE`/`truncLen` (a warning surfaces #12, not a
  silent default change).
- QC for non-dada2-r backends.

## Verified context

- The denoising-stats TSV is written by R `write.table(track, ..., col.names=NA)`
  → header `<empty>\tinput\tfiltered\t...`; rows are samples. Columns:
  paired = `input, filtered, denoised_f, denoised_r, merged, nonchim`;
  single = `input, filtered, denoised, nonchim`.
- microsuite has **no non-fatal notice channel** (methods raise
  `MicrobiomeSuiteError` for fatal only; CLI uses `typer.echo`). B introduces
  warnings via Python `warnings.warn`, which surfaces in CLI stderr and the SDK.
- Read length is not captured; the overlap check derives it from `truncLen` (if
  set) or by peeking the first record of an input FASTQ.
- P1/P2 established the post-run hook point: `denoise_dada2_r` runs `_run`, then
  (gated by `validate`) P2's `_validate_dada2_asv_samples`. B adds its summary
  there and its overlap check pre-flight.

## Design

### Component 1 — `src/microsuite/methods/dada2_qc.py`

```python
_RETENTION_THRESHOLDS = {"filtered": 0.5, "merged": 0.5, "nonchim": 0.4}

def summarize_dada2_stats(stats_path: Path) -> dict:
    """Parse the denoising-stats TSV and return per-sample + overall retention.

    Returns {"per_sample": {sample: {"input": int, "filtered_frac": float,
    "merged_frac": float|None, "nonchim_frac": float}, ...},
    "overall": {"input": int, "filtered_frac": float, "merged_frac": float|None,
    "nonchim_frac": float}, "bottleneck": "<step with the largest fractional
    drop>", "paired": bool}."""

def write_qc_summary(summary: dict, out_dir: Path) -> tuple[Path, Path]:
    """Write dada2_qc_summary.json and dada2_qc_summary.tsv into out_dir;
    return their paths."""

def retention_warnings(summary: dict) -> list[str]:
    """Return a message per overall fraction below its threshold, each naming the
    fraction, its value, and suggested checks (maxEE, truncLen, overlap)."""

def check_overlap(*, trunc_len_f: int, trunc_len_r: int, read_len_f: int,
                  read_len_r: int, amplicon_length: int, min_overlap: int) -> str | None:
    """retained_x = trunc_len_x if trunc_len_x > 0 else read_len_x; overlap =
    retained_f + retained_r - amplicon_length; return a warning string if overlap
    < min_overlap, else None."""

def first_read_length(fastq: Path) -> int:
    """Length of the first sequence record in a (optionally gzipped) FASTQ."""
```

Retention fractions use `input` as the denominator; `merged_frac` is `None` in
single mode. `bottleneck` = the transition (`input→filtered`, `filtered→merged`,
`merged→nonchim`, etc.) with the largest fractional loss overall.

### Component 2 — retention summary + warnings in `denoise_dada2_r` (#13, #15)

After a successful `_run(...)` and when `validate` (consistent with P2 — so
`--no-validate` skips QC too): `summary = summarize_dada2_stats(output_stats)`;
`write_qc_summary(summary, output_stats.parent)`; for each message in
`retention_warnings(summary)`: if `strict_qc` raise `MicrobiomeSuiteError`, else
`warnings.warn(message)`. Runs identically for local and docker (the stats TSV is
on the host both ways).

### Component 3 — pre-run overlap check (#16)

New optional `amplicon_length: int | None`. When set and the run is paired,
before `_run`: derive `read_len_f/r` (peek the first R1/R2 FASTQ if the matching
`trunc_len` is 0), call `check_overlap(...)` with the effective `min_overlap`
(the `--min-overlap` value or DADA2's default 12); if it returns a message,
`strict_qc` → raise, else `warnings.warn`. Skipped when `amplicon_length` is
None or the run is single-end.

### Component 4 — CLI + threading

`denoise` CLI gains `--amplicon-length INT` (optional) and `--strict-qc`
(default off). `denoise()` and `denoise_dada2_r()` gain `amplicon_length: int |
None = None` and `strict_qc: bool = False`, threaded to the checks.

## Testing (offline)

- `summarize_dada2_stats`: a paired fixture stats TSV → correct per-sample and
  overall fractions, correct `bottleneck`; a single fixture → `merged_frac` None.
- `write_qc_summary`: writes both files; JSON round-trips the summary.
- `retention_warnings`: a low-retention summary yields the expected messages; a
  healthy one yields `[]`.
- `check_overlap`: insufficient overlap → message; sufficient → None; respects
  `trunc_len` vs read length.
- `denoise_dada2_r` wiring (monkeypatched subprocess writing a fixture stats
  TSV): low retention emits a `warnings.warn` (assert via `pytest.warns`) and,
  with `strict_qc=True`, raises; `--no-validate` skips QC; a paired run with a
  too-short `amplicon_length` warns/raises pre-run.

## Success criteria

1. A successful `dada2-r` run writes `dada2_qc_summary.{json,tsv}` beside the
   stats table with per-sample + overall retention and the bottleneck step.
2. Overall retention below the documented thresholds emits actionable warnings;
   `--strict-qc` makes them fatal.
3. With `--amplicon-length`, a paired run warns pre-run when `truncLen` leaves
   insufficient overlap; `--strict-qc` makes it fatal; absent the flag it is
   skipped.
4. All QC is gated by `validate` (so `--no-validate` skips it) and runs the same
   for local and docker.
5. The full offline suite stays green.

## Open questions / follow-ups (not blocking B)

- Retention thresholds are fixed documented heuristics; making them configurable
  is a later addition if needed.
- Error-rate-plot facts beyond retention (from `plotErrors`) are not summarized
  here (they'd need R-side extraction); the stats-derived summary covers the
  bottleneck facts codex asked for.
