# DADA2 parameter-sensitivity sweep (Round-2 C) — Design

- **Date:** 2026-07-12
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** Codex's ERP120510 friction list, complaint **#10** (parameter-
  sensitivity reruns are still manual; microsuite should offer a small
  parameter-grid helper that compares retention, ASV count, sample depth,
  observed ASVs, chimera fraction, and correlation against a baseline run). Final
  DADA2 round-2 sub-project (**C**); B (QC) and A (provenance) merged. See
  [[dada2-improvement-roadmap]].

## Scope

C adds a DADA2 parameter-sweep helper: run the `dada2-r` backend across a small
grid of parameter sets, then emit a machine-readable `dada2_sweep_summary.tsv`
comparing each run's retention, ASV count, sample depth, observed ASVs, chimera
fraction, and similarity to a designated baseline run (shared-ASV abundance
correlation + per-sample metric correlation). The metrics/comparison logic is
pure and offline-testable; the orchestration is a thin wrapper over `denoise()`.

### Out of scope for C
- Automatic "best params" selection — C reports, the user decides.
- Non-dada2-r backends.
- Parallel execution of grid points (sequential; DADA2 is already multithreaded
  per run). A follow-up could parallelize.

## Verified context

- `denoise(*, backend="dada2-r", input_dir?, demux, output_table,
  output_rep_seqs, output_stats, mode/paired, max_ee_f, max_ee_r, trunc_len_f,
  trunc_len_r, trunc_q, min_overlap, runtime, dada2_image, threads, force, …)` is
  the single-run entry (`methods/denoise.py`).
- Per-run outputs: ASV table (features×samples TSV, first col = ASV id, written
  with `col.names=NA` leading empty header), rep-seqs FASTA (`>ASVn\n<seq>`), and
  the denoising-stats TSV. ASVs are matchable across runs by their **sequence**
  (rep-seqs).
- `dada2_qc.summarize_dada2_stats(stats_path)` (sub-project B) already returns
  overall `filtered_frac`/`merged_frac`/`nonchim_frac` + `bottleneck` + `paired`.
- `scipy.stats` is available (used elsewhere) for `pearsonr`/`spearmanr`.
- Fatal → `MicrobiomeSuiteError` (`microsuite._errors`); non-fatal →
  `warnings.warn`.

## Design

### Component 1 — grid builder (`methods/dada2_sweep.py`)

```python
@dataclass(frozen=True)
class GridPoint:
    name: str
    params: dict            # denoise() kwargs to override for this run
    is_baseline: bool

_SWEEP_AXES = ("max_ee_f", "max_ee_r", "trunc_len_f", "trunc_len_r",
               "trunc_q", "min_overlap")  # single-mode: max_ee, trunc_len also allowed

def build_grid(*, config: Path | None, axes: dict[str, list]) -> list[GridPoint]:
    """Config XOR axes. Exactly one GridPoint has is_baseline=True."""
```

- **Config (JSON)** — a list of objects `{"name": str, "params": {...}, "baseline":
  bool}`; exactly one `baseline: true` (else `MicrobiomeSuiteError`); names must be
  unique and filesystem-safe.
- **CLI axes** — `axes` maps a swept param → its list of values; the Cartesian
  product forms the grid; **grid point 0 (first value of every axis) is the
  baseline**. Each point's `name` is auto-derived from its non-baseline-differing
  params (e.g. `maxEEf3_truncLenf220`), with the baseline named `baseline`.
- Providing both config and axes, or neither, → `MicrobiomeSuiteError`.

### Component 2 — per-run metrics + baseline comparison (pure)

```python
def run_metrics(table_path: Path, stats_path: Path) -> dict:
    """n_asvs, filtered_frac, merged_frac, nonchim_frac, chimera_frac,
    mean_sample_depth, median_sample_depth, mean_observed_asvs."""

def compare_to_baseline(*, table_path, rep_seqs_path, baseline_table_path,
                        baseline_rep_seqs_path) -> dict:
    """shared_asv_count, frac_baseline_reads_shared, abundance_pearson,
    abundance_spearman, depth_pearson, observed_asv_pearson."""
```

- `run_metrics`: read the ASV table (features×samples) and stats. `n_asvs` = row
  count. Per-sample depth = column sums → `mean`/`median`. `mean_observed_asvs` =
  mean per-sample count of nonzero ASVs. Retention fracs from
  `summarize_dada2_stats`. `chimera_frac` = `(pre − nonchim) / pre` from the stats
  totals, where `pre` = `merged` (paired) or `denoised` (single); `0.0` when
  `pre == 0`.
- `compare_to_baseline`: parse both rep-seqs FASTAs to `{asv_id: sequence}`; map
  each run's ASV rows to sequences; `shared` = sequences present in both. Report
  `shared_asv_count`; `frac_baseline_reads_shared` = baseline reads in shared
  ASVs / total baseline reads; `abundance_pearson`/`abundance_spearman` on the
  per-sequence total-abundance vectors over `shared` (run vs baseline);
  `depth_pearson` and `observed_asv_pearson` on the per-sample vectors (samples
  are identical across runs; align by sample id). Correlations are `float("nan")`
  when undefined (< 2 shared points or zero variance).

### Component 3 — summary assembly (pure)

```python
def summarize_sweep(runs: list[SweepRun]) -> pd.DataFrame:
    """One row per run: name, is_baseline, the run's swept params, run_metrics,
    and compare_to_baseline vs the baseline run (baseline row compares to itself
    → shared=all, correlations 1.0). Ordered baseline-first."""

def write_sweep_summary(summary: pd.DataFrame, out_path: Path) -> Path
```

`SweepRun` bundles a `GridPoint` with its output paths (`table`, `rep_seqs`,
`stats`) and a `status` (`"ok"` / `"failed"`). Failed runs contribute a row with
metric/comparison columns as `NaN` and `status="failed"`.

### Component 4 — orchestration (`methods/dada2_sweep.py`)

```python
def run_dada2_sweep(*, input_dir: Path, mode: str, output_dir: Path,
                    grid: list[GridPoint], runtime="local", dada2_image=None,
                    threads=1, force=False, timeout=None) -> Path:
    """Run denoise(backend='dada2-r') for each GridPoint into
    output_dir/<name>/{table.tsv,rep_seqs.fasta,stats.tsv}; assemble and write
    output_dir/dada2_sweep_summary.tsv; return its path."""
```

- Runs the **baseline point first**; if the baseline run fails →
  `MicrobiomeSuiteError` (no baseline to compare against). Each non-baseline
  point runs next; a failed point emits `warnings.warn` and is recorded
  `status="failed"` (the sweep continues — an over-aggressive param set that
  yields no ASVs must not sink the whole sweep).
- Calls `denoise(backend="dada2-r", demux=input_dir, mode=mode, output_table=…,
  output_rep_seqs=…, output_stats=…, runtime=…, dada2_image=…, threads=…,
  force=…, **point.params)` (the reads dir is `denoise`'s `demux` argument).
  Reuses the existing per-run validation (`denoise(validate=True)`).

### Component 5 — CLI (`cli/method_features_cmd.py`)

`denoise-sweep` (beside `denoise`): `--input-dir`, `--mode`/`--paired`,
`--output-dir`, and the grid source — `--grid-config PATH` **or** the axis
options (`--max-ee-f`, `--max-ee-r`, `--trunc-len-f`, `--trunc-len-r`,
`--trunc-q`, `--min-overlap`; each a comma-separated list) — plus `--runtime`,
`--dada2-image`, `--threads`, `--force`. Providing both a config and axes, or
neither, errors. Calls `build_grid` then `run_dada2_sweep`.

### Data flow

`reads/ + grid (config or axes) → run_dada2_sweep → per-point denoise runs →
dada2_sweep_summary.tsv` (baseline-first, one row per param set, retention/ASV/
depth/observed/chimera + baseline-similarity columns).

## Testing (offline)

- `build_grid`: a JSON config → the expected `GridPoint`s with one baseline;
  duplicate/zero/multiple baselines → `MicrobiomeSuiteError`; CLI axes → the
  Cartesian product with point 0 as baseline; both/neither source → error.
- `run_metrics`: a fixture ASV table + stats TSV (paired and single) → exact
  `n_asvs`, depth mean/median, `mean_observed_asvs`, retention fracs, and
  `chimera_frac` (with a known merged/nonchim so the value is checkable);
  `pre == 0` → `chimera_frac == 0.0`.
- `compare_to_baseline`: fixtures where run and baseline share a known subset of
  sequences → exact `shared_asv_count`, `frac_baseline_reads_shared`, and
  correlations (use vectors with a known Pearson, e.g. identical → 1.0); disjoint
  sequences → `shared_asv_count == 0` and `nan` correlations.
- `summarize_sweep`: assembles a baseline-first DataFrame with the right columns;
  a `failed` run yields a NaN row.
- `run_dada2_sweep` orchestration: monkeypatch `denoise` to write fake per-point
  outputs (and to raise for one designated point) → `dada2_sweep_summary.tsv`
  exists with a row per point, the raising point marked `failed`, and a baseline
  failure raising `MicrobiomeSuiteError`.
- CLI smoke (`CliRunner`) with a monkeypatched `run_dada2_sweep`/`denoise`:
  `denoise-sweep` with `--grid-config` and with axis options both exit 0; both /
  neither source → non-zero exit.

## Success criteria

1. `microsuite denoise-sweep --input-dir reads --mode paired --output-dir out
   --grid-config grid.json` (or with `--max-ee-f 2,3,5 …`) runs each param set and
   writes `dada2_sweep_summary.tsv` with one baseline-first row per run.
2. Each row carries retention (filtered/merged/nonchim fracs), `n_asvs`, sample
   depth (mean/median), `mean_observed_asvs`, `chimera_frac`, and baseline
   similarity (`shared_asv_count`, `frac_baseline_reads_shared`, abundance
   Pearson/Spearman, per-sample depth/observed Pearson).
3. Exactly one baseline; CLI-axes baseline is grid point 0; a config must flag
   exactly one baseline. A non-baseline point that fails is recorded `failed` and
   does not abort the sweep; a baseline failure aborts with a clear error.
4. The metrics/comparison/grid logic is fully offline-tested; orchestration is
   tested with a stubbed `denoise`. Full offline suite green and both CI gates
   pass (`ruff check .`, `ruff format --check .`).

## Open questions / follow-ups (not blocking C)

- An opt-in real-DADA2 end-to-end sweep test (like
  `tests/integration/test_dada2_naming_contract_live.py`) would exercise the
  actual R runs; left as a follow-up (a maintainer runs it with real dada2).
- Parallel grid execution and a "recommended params" heuristic are later
  enhancements.
- Continuous vs discrete handling of axis values is out of scope — axes are
  explicit value lists.
