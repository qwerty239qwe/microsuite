# Taxonomy assignment QC (Round-3 E) — Design

- **Date:** 2026-07-12
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** Codex's ERP120510 friction list, round-3 complaints **#16**
  (assignment summary outputs are easy to plot incorrectly — provide canonical
  assignment plots) and **#17** (`taxonomy_unassigned_summary.tsv` should be a
  native taxonomy QC output, not only a project-wrapper artifact). Second
  sub-project (**E**) of round-3; **D** merged, **F** (metadata-aware viz #18)
  follows. See [[dada2-improvement-roadmap]].

## Scope

E makes taxonomy-assignment QC a first-class, hard-to-misplot output: a native
per-`(sample, rank)` `taxonomy_unassigned_summary.tsv` (#17) and three canonical
assignment plots (#16) — assigned-vs-unassigned ASV% by rank, a per-sample
assigned-read-fraction heatmap by rank, and a deepest-assigned-rank distribution
(each ASV counted once). Both derive from one computation over `adata.var`
(assignment) and `adata.X` (read counts).

### Out of scope for E
- Metadata-aware taxonomy visualizations (top taxa by time/phase, CLR by phase,
  subject-aware Bray-Curtis) — that is **F** (#18).
- Producing taxonomy itself (`tax_classify` already does); E consumes an `.h5ad`
  that already carries taxonomy in `var`.
- Changing the taxonomy parsing / `add_taxonomy_levels` contract.

## Verified context

- Taxonomy is feature-level in `adata.var`: `add_taxonomy_levels`
  (`io/taxonomy.py`) fills columns `LEVELS = [kingdom, phylum, class, order,
  family, genus, species]`; an **unassigned** rank is the empty string `""`.
- Read counts: `dense_counts(adata)` (`diversity/_matrix.py`) is samples×features
  (`obs_names`×`var_names`); reads-per-feature = column sum, per-sample reads =
  the sample's row.
- No native assignment summary exists in `src/` today (it is a SILVA-wrapper
  artifact — the exact gap #17 names).
- Plotting stack: `viz/barplot.py` uses matplotlib with `matplotlib.use("Agg")`
  and `plt.savefig`; E mirrors this in a new `viz/assignment.py`.
- Existing taxonomy commands are flat `tax_*` registered in
  `cli/method_taxonomy_cmd.py` (`tax_classify`, `tax_barplot`, `tax_collapse`,
  …); E adds two more there. Fatal → `MicrobiomeSuiteError`
  (`microsuite._errors`).

## Design

### Component 1 — `methods/assignment_qc.py` (pure, offline-testable)

```python
def summarize_assignment(adata: ad.AnnData) -> pd.DataFrame:
    """Long-format taxonomy-assignment summary. One row per (sample, rank) plus
    pooled rows with sample == '_all_samples_'. A feature counts toward a sample
    when it has >0 reads there; 'assigned at rank' means var[rank] != ''.
    Columns: sample, rank, assigned_features, unassigned_features,
    assigned_reads, unassigned_reads, assigned_feature_frac, assigned_read_frac.
    Ranks are LEVELS present in var, in LEVELS order."""

def write_assignment_summary(summary: pd.DataFrame, out_path: Path) -> Path:
    """Write taxonomy_unassigned_summary.tsv (the long-format summary)."""

def deepest_rank_distribution(adata: ad.AnnData) -> pd.Series:
    """For each feature, the deepest rank with a non-empty value (each ASV counted
    once); return counts indexed by rank in LEVELS order, with an 'Unassigned'
    entry for features assigned at no rank."""
```

- Ranks used = `[r for r in LEVELS if r in adata.var.columns]`; raise
  `MicrobiomeSuiteError` if none are present ("table has no taxonomy rank
  columns; run tax_classify or import with --taxonomy first").
- `assigned_mask[rank]` = `adata.var[rank].astype(str) != ""` (per feature).
- Per sample `s`: `present = counts[s] > 0`; for each rank —
  `assigned_features = (present & assigned_mask[rank]).sum()`,
  `unassigned_features = (present & ~assigned_mask[rank]).sum()`,
  `assigned_reads = counts[s][assigned_mask[rank]].sum()`,
  `unassigned_reads = counts[s][~assigned_mask[rank]].sum()`. Fracs use
  `assigned / (assigned + unassigned)`, `0.0` when the denominator is 0.
- Pooled `_all_samples_` rows: feature counts over `present_any = counts.sum(0) >
  0`; reads over the full column sums. Same columns.
- Deepest rank per feature: the last rank in LEVELS order where the value is
  non-empty; features with all-empty ranks → `Unassigned`. The returned Series is
  reindexed to `LEVELS (present) + ["Unassigned"]`, missing entries `0`.

### Component 2 — `viz/assignment.py` (matplotlib Agg, mirrors `viz/barplot.py`)

```python
def plot_assigned_asv_by_rank(adata: ad.AnnData, output: Path) -> None
def plot_assigned_reads_by_rank(adata: ad.AnnData, output: Path) -> None
def plot_deepest_rank(adata: ad.AnnData, output: Path) -> None
```

- `plot_assigned_asv_by_rank` — a stacked bar, x = rank, two series (assigned /
  unassigned **feature** fraction, summing to 1). Feature assignment is
  sample-independent, so this uses the pooled `_all_samples_` feature counts.
- `plot_assigned_reads_by_rank` — a **per-sample heatmap**: rows = samples,
  columns = ranks, cell = `assigned_read_frac` for that (sample, rank), from the
  per-sample summary rows (excludes `_all_samples_`). `imshow`/`pcolormesh` with a
  colorbar (0–1); rank labels on x, sample labels on y.
- `plot_deepest_rank` — a bar of `deepest_rank_distribution` (x = rank incl
  `Unassigned`, y = feature count).
- Each function computes via Component 1 (`summarize_assignment` /
  `deepest_rank_distribution`), so the plots and the TSV are guaranteed
  consistent. Each sets `matplotlib.use("Agg")` at import (as `barplot.py` does)
  and `plt.savefig(output, ...)`; `plt.close(fig)` after.

### Component 3 — CLI (`cli/method_taxonomy_cmd.py`)

Two flat commands (matching the existing `tax_*` pattern), each reading the h5ad
via `read_h5ad(ensure_input(table))`:

- `tax_assignment_summary`: `--table PATH` (h5ad, required), `--output/-o PATH`
  (required), `--force`. Calls `summarize_assignment` + `write_assignment_summary`.
- `tax_assignment_plots`: `--table PATH` (h5ad, required), `--output-dir PATH`
  (required; created if absent), `--force`. Writes
  `assigned_asv_by_rank.png`, `assigned_reads_by_rank.png`, `deepest_rank.png`
  into the directory via the three Component-2 functions (each through
  `prepare_output(dir / name, force=force)`).

### Data flow

`x.h5ad (var has taxonomy) → tax_assignment_summary → taxonomy_unassigned_summary.tsv`.
`x.h5ad → tax_assignment_plots → {assigned_asv_by_rank, assigned_reads_by_rank,
deepest_rank}.png`.

## Testing (offline)

- `summarize_assignment`: a small fixture AnnData (e.g. 3 features with known
  rank assignment — one assigned to species, one to genus only, one unassigned —
  across 2 samples with a known presence/read pattern) → assert exact
  `assigned_features`/`unassigned_features`/`assigned_reads`/`unassigned_reads`
  and fracs for a couple of `(sample, rank)` cells and the `_all_samples_` rows;
  assert monotonicity (assigned_features non-increasing from kingdom→species).
- `deepest_rank_distribution`: the same fixture → each ASV counted once, the
  genus-only feature lands in `genus`, the unassigned in `Unassigned`, sum ==
  n_features.
- Error path: an AnnData with no LEVELS columns → `MicrobiomeSuiteError`.
- `write_assignment_summary`: writes the TSV; round-trips via `read_csv`.
- Plots: each of the three functions writes a non-empty PNG (assert file exists
  and size > 0) on the fixture; the reads heatmap handles a single-sample table.
- CLI smoke (`CliRunner`): `tax_assignment_summary` writes the TSV;
  `tax_assignment_plots` writes all three PNGs into the dir; both exit 0.

## Success criteria

1. `microsuite tax_assignment_summary --table x.h5ad -o summary.tsv` writes a
   native long-format `taxonomy_unassigned_summary.tsv` with per-`(sample, rank)`
   assigned/unassigned feature and read counts (+ fracs) and pooled
   `_all_samples_` rows.
2. `microsuite tax_assignment_plots --table x.h5ad --output-dir d` writes the
   three canonical PNGs (assigned ASV% by rank; per-sample assigned-read-fraction
   heatmap by rank; deepest-assigned-rank distribution, each ASV once).
3. Plots and the TSV are computed from the same functions, so they never
   disagree; a table without taxonomy ranks errors clearly.
4. Full offline suite green and both CI gates pass
   (`ruff check .`, `ruff format --check .`).

## Open questions / follow-ups (not blocking E)

- A pooled (non-heatmap) reads-by-rank bar could be added later for a quick
  single-figure view; the per-sample heatmap plus the TSV cover the QC need now.
- Assignment confidence/score summaries (if a classifier emits per-feature
  confidences into `var`) could extend the summary later.
- Metadata-grouped assignment views (assignment by phase/time) belong to **F**.
