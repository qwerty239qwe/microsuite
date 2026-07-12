# Viz hardening for real-scale metadata plots (G) — Design

- **Date:** 2026-07-12
- **Status:** Proposed — pending user review
- **Origin:** Validating sub-project F's `viz` commands on a synthetic
  ERP120510-scale dataset (105 samples, 15 subjects × 7 timepoints, phases)
  surfaced four "weird plot" defects that would appear on real data. G fixes them
  so the standardized plots are usable before the oral-pipeline scripts adopt
  them. See [[metadata-aware-viz-design]] and the scripts-refactor spec in the
  oral repo.

## Evidence (from the synthetic 105-sample run)

1. `tax_assignment_plots` reads-heatmap rendered as a **1232×8560 px, ~53-inch
   tall** strip — height scales unbounded with sample count
   (`figsize height = len(samples)*0.5 + 1`).
2. `viz taxa-by-group` with 12 taxa + "Other" **reused colors** (Streptococcus ≡
   "Other" green, etc.) — the default cycle has only 10 colors.
3. Time axis rendered **lexicographically** (`0,14,21,28,35,7,B`) instead of the
   biological order (`B,0,7,14,21,28,35`) in `taxa-by-group`/`clr-by-group`.
4. `braycurtis-ordination --style trajectory` produced **spaghetti**: it orders
   each subject's path by the `--color-by` column, so a categorical `phase_code`
   (or a `time_code` containing the non-numeric `"B"`) yields arbitrary order.

## Design

### Fix 1 — explicit / natural group ordering (`viz/metadata.py`)

`plot_taxa_by_group` and `plot_clr_by_group` gain `group_order: list[str] | None
= None`.
- When given, the group axis uses exactly that order (data groups not listed are
  appended after, in natural order; listed-but-absent groups are skipped).
- When omitted, groups are ordered by a shared **natural/numeric-aware sort**
  (`_natural_group_order`): values parseable as numbers sort numerically, the
  rest sort as strings after the numbers. This alone fixes `0,7,14,21,28,35`; the
  `"B"`/`-7` case is handled by passing `--group-order "B,0,7,14,21,28,35"`.
- Applies to the taxa-by-group x-axis and to the clr-by-group group axis (heatmap
  columns/rows and boxplot/violin group order).
- CLI: both commands gain `--group-order TEXT` (comma-separated).

### Fix 2 — distinct colors for many categories (`viz/metadata.py`)

`plot_taxa_by_group` stops relying on the default 10-color cycle. Build an
explicit color list from a 20-entry qualitative map (`tab20`), cycling if a
plot has >20 taxa; `"Other"` is pinned to a neutral gray. `plot_clr_by_group`'s
per-group colors move to `tab20` too (harmless; groups are few). This removes the
collision at 12–13 taxa.

### Fix 3 — trajectory ordering decoupled from color (`viz/metadata.py`)

`plot_braycurtis_ordination` gains `order_by: str | None = None` and `order:
list[str] | None = None`.
- `order_by` (obs column) is the column whose values order each subject's
  trajectory; default = `color_by` (preserves current behavior).
- Each subject's points are sorted by `order_by` using the explicit `order` list
  if given, else the natural/numeric-aware sort. This fixes phase-colored,
  time-ordered trajectories and the `"B"` case (`--order "B,0,7,14,21,28,35"`).
- CLI: `braycurtis-ordination` gains `--order-by TEXT` and `--order TEXT`.

### Fix 4 — bounded reads-heatmap (`viz/assignment.py`)

`plot_assigned_reads_by_rank` caps the figure so it stays usable at large sample
counts:
- height = `min(0.22 * n_samples + 1.5, 16.0)` inches (cap ~16");
- when `n_samples > 60`, show a thinned subset of y tick labels (every
  `ceil(n_samples / 50)`) at a small font, so the axis is legible rather than
  105 overlapping labels;
- width unchanged (ranks ≤ 7).

## Testing

- **Unit (offline, in-repo):** extend `tests/test_viz_metadata.py` /
  `tests/test_assignment_plots_cli.py`:
  - `_natural_group_order` orders numeric-like values numerically and mixes
    sensibly; `--group-order` reorders and appends/omits correctly.
  - `taxa-by-group` with >20 taxa still assigns colors (no crash) and "Other" is
    distinct.
  - `braycurtis-ordination` with `order_by` ≠ `color_by` runs; `order` list
    respected; a fixture where color-by is categorical and order-by is numeric
    yields a monotone-by-order trajectory (assert point order used).
  - `plot_assigned_reads_by_rank` on a many-sample fixture writes a PNG and the
    figure height is ≤ the cap (assert via the returned/inspected figure size or
    that it doesn't raise / stays bounded).
- **Acceptance (empirical):** re-run all commands on the synthetic 105-sample
  dataset and visually confirm: reads-heatmap is a normal-proportioned figure;
  taxa colors are distinct; `--group-order "B,0,7,14,21,28,35"` gives biological
  time order; `braycurtis-ordination --color-by phase_code --order-by time_code
  --order "B,0,7,..." --style trajectory` follows time, not spaghetti.

## Success criteria

1. `tax_assignment_plots` reads-heatmap is a bounded, legible figure at 105
   samples.
2. `taxa-by-group`/`clr-by-group` render distinct colors for ≥13 categories and
   honor `--group-order` (and natural order by default).
3. `braycurtis-ordination` trajectories follow `--order-by`/`--order` (time),
   independent of `--color-by` (phase).
4. Existing viz tests stay green; new unit tests cover the four fixes; full
   offline suite + both CI gates pass.

## Out of scope
- Alpha-diversity plots, Bray-Curtis distance heatmap, baseline-dissimilarity
  boxplots (the oral scripts keep these; no microsuite command).
- Per-sample CLR heatmap with `B→-7` custom styling (the standardized per-group
  clr-by-group heatmap replaces it under the aggressive decision).
