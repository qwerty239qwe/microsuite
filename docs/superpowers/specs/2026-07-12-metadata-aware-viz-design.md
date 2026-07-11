# Metadata-aware taxonomy viz (Round-3 F) — Design

- **Date:** 2026-07-12
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** Codex's ERP120510 friction list, round-3 complaint **#18**
  (metadata-aware taxonomy visualizations are outside microsuite: top taxa by
  time, top taxa by phase, CLR by time, CLR by phase, subject/time-aware
  Bray-Curtis). Third and final round-3 sub-project (**F**); **D** and **E**
  merged. See [[dada2-improvement-roadmap]].

## Scope

F adds three parametrized metadata-aware `viz` commands. Because "by time" vs
"by phase" is just a different metadata column, each of #18's five views maps to
one of three commands driven by a `--group-by` / `--color-by` column read from
`adata.obs`:

| #18 view | Command |
|----------|---------|
| top taxa by time / by phase | `viz taxa-by-group --group-by <col>` |
| CLR by time / by phase | `viz clr-by-group --group-by <col>` |
| subject/time-aware Bray-Curtis | `viz braycurtis-ordination --color-by <col> [--subject <col>]` |

### Out of scope for F
- Statistical testing of group differences (PERMANOVA/ANCOM already exist as
  `diversity`/`diffab` commands); F is visualization only.
- New distance metrics or ordination methods (reuses `beta_diversity` +
  `pcoa`).
- The deferred cross-cutting taxonomy-NaN normalization (separate ticket).

## Verified context

- Metadata is in `adata.obs` (from `import tsv --metadata`). Read counts via
  `dense_counts` (samples×features).
- `abundance_native(adata, *, level, relative) -> pd.DataFrame`
  (`methods/abundance.py`) collapses to a taxonomy `level` (samples×taxa), values
  relative or counts; unclassified → `"Unclassified"`.
- `normalize_native(adata, *, method="clr", ...) -> ad.AnnData` returns an
  AnnData whose `X` is the CLR matrix (shape-preserving).
- `beta_diversity(adata, "bray-curtis") -> pd.DataFrame` (square, sample-indexed
  distance matrix); `pcoa(distance_matrix, *, dimensions) -> pd.DataFrame` with a
  `sample` column and `PC1..PCn` columns.
- Plotting stack: `viz/barplot.py` / `viz/assignment.py` use matplotlib with
  `matplotlib.use("Agg")` and `plt.savefig` + `plt.close(fig)`. F mirrors this.
- `viz` CLI group is `cli/viz_cmd.py` (`viz barplot`); F adds three commands
  there. Fatal → `MicrobiomeSuiteError` (`microsuite._errors`).

## Design

### Component 1 — `viz/metadata.py` (matplotlib Agg, mirrors existing viz)

Shared helpers + the plot functions. A small internal helper validates that a
named `obs` column exists (`MicrobiomeSuiteError` naming the column and listing
available columns) and that `--level` is a `var` column.

```python
def plot_taxa_by_group(adata, *, level: str, group_by: str, output: Path,
                       top_n: int = 15) -> None
def plot_clr_by_group(adata, *, level: str, group_by: str, output: Path,
                      top_n: int = 10, style: str = "boxplot") -> None
def plot_braycurtis_ordination(adata, *, color_by: str, output: Path,
                               subject: str | None = None,
                               style: str | None = None) -> None
```

**`plot_taxa_by_group`** — `frame = abundance_native(adata, level=level,
relative=True)` (samples×taxa); attach `group = adata.obs[group_by]`; mean
relative abundance per group (`frame.groupby(group).mean()`); keep the top-N taxa
by overall mean, fold the rest into `"Other"`; stacked bar, x = group, y =
mean relative abundance (0–1), legend = taxa.

**`plot_clr_by_group`** — collapse to counts:
`counts = abundance_native(adata, level=level, relative=False)` (samples×taxa);
build an AnnData `ad.AnnData(X=counts.values, obs=adata.obs.loc[counts.index],
var=pd.DataFrame(index=counts.columns))`; `clr = normalize_native(that,
method="clr").X`; select the top-N taxa by mean relative abundance; assemble a
long frame `(sample, taxon, clr, group)`. Render by `style`:
- `boxplot` (default): grouped boxplots — x = taxon, one box per group per taxon.
- `heatmap`: taxa (rows) × group (cols) mean CLR, `imshow` + colorbar (diverging
  cmap centered at 0).
- `violin`: per-taxon violins split by group (small multiples or grouped).
Unknown `style` → `MicrobiomeSuiteError` listing the three choices.

**`plot_braycurtis_ordination`** — `dist = beta_diversity(adata, "bray-curtis")`;
`coords = pcoa(dist, dimensions=2)` indexed by sample; join
`adata.obs[[color_by]]` (and `subject` if given). Resolve effective style:
`scatter` when `subject is None` (or `style == "scatter"`), else `trajectory`
default; explicit `style` overrides. Render:
- `scatter`: PC1/PC2 scatter colored by `color_by` (categorical legend or
  continuous colorbar — categorical for object/category dtype, continuous
  otherwise).
- `trajectory` (needs `subject`): scatter + for each subject a line connecting its
  points ordered by `color_by` (numeric sort if numeric, else category order),
  showing the subject's path through ordination space.
- `facet` (needs `subject`): one PC1/PC2 subplot per subject, points colored by
  `color_by`.
`trajectory`/`facet` without `--subject` → `MicrobiomeSuiteError`. Axis labels
carry the PCoA proportion-explained if `pcoa` exposes it, else `PC1`/`PC2`.

### Component 2 — CLI (`cli/viz_cmd.py`)

Three commands, each `read_h5ad(ensure_input(table))` then call Component 1 with
`prepare_output(output, force=force)`:

- `viz taxa-by-group`: `--table` (h5ad), `--level`, `--group-by`, `--output/-o`,
  `--top-n` (default 15, min 1), `--force`.
- `viz clr-by-group`: `--table`, `--level`, `--group-by`, `--output/-o`,
  `--top-n` (default 10, min 1), `--style` (default `boxplot`), `--force`.
- `viz braycurtis-ordination`: `--table`, `--color-by`, `--output/-o`,
  `--subject` (optional), `--style` (optional), `--force`.

### Data flow

`x.h5ad (obs has time/phase/subject) → viz taxa-by-group --group-by phase → PNG`;
`… → viz clr-by-group --group-by time --style heatmap → PNG`;
`… → viz braycurtis-ordination --color-by time --subject host → PNG`.

## Testing (offline)

Fixture: a small AnnData with taxonomy in `var` (a couple of ranks), read counts,
and `obs` columns `phase` (categorical, ≥2 groups), `time` (numeric), `subject`
(≥2 subjects, multiple timepoints each).

- `plot_taxa_by_group`: writes a non-empty PNG; a missing `group_by` column or a
  missing `level` → `MicrobiomeSuiteError`.
- `plot_clr_by_group`: each `style` in {boxplot, heatmap, violin} writes a
  non-empty PNG; unknown style → `MicrobiomeSuiteError`; CLR values used are
  finite.
- `plot_braycurtis_ordination`: `scatter` (no subject), `trajectory` (with
  subject), `facet` (with subject) each write a non-empty PNG; `trajectory`
  without `subject` → `MicrobiomeSuiteError`; a continuous `color_by` (time) and
  a categorical one (phase) both render.
- CLI smoke (`CliRunner`): all three commands exit 0 and write the PNG; a
  `--style` pass-through for `clr-by-group` and `braycurtis-ordination`.

## Success criteria

1. `viz taxa-by-group --group-by <col>` renders top-N taxa summarized per
   metadata group (covers "top taxa by time" and "by phase").
2. `viz clr-by-group --group-by <col> [--style …]` renders CLR of top taxa across
   groups in boxplot / heatmap / violin form (covers "CLR by time" and "by
   phase").
3. `viz braycurtis-ordination --color-by <col> [--subject <col>] [--style …]`
   renders a Bray-Curtis PCoA colored by group, with per-subject trajectories or
   facets when a subject column is given (covers "subject/time-aware
   Bray-Curtis").
4. Every command errors clearly when a named `obs` column, `--level`, or `--style`
   is invalid.
5. Full offline suite green and both CI gates pass
   (`ruff check .`, `ruff format --check .`).

## Open questions / follow-ups (not blocking F)

- Ordering of a categorical `--group-by` (e.g. phase) uses the data's category
  order / sorted unique; an explicit `--order` option could be added later.
- Exporting the underlying per-group tables (mean abundance, mean CLR, PCoA
  coords) as TSV alongside the plots could be a later addition; F focuses on the
  plots (#18 is a visualization gap).
- Continuous `--group-by` for the taxa/CLR bar/box plots is treated as discrete
  categories; binning is a later enhancement.
