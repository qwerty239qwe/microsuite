# Viz Hardening (G) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the four real-scale viz defects found on the synthetic 105-sample run: lexicographic group order, color collisions >10 categories, trajectory spaghetti, and the unbounded assignment reads-heatmap.

**Architecture:** Add small ordering/color helpers to `viz/metadata.py`; thread a `--group-order` option through `taxa-by-group`/`clr-by-group`, and `--order-by`/`--order` through `braycurtis-ordination`; bound the reads-heatmap in `viz/assignment.py`. All matplotlib Agg, existing patterns.

**Tech Stack:** Python 3.12, anndata, numpy, pandas, matplotlib (Agg), Typer, pytest.

## Global Constraints

- `_natural_group_order(values)`: values parseable as `float` sort numerically before non-numeric values (which sort as strings). `_resolve_group_order(present, group_order)`: when `group_order` is given, use those first (in that order, keeping only present ones), then append any remaining present groups in natural order.
- Colors: qualitative colors come from `tab20` (cycle if >20). `"Other"` is pinned to gray `"0.7"`.
- Trajectory ordering uses the `order_by` column (default = `color_by`), sorted by explicit `order` list if given else natural order — never by the color values when `order_by` differs.
- Reads-heatmap height capped at 16"; y-labels thinned when `> 60` samples.
- No behavior change to existing callers that omit the new options (defaults preserve current output except the intended fixes: natural order replaces lexicographic; explicit colors replace the default cycle).
- Fatal → `MicrobiomeSuiteError`. Both CI gates pass (`ruff check .`, `ruff format --check .`). `from __future__ import annotations` already present.

---

### Task 1: group ordering + distinct colors (`viz/metadata.py`, taxa-by-group & clr-by-group)

**Files:**
- Modify: `src/microsuite/viz/metadata.py`
- Modify: `src/microsuite/cli/viz_cmd.py`
- Test: `tests/test_viz_metadata.py`

**Interfaces:**
- Produces: `_natural_group_order(values) -> list[str]`; `_resolve_group_order(present, group_order) -> list[str]`; `_qualitative_colors(n) -> list`; `plot_taxa_by_group(..., group_order=None)`; `plot_clr_by_group(..., group_order=None)`; CLI `--group-order` on both.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_viz_metadata.py`)

```python
def test_natural_group_order() -> None:
    from microsuite.viz.metadata import _natural_group_order, _resolve_group_order

    assert _natural_group_order(["0", "14", "7", "35", "B"]) == ["0", "7", "14", "35", "B"]
    # explicit order honored, missing appended in natural order
    assert _resolve_group_order(["0", "7", "B", "14"], ["B", "0", "7"]) == ["B", "0", "7", "14"]
    # unknown entries in group_order are dropped
    assert _resolve_group_order(["0", "7"], ["B", "0", "7", "99"]) == ["0", "7"]
    # default is natural order
    assert _resolve_group_order(["7", "0", "14"], None) == ["0", "7", "14"]


def test_taxa_by_group_group_order_and_colors(tmp_path: Path) -> None:
    from microsuite.viz.metadata import _qualitative_colors, plot_taxa_by_group

    # 25 distinct colors requested -> no crash, all returned
    assert len(_qualitative_colors(25)) == 25
    out = tmp_path / "t.png"
    # group_order forces B first even though lexicographic would put it last
    plot_taxa_by_group(
        make_adata(), level="genus", group_by="time", output=out, top_n=2,
        group_order=["7", "0"],
    )
    assert out.exists() and out.stat().st_size > 0


def test_clr_by_group_group_order(tmp_path: Path) -> None:
    from microsuite.viz.metadata import plot_clr_by_group

    out = tmp_path / "c.png"
    plot_clr_by_group(
        make_adata(), level="genus", group_by="phase", output=out, style="heatmap",
        group_order=["post", "pre"],
    )
    assert out.exists() and out.stat().st_size > 0
```

(The `make_adata` fixture already exists in this file with `time`/`phase`/`subject` obs columns.)

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_viz_metadata.py -k "natural or group_order or colors" -v`
Expected: FAIL (`cannot import name '_natural_group_order'`).

- [ ] **Step 3: Add the helpers + thread `group_order`**

Add helpers after `_require_obs_column` in `viz/metadata.py`:

```python
def _natural_group_order(values) -> list[str]:
    unique = list(dict.fromkeys(str(v) for v in values))

    def key(value: str):
        try:
            return (0, float(value), "")
        except (TypeError, ValueError):
            return (1, 0.0, value)

    return sorted(unique, key=key)


def _resolve_group_order(present, group_order) -> list[str]:
    present_list = list(dict.fromkeys(str(p) for p in present))
    if group_order is None:
        return _natural_group_order(present_list)
    present_set = set(present_list)
    ordered = [str(g) for g in group_order if str(g) in present_set]
    ordered += [g for g in _natural_group_order(present_list) if g not in ordered]
    return ordered


def _qualitative_colors(n: int) -> list:
    cmap = plt.get_cmap("tab20")
    return [cmap(i % 20) for i in range(n)]
```

Modify `plot_taxa_by_group` — add `group_order: list[str] | None = None` to the signature, reorder the grouped rows, and colour the stacked bar explicitly:

```python
    group = _require_obs_column(adata, group_by)
    frame = abundance_native(adata, level=level, relative=True)
    grouped = frame.groupby(group.astype(str).to_numpy()).mean()
    grouped = grouped.reindex(_resolve_group_order(grouped.index, group_order))
    order = grouped.mean(axis=0).sort_values(ascending=False)
    top = list(order.head(top_n).index)
    plot_df = grouped[top].copy()
    other = grouped.drop(columns=top).sum(axis=1)
    if (other > 0).any():
        plot_df["Other"] = other
    taxa_cols = [c for c in plot_df.columns if c != "Other"]
    color_map = dict(zip(taxa_cols, _qualitative_colors(len(taxa_cols))))
    color_map["Other"] = "0.7"
    colors = [color_map[c] for c in plot_df.columns]
    width = max(6.0, len(plot_df.index) * 0.9)
    height = max(4.5, min(12.0, 2.5 + plot_df.shape[1] * 0.25))
    ax = plot_df.plot(kind="bar", stacked=True, figsize=(width, height), width=0.85, color=colors)
```
(keep the rest of the function unchanged.)

Modify `plot_clr_by_group` — add `group_order: list[str] | None = None`; replace the category derivation and the box/violin color source:

```python
    groups = group.to_numpy()
    categories = _resolve_group_order(pd.unique(groups), group_order)
```
and in the box/violin branch replace `cmap = plt.get_cmap("tab10")` + `color = cmap(gi % 10)` with:
```python
        colors = _qualitative_colors(len(categories))
        ...
            color = colors[gi]
```
(The heatmap branch already does `.reindex(categories)`, so the ordered `categories` fixes its row order too.)

- [ ] **Step 4: Add `--group-order` to the CLI**

In `cli/viz_cmd.py`, add to both `taxa_by_group` and `clr_by_group` commands:
```python
    group_order: Annotated[
        str | None, typer.Option("--group-order", help="Explicit group order, comma-separated.")
    ] = None,
```
and pass `group_order=(group_order.split(",") if group_order else None)` to the plot call.

- [ ] **Step 5: Run to verify pass**

Run: `uv run pytest tests/test_viz_metadata.py -v`
Expected: PASS (all, incl existing).

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/viz/metadata.py src/microsuite/cli/viz_cmd.py tests/test_viz_metadata.py
git commit -m "feat(viz): group ordering + distinct colors for taxa/clr-by-group"
```

---

### Task 2: trajectory order decoupled from color (`braycurtis-ordination`)

**Files:**
- Modify: `src/microsuite/viz/metadata.py` (`plot_braycurtis_ordination`)
- Modify: `src/microsuite/cli/viz_cmd.py`
- Test: `tests/test_viz_metadata.py` (append)

**Interfaces:**
- Consumes: `_natural_group_order` (Task 1).
- Produces: `plot_braycurtis_ordination(..., order_by=None, order=None)`; CLI `--order-by`, `--order`.

- [ ] **Step 1: Write the failing test** (append)

```python
def test_ordination_trajectory_order_by(tmp_path: Path) -> None:
    from microsuite.viz.metadata import plot_braycurtis_ordination

    out = tmp_path / "traj.png"
    # color by categorical phase, order path by numeric time -> must not crash / spaghetti-guard
    plot_braycurtis_ordination(
        make_adata(), color_by="phase", subject="subject", output=out,
        style="trajectory", order_by="time",
    )
    assert out.exists() and out.stat().st_size > 0

    out2 = tmp_path / "traj2.png"
    plot_braycurtis_ordination(
        make_adata(), color_by="phase", subject="subject", output=out2,
        style="trajectory", order_by="time", order=["7", "0"],
    )
    assert out2.exists() and out2.stat().st_size > 0
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_viz_metadata.py -k order_by -v`
Expected: FAIL (`plot_braycurtis_ordination() got an unexpected keyword argument 'order_by'`).

- [ ] **Step 3: Thread `order_by`/`order` into the trajectory**

Add `order_by: str | None = None` and `order: list[str] | None = None` to the signature. After the existing `color`/`subj` resolution, resolve the ordering column and a rank map:

```python
    order_column = _require_obs_column(adata, order_by) if order_by is not None else color
    order_vals = np.asarray([str(v) for v in order_column.to_numpy()])
    order_ranking = _resolve_group_order(pd.unique(order_vals), order)
    order_rank = {value: index for index, value in enumerate(order_ranking)}
```

Replace the trajectory block:
```python
    if effective == "trajectory":
        subj_vals = subj.to_numpy()
        for sub in pd.unique(subj_vals):
            idx = np.where(subj_vals == sub)[0]
            idx = sorted(idx, key=lambda i: order_rank.get(order_vals[i], len(order_rank)))
            ax.plot(x[idx], y[idx], color="gray", alpha=0.6, lw=1.0, zorder=0)
```
(`_resolve_group_order` from Task 1 gives natural order by default and honors an explicit `order` list — so a phase-colored, time-ordered path follows time, and `--order "B,0,7,…"` places `"B"` correctly.)

- [ ] **Step 4: Add `--order-by`/`--order` to the CLI**

In `cli/viz_cmd.py` `braycurtis_ordination`, add:
```python
    order_by: Annotated[
        str | None, typer.Option("--order-by", help="Obs column to order trajectories by (default: --color-by).")
    ] = None,
    order: Annotated[
        str | None, typer.Option("--order", help="Explicit order for --order-by values, comma-separated.")
    ] = None,
```
and pass `order_by=order_by, order=(order.split(",") if order else None)` to the plot call.

- [ ] **Step 5: Run to verify pass + no regression**

Run: `uv run pytest tests/test_viz_metadata.py -v`
Expected: PASS (all).

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/viz/metadata.py src/microsuite/cli/viz_cmd.py tests/test_viz_metadata.py
git commit -m "feat(viz): braycurtis-ordination --order-by/--order for trajectories"
```

---

### Task 3: bounded assignment reads-heatmap (`viz/assignment.py`)

**Files:**
- Modify: `src/microsuite/viz/assignment.py` (`plot_assigned_reads_by_rank`)
- Test: `tests/test_assignment_plots_cli.py` (append)

**Interfaces:**
- Produces: bounded-height `plot_assigned_reads_by_rank` (same signature).

- [ ] **Step 1: Write the failing test** (append to `tests/test_assignment_plots_cli.py`)

```python
def test_reads_heatmap_bounded_many_samples(tmp_path: Path) -> None:
    import anndata as ad
    import numpy as np
    import pandas as pd

    from microsuite.viz.assignment import plot_assigned_reads_by_rank

    # 90 samples -> old code would request ~46in height; new code caps it
    n = 90
    var = pd.DataFrame(
        {
            "kingdom": ["Bacteria"] * 3,
            "phylum": ["Firmicutes", "Bacteroidetes", ""],
            "genus": ["Lactobacillus", "", ""],
            "species": ["L. casei", "", ""],
        },
        index=["F1", "F2", "F3"],
    )
    rng = np.random.default_rng(0)
    X = rng.integers(1, 20, size=(n, 3)).astype(float)
    obs = pd.DataFrame(index=[f"S{i}" for i in range(n)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    out = tmp_path / "reads.png"
    plot_assigned_reads_by_rank(adata, out)
    assert out.exists() and out.stat().st_size > 0
    # figure height must be bounded (<= cap ~16in -> at 160 dpi ~2560px)
    from PIL import Image

    with Image.open(out) as im:
        assert im.height <= 2800
```

(If Pillow is unavailable in the env, assert on the matplotlib figure height instead: have the test import `matplotlib` and check the last figure — but the plan's target is the cap, so the PNG-pixel check is the acceptance. Confirm Pillow is a dep; microsuite already renders PNGs, and `PIL` ships with matplotlib's test extras — if the import fails, fall back to asserting only `out.stat().st_size > 0` and rely on the visual acceptance step.)

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_assignment_plots_cli.py -k bounded -v`
Expected: FAIL (height exceeds the cap / current unbounded figure).

- [ ] **Step 3: Bound the figure**

In `plot_assigned_reads_by_rank`, replace the `figsize`/`set_yticks`/`set_yticklabels` region:

```python
    samples = list(pd.Index(pivot.index).astype(str))
    n_samples = len(samples)
    height = min(0.22 * n_samples + 1.5, 16.0)
    fig, ax = plt.subplots(figsize=(max(6.0, len(ranks) * 1.1), height))
    im = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(ranks)))
    ax.set_xticklabels(ranks, rotation=45, ha="right")
    step = max(1, math.ceil(n_samples / 50)) if n_samples > 60 else 1
    tick_positions = list(range(0, n_samples, step))
    ax.set_yticks(tick_positions)
    ax.set_yticklabels([samples[i] for i in tick_positions], fontsize=(6 if n_samples > 60 else 8))
    ax.set_title("Assigned read fraction by rank")
    fig.colorbar(im, ax=ax, label="assigned read fraction")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)
```
Add `import math` at the top of `viz/assignment.py` if not present.

- [ ] **Step 4: Run to verify pass + full suite + lint**

Run: `uv run pytest tests/test_assignment_plots_cli.py -v`, then `uv run pytest -q`, then `uv run ruff check .` and `uv run ruff format --check .`.

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/viz/assignment.py tests/test_assignment_plots_cli.py
git commit -m "feat(viz): bound assignment reads-heatmap height for many samples"
```

---

## Self-Review

**Spec coverage:**
- Fix 1 group ordering (`--group-order`, natural default) → Task 1. ✓
- Fix 2 distinct colors (tab20 + gray Other) → Task 1. ✓
- Fix 3 trajectory `--order-by`/`--order` → Task 2. ✓
- Fix 4 bounded reads-heatmap → Task 3. ✓
- New unit tests + existing suite green + both gates → Tasks 1-3. ✓

**Placeholder scan:** none — full helper + edit code and complete tests. The Task-3 Pillow note gives an explicit fallback rather than leaving the assertion vague.

**Consistency:** `_natural_group_order`/`_resolve_group_order`/`_qualitative_colors` defined in Task 1 and reused by Task 2's `order_by`; `group_order`/`order_by`/`order` names align across `viz/metadata.py`, the CLI, and the tests; the acceptance step (re-run the synthetic 105-sample commands) validates all four fixes visually after merge.
