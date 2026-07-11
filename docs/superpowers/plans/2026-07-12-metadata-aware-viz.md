# Metadata-aware Taxonomy Viz (Round-3 F) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Three metadata-aware `viz` commands — `taxa-by-group`, `clr-by-group`, `braycurtis-ordination` — that render #18's five views via a `--group-by`/`--color-by` obs column.

**Architecture:** A `viz/metadata.py` (matplotlib Agg, mirrors `viz/barplot.py`) holds the plot functions and a shared obs-column validator; they reuse `abundance_native`, `normalize_native`, `beta_diversity`, `pcoa`. Three commands in `cli/viz_cmd.py` wire them to `.h5ad` input.

**Tech Stack:** Python 3.12, anndata, numpy, pandas, matplotlib (Agg), Typer, pytest (`CliRunner`).

## Global Constraints

- Metadata is in `adata.obs`; a named column that is absent → `MicrobiomeSuiteError` listing available columns (shared `_require_obs_column`).
- `abundance_native(adata, *, level, relative)` returns samples×taxa (index `sample_id`, columns = taxon names), raising on an unknown `level`. `normalize_native(adata, method="clr").X` is the CLR matrix (samples×taxa). `beta_diversity(adata, "bray-curtis")` is a square sample distance frame. `pcoa(dist, dimensions=2)` returns columns `sample_id, PC1, PC2, PC1_variance, PC2_variance`.
- Sample order: `abundance_native`/`dense_counts` preserve `adata.obs_names` order, so grouping by `adata.obs[col]` aligns positionally; align explicitly by `sample_id` where a function returns its own index.
- matplotlib: `matplotlib.use("Agg")` at import, `# noqa: E402` on the `import matplotlib.pyplot as plt` line; every plot does `fig.savefig(output, dpi=160)` then `plt.close(fig)`.
- Fatal → `MicrobiomeSuiteError` (`microsuite._errors`). Both CI gates pass (`ruff check .`, `ruff format --check .`).
- `from __future__ import annotations` at the top of new modules.

---

### Task 1: `viz/metadata.py` scaffold + `taxa-by-group`

**Files:**
- Create: `src/microsuite/viz/metadata.py`
- Modify: `src/microsuite/cli/viz_cmd.py`
- Test: `tests/test_viz_metadata.py`

**Interfaces:**
- Produces: `_require_obs_column(adata, column) -> pd.Series`; `plot_taxa_by_group(adata, *, level, group_by, output, top_n=15) -> None`; CLI `viz taxa-by-group`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_viz_metadata.py
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.io.h5ad import write_h5ad
from microsuite.viz.metadata import plot_taxa_by_group


def make_adata() -> ad.AnnData:
    # 6 samples, 4 features across 2 genera; obs has phase/time/subject
    rng = np.random.default_rng(0)
    X = rng.integers(1, 50, size=(6, 4)).astype(float)
    var = pd.DataFrame(
        {"genus": ["Bacteroides", "Bacteroides", "Prevotella", "Faecalibacterium"]},
        index=["F1", "F2", "F3", "F4"],
    )
    obs = pd.DataFrame(
        {
            "phase": ["pre", "pre", "post", "post", "pre", "post"],
            "time": [0, 0, 7, 7, 0, 7],
            "subject": ["s1", "s2", "s1", "s2", "s3", "s3"],
        },
        index=[f"S{i}" for i in range(6)],
    )
    return ad.AnnData(X=X, obs=obs, var=var)


def test_plot_taxa_by_group_writes_png(tmp_path: Path) -> None:
    out = tmp_path / "taxa.png"
    plot_taxa_by_group(make_adata(), level="genus", group_by="phase", output=out, top_n=2)
    assert out.exists() and out.stat().st_size > 0


def test_plot_taxa_by_group_bad_column(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError, match="nope"):
        plot_taxa_by_group(make_adata(), level="genus", group_by="nope", output=tmp_path / "x.png")


def test_plot_taxa_by_group_bad_level(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        plot_taxa_by_group(make_adata(), level="species", group_by="phase", output=tmp_path / "x.png")


def test_cli_taxa_by_group(tmp_path: Path) -> None:
    src = tmp_path / "t.h5ad"
    write_h5ad(make_adata(), src)
    out = tmp_path / "taxa.png"
    r = CliRunner().invoke(
        app,
        ["viz", "taxa-by-group", "--table", str(src), "--level", "genus",
         "--group-by", "phase", "-o", str(out)],
    )
    assert r.exit_code == 0, r.stdout
    assert out.exists()
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_viz_metadata.py -v`
Expected: FAIL (`ModuleNotFoundError: microsuite.viz.metadata`).

- [ ] **Step 3: Create `src/microsuite/viz/metadata.py`**

```python
from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.abundance import abundance_native

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def _require_obs_column(adata: ad.AnnData, column: str) -> pd.Series:
    if column not in adata.obs.columns:
        available = ", ".join(str(c) for c in adata.obs.columns) or "(none)"
        raise MicrobiomeSuiteError(
            f"Metadata column '{column}' not found in obs; available: {available}."
        )
    return adata.obs[column]


def plot_taxa_by_group(
    adata: ad.AnnData, *, level: str, group_by: str, output: Path, top_n: int = 15
) -> None:
    group = _require_obs_column(adata, group_by)
    frame = abundance_native(adata, level=level, relative=True)  # samples x taxa (raises on level)
    grouped = frame.groupby(group.astype(str).to_numpy()).mean()  # group x taxa
    order = grouped.mean(axis=0).sort_values(ascending=False)
    top = list(order.head(top_n).index)
    plot_df = grouped[top].copy()
    other = grouped.drop(columns=top).sum(axis=1)
    if (other > 0).any():
        plot_df["Other"] = other
    width = max(6.0, len(plot_df.index) * 0.9)
    height = max(4.5, min(12.0, 2.5 + plot_df.shape[1] * 0.25))
    ax = plot_df.plot(kind="bar", stacked=True, figsize=(width, height), width=0.85)
    ax.set_xlabel(group_by)
    ax.set_ylabel("Mean relative abundance")
    ax.set_ylim(0, 1)
    ax.set_title(f"Top {top_n} {level} by {group_by}")
    ax.legend(title=level, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.setp(ax.get_xticklabels(), rotation=0)
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)
```

- [ ] **Step 4: Add the CLI command in `cli/viz_cmd.py`**

Add imports at the top (keep existing):
```python
from microsuite.viz.metadata import plot_taxa_by_group
```
Add the command (mirror the existing `barplot` command's option style):
```python
@app.command("taxa-by-group")
def taxa_by_group(
    table: Annotated[Path, typer.Option("--table", help="Input .h5ad table.")],
    level: Annotated[str, typer.Option("--level", help="Taxonomy level.")],
    group_by: Annotated[str, typer.Option("--group-by", help="Metadata (obs) column to group by.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output PNG.")],
    top_n: Annotated[int, typer.Option("--top-n", min=1)] = 15,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    adata = read_h5ad(ensure_input(table))
    plot_taxa_by_group(
        adata, level=level, group_by=group_by,
        output=prepare_output(output, force=force), top_n=top_n,
    )
```
(`read_h5ad`, `ensure_input`, `prepare_output`, `Annotated`, `typer`, `Path` are already imported in `viz_cmd.py`; confirm and add any missing.)

- [ ] **Step 5: Run to verify pass**

Run: `uv run pytest tests/test_viz_metadata.py -v`
Expected: PASS (4).

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/viz/metadata.py src/microsuite/cli/viz_cmd.py tests/test_viz_metadata.py
git commit -m "feat(viz): taxa-by-group metadata bar plot"
```

---

### Task 2: `clr-by-group` (boxplot / heatmap / violin)

**Files:**
- Modify: `src/microsuite/viz/metadata.py`
- Modify: `src/microsuite/cli/viz_cmd.py`
- Test: `tests/test_viz_metadata.py` (append)

**Interfaces:**
- Consumes: `_require_obs_column`, `make_adata` fixture (Task 1); `normalize_native`, `abundance_native`.
- Produces: `plot_clr_by_group(adata, *, level, group_by, output, top_n=10, style="boxplot") -> None`; CLI `viz clr-by-group`.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_viz_metadata.py`)

```python
def test_plot_clr_by_group_styles(tmp_path: Path) -> None:
    from microsuite.viz.metadata import plot_clr_by_group

    for style in ("boxplot", "heatmap", "violin"):
        out = tmp_path / f"clr_{style}.png"
        plot_clr_by_group(
            make_adata(), level="genus", group_by="phase", output=out, top_n=3, style=style
        )
        assert out.exists() and out.stat().st_size > 0


def test_plot_clr_by_group_bad_style(tmp_path: Path) -> None:
    from microsuite.viz.metadata import plot_clr_by_group

    with pytest.raises(MicrobiomeSuiteError, match="style"):
        plot_clr_by_group(
            make_adata(), level="genus", group_by="phase", output=tmp_path / "x.png", style="pie"
        )


def test_cli_clr_by_group_heatmap(tmp_path: Path) -> None:
    src = tmp_path / "t.h5ad"
    write_h5ad(make_adata(), src)
    out = tmp_path / "clr.png"
    r = CliRunner().invoke(
        app,
        ["viz", "clr-by-group", "--table", str(src), "--level", "genus",
         "--group-by", "phase", "--style", "heatmap", "-o", str(out)],
    )
    assert r.exit_code == 0, r.stdout
    assert out.exists()
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_viz_metadata.py -k clr -v`
Expected: FAIL (`cannot import name 'plot_clr_by_group'`).

- [ ] **Step 3: Add `plot_clr_by_group` to `viz/metadata.py`**

Add the import `from microsuite.methods.normalize import normalize_native` at the top. Then:

```python
_CLR_STYLES = ("boxplot", "heatmap", "violin")


def plot_clr_by_group(
    adata: ad.AnnData, *, level: str, group_by: str, output: Path,
    top_n: int = 10, style: str = "boxplot",
) -> None:
    if style not in _CLR_STYLES:
        raise MicrobiomeSuiteError(
            f"Unknown style '{style}'. Choose one of: {', '.join(_CLR_STYLES)}."
        )
    group = _require_obs_column(adata, group_by).astype(str)
    counts = abundance_native(adata, level=level, relative=False)  # samples x taxa (raises on level)
    rel = abundance_native(adata, level=level, relative=True)
    collapsed = ad.AnnData(
        X=counts.to_numpy(dtype=float),
        obs=pd.DataFrame(index=counts.index),
        var=pd.DataFrame(index=counts.columns),
    )
    clr = normalize_native(collapsed, method="clr").X
    clr_df = pd.DataFrame(clr, index=counts.index, columns=counts.columns)
    top = list(rel.mean(axis=0).sort_values(ascending=False).head(top_n).index)
    clr_top = clr_df[top]
    groups = group.to_numpy()
    categories = sorted(pd.unique(groups))

    if style == "heatmap":
        mean_by_group = clr_top.groupby(groups).mean().reindex(categories)
        matrix = mean_by_group.to_numpy(dtype=float)
        vmax = float(np.nanmax(np.abs(matrix))) or 1.0
        fig, ax = plt.subplots(
            figsize=(max(6.0, len(top) * 0.8), max(3.0, len(categories) * 0.5 + 1.0))
        )
        im = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(len(top)))
        ax.set_xticklabels(top, rotation=45, ha="right")
        ax.set_yticks(range(len(categories)))
        ax.set_yticklabels([str(c) for c in categories])
        ax.set_title(f"Mean CLR of top {top_n} {level} by {group_by}")
        fig.colorbar(im, ax=ax, label="mean CLR")
    else:
        fig, ax = plt.subplots(figsize=(max(7.0, len(top) * 1.3), 5.0))
        n_groups = len(categories)
        slot = 0.8 / max(n_groups, 1)
        cmap = plt.get_cmap("tab10")
        handles = []
        for gi, cat in enumerate(categories):
            series = [clr_top.loc[groups == cat, taxon].to_numpy() for taxon in top]
            series = [d if len(d) else np.array([np.nan]) for d in series]
            positions = np.arange(len(top)) + (gi - (n_groups - 1) / 2) * slot
            color = cmap(gi % 10)
            if style == "boxplot":
                bp = ax.boxplot(
                    series, positions=positions, widths=slot * 0.9,
                    patch_artist=True, manage_ticks=False,
                )
                for box in bp["boxes"]:
                    box.set_facecolor(color)
                    box.set_alpha(0.7)
            else:
                vp = ax.violinplot(series, positions=positions, widths=slot * 0.9, showmeans=True)
                for body in vp["bodies"]:
                    body.set_facecolor(color)
                    body.set_alpha(0.7)
            handles.append(plt.Line2D([0], [0], color=color, lw=6, label=str(cat)))
        ax.set_xticks(range(len(top)))
        ax.set_xticklabels(top, rotation=45, ha="right")
        ax.set_ylabel("CLR")
        ax.set_title(f"CLR of top {top_n} {level} by {group_by}")
        ax.legend(handles=handles, title=group_by, bbox_to_anchor=(1.02, 1), loc="upper left")

    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)
```

- [ ] **Step 4: Add the CLI command in `cli/viz_cmd.py`**

Add `from microsuite.viz.metadata import plot_clr_by_group` (or extend the existing import), then:
```python
@app.command("clr-by-group")
def clr_by_group(
    table: Annotated[Path, typer.Option("--table", help="Input .h5ad table.")],
    level: Annotated[str, typer.Option("--level", help="Taxonomy level.")],
    group_by: Annotated[str, typer.Option("--group-by", help="Metadata (obs) column to group by.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output PNG.")],
    top_n: Annotated[int, typer.Option("--top-n", min=1)] = 10,
    style: Annotated[str, typer.Option("--style", help="boxplot, heatmap, or violin.")] = "boxplot",
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    adata = read_h5ad(ensure_input(table))
    plot_clr_by_group(
        adata, level=level, group_by=group_by,
        output=prepare_output(output, force=force), top_n=top_n, style=style,
    )
```

- [ ] **Step 5: Run to verify pass**

Run: `uv run pytest tests/test_viz_metadata.py -k clr -v`
Expected: PASS (3).

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/viz/metadata.py src/microsuite/cli/viz_cmd.py tests/test_viz_metadata.py
git commit -m "feat(viz): clr-by-group boxplot/heatmap/violin"
```

---

### Task 3: `braycurtis-ordination` (scatter / trajectory / facet)

**Files:**
- Modify: `src/microsuite/viz/metadata.py`
- Modify: `src/microsuite/cli/viz_cmd.py`
- Test: `tests/test_viz_metadata.py` (append)

**Interfaces:**
- Consumes: `_require_obs_column`, `make_adata` fixture; `beta_diversity`, `pcoa`.
- Produces: `plot_braycurtis_ordination(adata, *, color_by, output, subject=None, style=None) -> None`; CLI `viz braycurtis-ordination`.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_viz_metadata.py`)

```python
def test_ordination_scatter_and_trajectory(tmp_path: Path) -> None:
    from microsuite.viz.metadata import plot_braycurtis_ordination

    s = tmp_path / "scatter.png"
    plot_braycurtis_ordination(make_adata(), color_by="phase", output=s)  # scatter default
    assert s.exists() and s.stat().st_size > 0

    t = tmp_path / "traj.png"
    plot_braycurtis_ordination(  # numeric color + subject -> trajectory default
        make_adata(), color_by="time", subject="subject", output=t
    )
    assert t.exists() and t.stat().st_size > 0

    f = tmp_path / "facet.png"
    plot_braycurtis_ordination(
        make_adata(), color_by="phase", subject="subject", output=f, style="facet"
    )
    assert f.exists() and f.stat().st_size > 0


def test_ordination_trajectory_requires_subject(tmp_path: Path) -> None:
    from microsuite.viz.metadata import plot_braycurtis_ordination

    with pytest.raises(MicrobiomeSuiteError, match="subject"):
        plot_braycurtis_ordination(
            make_adata(), color_by="phase", output=tmp_path / "x.png", style="trajectory"
        )


def test_cli_braycurtis_ordination(tmp_path: Path) -> None:
    src = tmp_path / "t.h5ad"
    write_h5ad(make_adata(), src)
    out = tmp_path / "ord.png"
    r = CliRunner().invoke(
        app,
        ["viz", "braycurtis-ordination", "--table", str(src), "--color-by", "time",
         "--subject", "subject", "-o", str(out)],
    )
    assert r.exit_code == 0, r.stdout
    assert out.exists()
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_viz_metadata.py -k ordination -v`
Expected: FAIL (`cannot import name 'plot_braycurtis_ordination'`).

- [ ] **Step 3: Add `plot_braycurtis_ordination` to `viz/metadata.py`**

Add imports at the top: `from microsuite.diversity.beta import beta_diversity` and `from microsuite.ordination.pcoa import pcoa`. Then:

```python
_ORDINATION_STYLES = ("scatter", "trajectory", "facet")


def _pc_coords(adata: ad.AnnData) -> pd.DataFrame:
    dist = beta_diversity(adata, "bray-curtis")
    coords = pcoa(dist, dimensions=2).set_index("sample_id")
    return coords.loc[[str(s) for s in adata.obs_names]]


def plot_braycurtis_ordination(
    adata: ad.AnnData, *, color_by: str, output: Path,
    subject: str | None = None, style: str | None = None,
) -> None:
    color = _require_obs_column(adata, color_by)
    subj = _require_obs_column(adata, subject) if subject is not None else None
    effective = style or ("scatter" if subject is None else "trajectory")
    if effective not in _ORDINATION_STYLES:
        raise MicrobiomeSuiteError(
            f"Unknown style '{effective}'. Choose one of: {', '.join(_ORDINATION_STYLES)}."
        )
    if effective in ("trajectory", "facet") and subj is None:
        raise MicrobiomeSuiteError(f"--subject is required for style '{effective}'.")

    coords = _pc_coords(adata)
    x = coords["PC1"].to_numpy()
    y = coords["PC2"].to_numpy()
    xlab = f"PC1 ({coords['PC1_variance'].iloc[0] * 100:.1f}%)"
    ylab = f"PC2 ({coords['PC2_variance'].iloc[0] * 100:.1f}%)"
    color_vals = color.to_numpy()
    is_numeric = pd.api.types.is_numeric_dtype(color)

    if effective == "facet":
        subjects = list(pd.unique(subj.to_numpy()))
        ncols = min(3, len(subjects))
        nrows = int(np.ceil(len(subjects) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.5 * nrows), squeeze=False)
        for idx, sub in enumerate(subjects):
            ax = axes[idx // ncols][idx % ncols]
            mask = subj.to_numpy() == sub
            sc = ax.scatter(x[mask], y[mask], c=_color_arg(color_vals[mask], is_numeric), s=40)
            ax.set_title(str(sub))
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)
        for j in range(len(subjects), nrows * ncols):
            axes[j // ncols][j % ncols].axis("off")
        fig.tight_layout()
        fig.savefig(output, dpi=160)
        plt.close(fig)
        return

    fig, ax = plt.subplots(figsize=(7.0, 5.5))
    if is_numeric:
        sc = ax.scatter(x, y, c=color_vals.astype(float), cmap="viridis", s=45)
        fig.colorbar(sc, ax=ax, label=color_by)
    else:
        categories = list(pd.unique(color_vals))
        cmap = plt.get_cmap("tab10")
        cat_color = {c: cmap(i % 10) for i, c in enumerate(categories)}
        for cat in categories:
            mask = color_vals == cat
            ax.scatter(x[mask], y[mask], color=cat_color[cat], label=str(cat), s=45)
        ax.legend(title=color_by, bbox_to_anchor=(1.02, 1), loc="upper left")

    if effective == "trajectory":
        subj_vals = subj.to_numpy()
        for sub in pd.unique(subj_vals):
            mask = subj_vals == sub
            order = np.argsort(color_vals[mask]) if is_numeric else np.arange(mask.sum())
            ax.plot(x[mask][order], y[mask][order], color="gray", alpha=0.6, lw=1.0, zorder=0)

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(f"Bray-Curtis PCoA by {color_by}")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def _color_arg(values: np.ndarray, is_numeric: bool):
    if is_numeric:
        return values.astype(float)
    categories = list(pd.unique(values))
    lookup = {c: i for i, c in enumerate(categories)}
    return np.array([lookup[v] for v in values], dtype=float)
```

- [ ] **Step 4: Add the CLI command in `cli/viz_cmd.py`**

Add `from microsuite.viz.metadata import plot_braycurtis_ordination`, then:
```python
@app.command("braycurtis-ordination")
def braycurtis_ordination(
    table: Annotated[Path, typer.Option("--table", help="Input .h5ad table.")],
    color_by: Annotated[str, typer.Option("--color-by", help="Metadata (obs) column to color by.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output PNG.")],
    subject: Annotated[str | None, typer.Option("--subject", help="Obs column of subject IDs.")] = None,
    style: Annotated[str | None, typer.Option("--style", help="scatter, trajectory, or facet.")] = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    adata = read_h5ad(ensure_input(table))
    plot_braycurtis_ordination(
        adata, color_by=color_by, output=prepare_output(output, force=force),
        subject=subject, style=style,
    )
```

- [ ] **Step 5: Run to verify pass + full suite + lint gates**

Run: `uv run pytest tests/test_viz_metadata.py -v` (all pass), then `uv run pytest -q` (all green), then `uv run ruff check .` and `uv run ruff format --check .` (both clean; `uv run ruff format .` on the files you created/changed if needed, then re-check — do not reformat unrelated files). Sanity-check: `uv run python -c "from typer.testing import CliRunner; from microsuite.cli.app import app; print(CliRunner().invoke(app,['viz','--help']).stdout)"` lists `taxa-by-group`, `clr-by-group`, `braycurtis-ordination`.

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/viz/metadata.py src/microsuite/cli/viz_cmd.py tests/test_viz_metadata.py
git commit -m "feat(viz): braycurtis-ordination scatter/trajectory/facet"
```

---

## Self-Review

**Spec coverage:**
- top taxa by time/phase (#18) → Task 1 `taxa-by-group` (`--group-by`). ✓
- CLR by time/phase, boxplot/heatmap/violin (#18) → Task 2 `clr-by-group` (`--group-by`, `--style`). ✓
- subject/time-aware Bray-Curtis (#18) → Task 3 `braycurtis-ordination` (`--color-by`, `--subject`, `--style` scatter/trajectory/facet). ✓
- Clear errors on bad obs column / level / style → `_require_obs_column`, `abundance_native` level check, style guards. ✓
- Both CI gates → Task 3 Step 5. ✓

**Placeholder scan:** none — full plot code for all three commands and every style, exact CLI additions, and complete tests with a concrete fixture.

**Consistency:** `_require_obs_column`, `make_adata`, and the three `plot_*` signatures match between `viz/metadata.py`, its tests, and the CLI; samples×taxa orientation from `abundance_native` throughout; `pcoa` columns (`sample_id`, `PC1/PC2`, `PC*_variance`) match `_pc_coords`'s use; `--style` choice sets (`_CLR_STYLES`, `_ORDINATION_STYLES`) match the CLI help and the tests' style values.
