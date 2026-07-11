# Taxonomy Assignment QC (Round-3 E) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Native taxonomy-assignment QC: a per-`(sample, rank)` `taxonomy_unassigned_summary.tsv` (#17) and three canonical assignment plots (#16), both derived from one computation over `adata.var` (assignment) and `adata.X` (reads).

**Architecture:** A pure `methods/assignment_qc.py` computes the long-format summary + deepest-rank distribution. A `viz/assignment.py` (matplotlib Agg, mirrors `viz/barplot.py`) renders three PNGs from those same functions. Two flat `tax_*` CLI commands in `method_taxonomy_cmd` wire them to `.h5ad` input.

**Tech Stack:** Python 3.12, anndata, numpy, pandas, matplotlib (Agg), Typer, pytest (`CliRunner`).

## Global Constraints

- Taxonomy is feature-level in `adata.var`; ranks are `LEVELS = [kingdom, phylum, class, order, family, genus, species]` (`io/taxonomy.py`); **unassigned at a rank = empty string `""`**. Ranks used = `[r for r in LEVELS if r in adata.var.columns]`; none present → `MicrobiomeSuiteError`.
- `dense_counts(adata)` is samples×features (`obs_names`×`var_names`). A feature "counts" for a sample when its reads there are `> 0`.
- Pooled rows use `sample == "_all_samples_"` (constant `POOLED_LABEL`); pooled feature presence = column sum `> 0`, pooled reads = column sums.
- Fracs = `assigned / (assigned + unassigned)`, `0.0` when denominator is 0.
- Plots and TSV both call the Component-1 functions so they never disagree. matplotlib uses the `Agg` backend; each plot `savefig` then `plt.close(fig)`.
- Fatal → `MicrobiomeSuiteError` (`microsuite._errors`). Both CI gates pass (`ruff check .`, `ruff format --check .`).
- `from __future__ import annotations` at the top of new modules.

---

### Task 1: `methods/assignment_qc.py` — summary + deepest-rank

**Files:**
- Create: `src/microsuite/methods/assignment_qc.py`
- Test: `tests/test_assignment_qc.py`

**Interfaces:**
- Produces:
  - `summarize_assignment(adata: ad.AnnData) -> pd.DataFrame`
  - `write_assignment_summary(summary: pd.DataFrame, out_path: Path) -> Path`
  - `deepest_rank_distribution(adata: ad.AnnData) -> pd.Series`
  - `POOLED_LABEL = "_all_samples_"`, `SUMMARY_COLUMNS` (list)

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_assignment_qc.py
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.io.taxonomy import LEVELS
from microsuite.methods.assignment_qc import (
    POOLED_LABEL,
    deepest_rank_distribution,
    summarize_assignment,
    write_assignment_summary,
)


def _fixture() -> ad.AnnData:
    # F1 assigned to species (all ranks), F2 assigned to genus only, F3 unassigned.
    rank_values = {
        "kingdom": ["Bacteria", "Bacteria", ""],
        "phylum": ["Firmicutes", "Bacteroidetes", ""],
        "class": ["Bacilli", "Bacteroidia", ""],
        "order": ["Lactobacillales", "Bacteroidales", ""],
        "family": ["Lactobacillaceae", "Prevotellaceae", ""],
        "genus": ["Lactobacillus", "Prevotella", ""],
        "species": ["L. casei", "", ""],
    }
    var = pd.DataFrame(rank_values, index=["F1", "F2", "F3"])
    # samples s1 (all present), s2 (F1 absent)
    X = np.array([[10.0, 5.0, 1.0], [0.0, 3.0, 2.0]])
    return ad.AnnData(X=X, obs=pd.DataFrame(index=["s1", "s2"]), var=var)


def test_summarize_assignment_overall_counts() -> None:
    df = summarize_assignment(_fixture())
    pooled = df[df["sample"] == POOLED_LABEL].set_index("rank")
    # species: only F1 assigned -> 1 assigned, 2 unassigned features
    assert pooled.loc["species", "assigned_features"] == 1
    assert pooled.loc["species", "unassigned_features"] == 2
    # genus: F1 + F2 assigned -> 2 assigned, 1 unassigned
    assert pooled.loc["genus", "assigned_features"] == 2
    assert pooled.loc["genus", "unassigned_features"] == 1
    # pooled reads at species: assigned = F1 reads (10) ; unassigned = F2+F3 (5+3+1+2)=11
    assert pooled.loc["species", "assigned_reads"] == 10.0
    assert pooled.loc["species", "unassigned_reads"] == 11.0
    # assigned_features is non-increasing kingdom -> species (nested assignment)
    ranks = [r for r in LEVELS]
    seq = [pooled.loc[r, "assigned_features"] for r in ranks]
    assert all(seq[i] >= seq[i + 1] for i in range(len(seq) - 1))


def test_summarize_assignment_per_sample_presence() -> None:
    df = summarize_assignment(_fixture())
    s2 = df[df["sample"] == "s2"].set_index("rank")
    # in s2, F1 absent (0 reads); present = {F2, F3}
    # genus: F2 assigned, F3 not -> 1 assigned, 1 unassigned
    assert s2.loc["genus", "assigned_features"] == 1
    assert s2.loc["genus", "unassigned_features"] == 1
    # s2 reads: total 5? no -> F2=3, F3=2 => genus assigned_reads=3, unassigned=2
    assert s2.loc["genus", "assigned_reads"] == 3.0
    assert s2.loc["genus", "unassigned_reads"] == 2.0
    assert s2.loc["genus", "assigned_read_frac"] == pytest.approx(3.0 / 5.0)


def test_deepest_rank_distribution() -> None:
    dist = deepest_rank_distribution(_fixture())
    assert dist["species"] == 1  # F1
    assert dist["genus"] == 1  # F2
    assert dist["Unassigned"] == 1  # F3
    assert int(dist.sum()) == 3  # each ASV counted once


def test_no_taxonomy_ranks_raises() -> None:
    adata = ad.AnnData(
        X=np.array([[1.0, 2.0]]), obs=pd.DataFrame(index=["s1"]),
        var=pd.DataFrame(index=["F1", "F2"]),
    )
    with pytest.raises(MicrobiomeSuiteError):
        summarize_assignment(adata)


def test_write_assignment_summary_roundtrip(tmp_path: Path) -> None:
    df = summarize_assignment(_fixture())
    out = write_assignment_summary(df, tmp_path / "s.tsv")
    assert out.exists()
    back = pd.read_csv(out, sep="\t")
    assert set(["sample", "rank", "assigned_features", "assigned_read_frac"]).issubset(back.columns)
    assert (back["sample"] == POOLED_LABEL).any()
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_assignment_qc.py -v`
Expected: FAIL (`ModuleNotFoundError: microsuite.methods.assignment_qc`).

- [ ] **Step 3: Create `src/microsuite/methods/assignment_qc.py`**

```python
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity._matrix import dense_counts
from microsuite.io.taxonomy import LEVELS

POOLED_LABEL = "_all_samples_"
SUMMARY_COLUMNS = [
    "sample",
    "rank",
    "assigned_features",
    "unassigned_features",
    "assigned_reads",
    "unassigned_reads",
    "assigned_feature_frac",
    "assigned_read_frac",
]


def _ranks(adata: ad.AnnData) -> list[str]:
    ranks = [r for r in LEVELS if r in adata.var.columns]
    if not ranks:
        raise MicrobiomeSuiteError(
            "Table has no taxonomy rank columns; run tax_classify or import with "
            "--taxonomy first."
        )
    return ranks


def _assigned_masks(adata: ad.AnnData, ranks: list[str]) -> dict[str, np.ndarray]:
    return {r: (adata.var[r].astype(str).to_numpy() != "") for r in ranks}


def _row(sample: str, present: np.ndarray, reads: np.ndarray, assigned: np.ndarray, rank: str) -> list:
    af = int((present & assigned).sum())
    uf = int((present & ~assigned).sum())
    ar = float(reads[assigned].sum())
    ur = float(reads[~assigned].sum())
    ff = af / (af + uf) if (af + uf) else 0.0
    rf = ar / (ar + ur) if (ar + ur) else 0.0
    return [sample, rank, af, uf, ar, ur, round(ff, 6), round(rf, 6)]


def summarize_assignment(adata: ad.AnnData) -> pd.DataFrame:
    ranks = _ranks(adata)
    counts = dense_counts(adata)  # samples x features
    assigned = _assigned_masks(adata, ranks)
    rows: list[list] = []
    for i, sample in enumerate([str(s) for s in adata.obs_names]):
        reads = counts[i]
        present = reads > 0
        for rank in ranks:
            rows.append(_row(sample, present, reads, assigned[rank], rank))
    pooled_reads = counts.sum(axis=0)
    pooled_present = pooled_reads > 0
    for rank in ranks:
        rows.append(_row(POOLED_LABEL, pooled_present, pooled_reads, assigned[rank], rank))
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


def write_assignment_summary(summary: pd.DataFrame, out_path: Path) -> Path:
    summary.to_csv(out_path, sep="\t", index=False)
    return out_path


def deepest_rank_distribution(adata: ad.AnnData) -> pd.Series:
    ranks = _ranks(adata)
    deepest = pd.Series("Unassigned", index=adata.var.index)
    for rank in ranks:  # shallow -> deep; later ranks overwrite
        mask = adata.var[rank].astype(str).to_numpy() != ""
        deepest[mask] = rank
    return deepest.value_counts().reindex([*ranks, "Unassigned"], fill_value=0)
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_assignment_qc.py -v`
Expected: PASS (5).

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/methods/assignment_qc.py tests/test_assignment_qc.py
git commit -m "feat(taxonomy): assignment QC summary + deepest-rank helpers"
```

---

### Task 2: plots (`viz/assignment.py`) + CLI (`method_taxonomy_cmd`)

**Files:**
- Create: `src/microsuite/viz/assignment.py`
- Modify: `src/microsuite/cli/method_taxonomy_cmd.py`
- Test: `tests/test_assignment_plots_cli.py`

**Interfaces:**
- Consumes: `summarize_assignment`, `deepest_rank_distribution`, `POOLED_LABEL` (Task 1); `read_h5ad`, `ensure_input`, `prepare_output`.
- Produces: `plot_assigned_asv_by_rank`, `plot_assigned_reads_by_rank`, `plot_deepest_rank` (each `(adata, output: Path) -> None`); CLI `tax_assignment_summary`, `tax_assignment_plots`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_assignment_plots_cli.py
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from typer.testing import CliRunner

from microsuite.cli.app import app
from microsuite.io.h5ad import write_h5ad
from microsuite.viz.assignment import (
    plot_assigned_asv_by_rank,
    plot_assigned_reads_by_rank,
    plot_deepest_rank,
)


def _fixture() -> ad.AnnData:
    var = pd.DataFrame(
        {
            "kingdom": ["Bacteria", "Bacteria", ""],
            "phylum": ["Firmicutes", "Bacteroidetes", ""],
            "class": ["Bacilli", "Bacteroidia", ""],
            "order": ["Lactobacillales", "Bacteroidales", ""],
            "family": ["Lactobacillaceae", "Prevotellaceae", ""],
            "genus": ["Lactobacillus", "Prevotella", ""],
            "species": ["L. casei", "", ""],
        },
        index=["F1", "F2", "F3"],
    )
    X = np.array([[10.0, 5.0, 1.0], [0.0, 3.0, 2.0]])
    return ad.AnnData(X=X, obs=pd.DataFrame(index=["s1", "s2"]), var=var)


def test_plot_functions_write_pngs(tmp_path: Path) -> None:
    adata = _fixture()
    for fn, name in (
        (plot_assigned_asv_by_rank, "asv.png"),
        (plot_assigned_reads_by_rank, "reads.png"),
        (plot_deepest_rank, "deepest.png"),
    ):
        out = tmp_path / name
        fn(adata, out)
        assert out.exists() and out.stat().st_size > 0


def test_plot_reads_single_sample(tmp_path: Path) -> None:
    adata = _fixture()[[0]].copy()  # one sample
    out = tmp_path / "reads.png"
    plot_assigned_reads_by_rank(adata, out)
    assert out.exists() and out.stat().st_size > 0


def test_cli_assignment_summary_and_plots(tmp_path: Path) -> None:
    src = tmp_path / "t.h5ad"
    write_h5ad(_fixture(), src)
    runner = CliRunner()

    summary = tmp_path / "summary.tsv"
    r1 = runner.invoke(app, ["tax_assignment_summary", "--table", str(src), "-o", str(summary)])
    assert r1.exit_code == 0, r1.stdout
    assert summary.exists()

    plots = tmp_path / "plots"
    r2 = runner.invoke(app, ["tax_assignment_plots", "--table", str(src), "--output-dir", str(plots)])
    assert r2.exit_code == 0, r2.stdout
    for name in ("assigned_asv_by_rank.png", "assigned_reads_by_rank.png", "deepest_rank.png"):
        assert (plots / name).exists()
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_assignment_plots_cli.py -v`
Expected: FAIL (`ModuleNotFoundError: microsuite.viz.assignment`).

- [ ] **Step 3: Create `src/microsuite/viz/assignment.py`**

```python
from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np

from microsuite.io.taxonomy import LEVELS
from microsuite.methods.assignment_qc import (
    POOLED_LABEL,
    deepest_rank_distribution,
    summarize_assignment,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def plot_assigned_asv_by_rank(adata: ad.AnnData, output: Path) -> None:
    summary = summarize_assignment(adata)
    pooled = summary[summary["sample"] == POOLED_LABEL]
    ranks = list(pooled["rank"])
    total = (pooled["assigned_features"] + pooled["unassigned_features"]).to_numpy(dtype=float)
    assigned_frac = np.divide(
        pooled["assigned_features"].to_numpy(dtype=float), total,
        out=np.zeros_like(total), where=total > 0,
    )
    unassigned_frac = 1.0 - assigned_frac
    fig, ax = plt.subplots(figsize=(max(6.0, len(ranks) * 1.1), 4.5))
    ax.bar(ranks, assigned_frac, label="assigned")
    ax.bar(ranks, unassigned_frac, bottom=assigned_frac, label="unassigned")
    ax.set_ylim(0, 1)
    ax.set_xlabel("Rank")
    ax.set_ylabel("Feature fraction")
    ax.set_title("Assigned vs unassigned ASVs by rank")
    ax.legend()
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def plot_assigned_reads_by_rank(adata: ad.AnnData, output: Path) -> None:
    summary = summarize_assignment(adata)
    per_sample = summary[summary["sample"] != POOLED_LABEL]
    pivot = per_sample.pivot(index="sample", columns="rank", values="assigned_read_frac")
    ranks = [r for r in LEVELS if r in pivot.columns]
    pivot = pivot[ranks]
    fig, ax = plt.subplots(
        figsize=(max(6.0, len(ranks) * 1.1), max(3.0, len(pivot.index) * 0.5 + 1.0))
    )
    im = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(ranks)))
    ax.set_xticklabels(ranks, rotation=45, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels([str(s) for s in pivot.index])
    ax.set_title("Assigned read fraction by rank")
    fig.colorbar(im, ax=ax, label="assigned read fraction")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def plot_deepest_rank(adata: ad.AnnData, output: Path) -> None:
    dist = deepest_rank_distribution(adata)
    fig, ax = plt.subplots(figsize=(max(6.0, len(dist) * 1.1), 4.5))
    ax.bar([str(i) for i in dist.index], dist.to_numpy(dtype=float))
    ax.set_xlabel("Deepest assigned rank")
    ax.set_ylabel("Feature count")
    ax.set_title("Deepest assigned rank per ASV")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)
```

- [ ] **Step 4: Add the two CLI commands**

At the top of `cli/method_taxonomy_cmd.py`, add imports (keep existing ones):
```python
from microsuite._paths import ensure_input, prepare_output
from microsuite.io.h5ad import read_h5ad
from microsuite.methods.assignment_qc import summarize_assignment, write_assignment_summary
from microsuite.viz.assignment import (
    plot_assigned_asv_by_rank,
    plot_assigned_reads_by_rank,
    plot_deepest_rank,
)
```

Inside `register(app)`, add:
```python
    @app.command("tax_assignment_summary")
    def tax_assignment_summary_cmd(
        table: Annotated[Path, typer.Option("--table", help="Input .h5ad with taxonomy.")],
        output: Annotated[Path, typer.Option("--output", "-o", help="Output summary TSV.")],
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
    ) -> None:
        adata = read_h5ad(ensure_input(table))
        summary = summarize_assignment(adata)
        write_assignment_summary(summary, prepare_output(output, force=force))

    @app.command("tax_assignment_plots")
    def tax_assignment_plots_cmd(
        table: Annotated[Path, typer.Option("--table", help="Input .h5ad with taxonomy.")],
        output_dir: Annotated[Path, typer.Option("--output-dir", help="Directory for the 3 PNGs.")],
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
    ) -> None:
        adata = read_h5ad(ensure_input(table))
        output_dir.mkdir(parents=True, exist_ok=True)
        plot_assigned_asv_by_rank(
            adata, prepare_output(output_dir / "assigned_asv_by_rank.png", force=force)
        )
        plot_assigned_reads_by_rank(
            adata, prepare_output(output_dir / "assigned_reads_by_rank.png", force=force)
        )
        plot_deepest_rank(adata, prepare_output(output_dir / "deepest_rank.png", force=force))
```

(These commands are self-contained — they don't go through `_method_api`, since there is no backend dispatch. `method_taxonomy_cmd.register(app)` is already wired into the app.)

- [ ] **Step 5: Run to verify pass**

Run: `uv run pytest tests/test_assignment_plots_cli.py -v`
Expected: PASS (3). Sanity-check the group lists the new commands:
`uv run python -c "from typer.testing import CliRunner; from microsuite.cli.app import app; print(CliRunner().invoke(app, ['--help']).stdout)"` (or the method group help) shows `tax_assignment_summary` / `tax_assignment_plots`.

- [ ] **Step 6: Full suite + lint gates**

Run: `uv run pytest -q` (all green), then `uv run ruff check .` and `uv run ruff format --check .` (both clean; `uv run ruff format .` on the files you created if needed, then re-check — do not reformat unrelated files).

- [ ] **Step 7: Commit**

```bash
git add src/microsuite/viz/assignment.py src/microsuite/cli/method_taxonomy_cmd.py tests/test_assignment_plots_cli.py
git commit -m "feat(taxonomy): tax_assignment_summary + tax_assignment_plots commands"
```

---

## Self-Review

**Spec coverage:**
- Native `taxonomy_unassigned_summary.tsv`, per-(sample,rank) + pooled (#17) → Task 1 `summarize_assignment`/`write_assignment_summary` + Task 2 `tax_assignment_summary`. ✓
- Three canonical plots to `--output-dir` (#16): assigned ASV% by rank, per-sample assigned-read heatmap, deepest-rank distribution → Task 2 `viz/assignment.py` + `tax_assignment_plots`. ✓
- Plots and TSV from the same functions → both call Component 1. ✓
- Assignment = non-empty rank; ranks = LEVELS present; no-taxonomy error → Task 1 `_ranks`/`_assigned_masks`. ✓
- Both CI gates → Task 2 Step 6. ✓

**Placeholder scan:** none — full module code, full plot code, exact CLI additions, and complete tests with a concrete fixture and asserted numbers.

**Consistency:** `POOLED_LABEL`, `SUMMARY_COLUMNS`, and the `summarize_assignment`/`deepest_rank_distribution` signatures match between Task 1's module, its tests, Task 2's `viz/assignment.py`, and the CLI; matrices are samples×features from `dense_counts` throughout; the three PNG basenames (`assigned_asv_by_rank.png`, `assigned_reads_by_rank.png`, `deepest_rank.png`) match between the CLI command and its test.
