# DADA2 Parameter-Sensitivity Sweep (Round-2 C) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A `denoise-sweep` helper that runs the `dada2-r` backend across a small parameter grid and writes `dada2_sweep_summary.tsv` comparing retention, ASV count, sample depth, observed ASVs, chimera fraction, and baseline similarity.

**Architecture:** `methods/dada2_sweep.py` holds pure metrics/comparison/summary logic + a grid builder (JSON config or CLI axes) + a thin orchestration wrapping `denoise()`. A `denoise-sweep` CLI command drives it. Metrics/grid are fully offline-tested; orchestration is tested with a stubbed `denoise`.

**Tech Stack:** Python 3.12, pandas, numpy, scipy.stats, anndata-free, Typer, pytest (`CliRunner`).

## Global Constraints

- ASV table on disk is features×samples: first column = ASV id (written with `col.names=NA`, leading empty header cell), header row = sample ids, values = counts. Read via `pd.read_csv(path, sep="\t", index_col=0)`.
- rep-seqs FASTA: `>ASVid\n<sequence>`; ASVs match across runs by **sequence**.
- Reuse `dada2_qc.summarize_dada2_stats(stats_path)` for `filtered_frac`/`merged_frac`/`nonchim_frac`. `chimera_frac = (pre − nonchim)/pre` from stats column totals, `pre` = `merged` (paired) or `denoised` (single), `0.0` when `pre == 0`.
- Correlations (`scipy.stats.pearsonr`/`spearmanr`) return `float("nan")` when undefined (< 2 points or zero variance in either vector).
- Exactly one baseline per grid. CLI-axes baseline = grid point 0 (first value of every axis). Config must flag exactly one `baseline: true`.
- `denoise()` reads from its `demux` argument (the sweep passes `demux=input_dir`). A non-baseline point that raises `MicrobiomeSuiteError` is recorded `status="failed"` (warn, continue); a baseline failure raises.
- Fatal → `MicrobiomeSuiteError` (`microsuite._errors`); non-fatal → `warnings.warn`. Both CI gates pass (`ruff check .`, `ruff format --check .`).
- `from __future__ import annotations` at the top of the new module.

---

### Task 1: analytical core (`methods/dada2_sweep.py`)

**Files:**
- Create: `src/microsuite/methods/dada2_sweep.py`
- Test: `tests/test_dada2_sweep_metrics.py`

**Interfaces:**
- Produces: `GridPoint`, `SweepRun` dataclasses; `run_metrics(table_path, stats_path) -> dict`; `compare_to_baseline(*, table_path, rep_seqs_path, baseline_table_path, baseline_rep_seqs_path) -> dict`; `summarize_sweep(runs: list[SweepRun]) -> pd.DataFrame`; `write_sweep_summary(summary, out_path) -> Path`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_dada2_sweep_metrics.py
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.dada2_sweep import (
    GridPoint,
    SweepRun,
    compare_to_baseline,
    run_metrics,
    summarize_sweep,
    write_sweep_summary,
)


def _write_table(path: Path, asvs: dict[str, list[int]], samples: list[str]) -> Path:
    # asvs: {asv_id: [counts per sample]}; features x samples with leading empty header
    frame = pd.DataFrame(asvs, index=samples).T  # asv x sample
    frame.columns = samples
    frame.to_csv(path, sep="\t")  # index label empty-ish; index_col=0 reads it back
    return path


def _write_fasta(path: Path, seqs: dict[str, str]) -> Path:
    path.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()), encoding="utf-8")
    return path


def _write_stats(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


PAIRED_STATS = (
    "\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\n"
    "s1\t1000\t900\t880\t870\t800\t700\n"
    "s2\t2000\t1800\t1700\t1690\t1600\t1500\n"
)


def test_run_metrics(tmp_path: Path) -> None:
    table = _write_table(tmp_path / "t.tsv", {"ASV1": [5, 1], "ASV2": [0, 3]}, ["s1", "s2"])
    stats = _write_stats(tmp_path / "s.tsv", PAIRED_STATS)
    m = run_metrics(table, stats)
    assert m["n_asvs"] == 2
    # depth: s1=5, s2=4 -> mean 4.5
    assert m["mean_sample_depth"] == pytest.approx(4.5)
    # observed: s1 has ASV1 (1), s2 has ASV1+ASV2 (2) -> mean 1.5
    assert m["mean_observed_asvs"] == pytest.approx(1.5)
    # chimera: pre=merged total 2400, nonchim 2200 -> (2400-2200)/2400
    assert m["chimera_frac"] == pytest.approx((2400 - 2200) / 2400, abs=1e-4)
    assert m["nonchim_frac"] == pytest.approx(2200 / 3000, abs=1e-4)


def test_compare_to_baseline_shared_sequences(tmp_path: Path) -> None:
    # baseline ASVs: seqAAA, seqCCC ; run ASVs: seqAAA (shared), seqGGG (new)
    bt = _write_table(tmp_path / "bt.tsv", {"B1": [10, 0], "B2": [0, 20]}, ["s1", "s2"])
    bf = _write_fasta(tmp_path / "bf.fasta", {"B1": "AAA", "B2": "CCC"})
    rt = _write_table(tmp_path / "rt.tsv", {"R1": [8, 0], "R2": [0, 5]}, ["s1", "s2"])
    rf = _write_fasta(tmp_path / "rf.fasta", {"R1": "AAA", "R2": "GGG"})
    c = compare_to_baseline(
        table_path=rt, rep_seqs_path=rf, baseline_table_path=bt, baseline_rep_seqs_path=bf
    )
    assert c["shared_asv_count"] == 1  # only seqAAA
    # baseline reads in shared (seqAAA total=10) / total baseline reads (30)
    assert c["frac_baseline_reads_shared"] == pytest.approx(10 / 30, abs=1e-4)


def test_compare_to_baseline_disjoint(tmp_path: Path) -> None:
    bt = _write_table(tmp_path / "bt.tsv", {"B1": [10]}, ["s1"])
    bf = _write_fasta(tmp_path / "bf.fasta", {"B1": "AAA"})
    rt = _write_table(tmp_path / "rt.tsv", {"R1": [8]}, ["s1"])
    rf = _write_fasta(tmp_path / "rf.fasta", {"R1": "TTT"})
    c = compare_to_baseline(
        table_path=rt, rep_seqs_path=rf, baseline_table_path=bt, baseline_rep_seqs_path=bf
    )
    assert c["shared_asv_count"] == 0
    assert np.isnan(c["abundance_pearson"])


def _run(tmp_path, name, asvs, seqs, stats_text, baseline=False) -> SweepRun:
    t = _write_table(tmp_path / f"{name}_t.tsv", asvs, ["s1", "s2"])
    f = _write_fasta(tmp_path / f"{name}_f.fasta", seqs)
    s = _write_stats(tmp_path / f"{name}_s.tsv", stats_text)
    return SweepRun(point=GridPoint(name=name, params={"max_ee_f": 2}, is_baseline=baseline),
                    table=t, rep_seqs=f, stats=s)


def test_summarize_sweep_baseline_first(tmp_path: Path) -> None:
    base = _run(tmp_path, "baseline", {"B1": [10, 5]}, {"B1": "AAA"}, PAIRED_STATS, baseline=True)
    variant = _run(tmp_path, "relaxed", {"R1": [8, 4]}, {"R1": "AAA"}, PAIRED_STATS)
    df = summarize_sweep([variant, base])  # order shouldn't matter
    assert list(df["name"]) == ["baseline", "relaxed"]  # baseline first
    assert df.iloc[0]["is_baseline"]
    assert "n_asvs" in df.columns and "shared_asv_count" in df.columns
    out = write_sweep_summary(df, tmp_path / "summary.tsv")
    assert out.exists() and "name" in out.read_text().splitlines()[0]


def test_summarize_sweep_no_baseline_raises(tmp_path: Path) -> None:
    v = _run(tmp_path, "a", {"R1": [1, 1]}, {"R1": "AAA"}, PAIRED_STATS)
    with pytest.raises(MicrobiomeSuiteError):
        summarize_sweep([v])


def test_summarize_sweep_failed_run_is_nan_row(tmp_path: Path) -> None:
    base = _run(tmp_path, "baseline", {"B1": [10, 5]}, {"B1": "AAA"}, PAIRED_STATS, baseline=True)
    failed = SweepRun(
        point=GridPoint(name="bad", params={"max_ee_f": 1}, is_baseline=False),
        table=tmp_path / "missing_t.tsv", rep_seqs=tmp_path / "missing_f.fasta",
        stats=tmp_path / "missing_s.tsv", status="failed",
    )
    df = summarize_sweep([base, failed])
    bad_row = df[df["name"] == "bad"].iloc[0]
    assert bad_row["status"] == "failed"
    assert pd.isna(bad_row["n_asvs"])
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_dada2_sweep_metrics.py -v`
Expected: FAIL (`ModuleNotFoundError: microsuite.methods.dada2_sweep`).

- [ ] **Step 3: Create `src/microsuite/methods/dada2_sweep.py`**

```python
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.dada2_qc import summarize_dada2_stats


@dataclass(frozen=True)
class GridPoint:
    name: str
    params: dict = field(default_factory=dict)
    is_baseline: bool = False


@dataclass
class SweepRun:
    point: GridPoint
    table: Path
    rep_seqs: Path
    stats: Path
    status: str = "ok"


def _read_asv_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    return frame


def _read_rep_seqs(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    current: str | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            seqs.setdefault(current, "")
        elif line.strip() and current is not None:
            seqs[current] += line.strip()
    return seqs


def _corr(a: np.ndarray, b: np.ndarray, *, method: str) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or len(b) < 2 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    fn = scipy_stats.pearsonr if method == "pearson" else scipy_stats.spearmanr
    return float(fn(a, b)[0])


def run_metrics(table_path: Path, stats_path: Path) -> dict:
    table = _read_asv_table(table_path)
    summary = summarize_dada2_stats(stats_path)
    overall = summary["overall"]
    depths = table.sum(axis=0)
    observed = (table > 0).sum(axis=0)

    stats = pd.read_csv(stats_path, sep="\t", index_col=0)
    totals = stats.sum(axis=0)
    pre_col = "merged" if "merged" in totals.index else "denoised"
    pre = float(totals.get(pre_col, 0.0))
    nonchim = float(totals.get("nonchim", 0.0))
    chimera_frac = (pre - nonchim) / pre if pre > 0 else 0.0

    return {
        "n_asvs": int(table.shape[0]),
        "filtered_frac": overall["filtered_frac"],
        "merged_frac": overall["merged_frac"],
        "nonchim_frac": overall["nonchim_frac"],
        "chimera_frac": round(chimera_frac, 6),
        "mean_sample_depth": round(float(depths.mean()), 3),
        "median_sample_depth": round(float(depths.median()), 3),
        "mean_observed_asvs": round(float(observed.mean()), 3),
    }


def _seq_abundance(table_path: Path, rep_seqs_path: Path) -> pd.DataFrame:
    table = _read_asv_table(table_path)
    seqs = _read_rep_seqs(rep_seqs_path)
    keep = [asv for asv in table.index if asv in seqs]
    table = table.loc[keep]
    table.index = [seqs[asv] for asv in keep]
    return table.groupby(level=0).sum()  # sequence x sample


def compare_to_baseline(
    *, table_path: Path, rep_seqs_path: Path,
    baseline_table_path: Path, baseline_rep_seqs_path: Path,
) -> dict:
    run_seq = _seq_abundance(table_path, rep_seqs_path)
    base_seq = _seq_abundance(baseline_table_path, baseline_rep_seqs_path)
    base_set = set(base_seq.index)
    shared = [s for s in run_seq.index if s in base_set]

    total_base = base_seq.sum(axis=1)
    base_reads = float(total_base.sum())
    shared_reads = float(total_base.loc[shared].sum()) if shared else 0.0
    frac_shared = shared_reads / base_reads if base_reads > 0 else 0.0

    if shared:
        run_vec = run_seq.loc[shared].sum(axis=1).to_numpy()
        base_vec = base_seq.loc[shared].sum(axis=1).to_numpy()
        ab_p = _corr(run_vec, base_vec, method="pearson")
        ab_s = _corr(run_vec, base_vec, method="spearman")
    else:
        ab_p = ab_s = float("nan")

    run_full = _read_asv_table(table_path)
    base_full = _read_asv_table(baseline_table_path)
    samples = [c for c in run_full.columns if c in set(base_full.columns)]
    depth_p = _corr(
        run_full[samples].sum(axis=0).to_numpy(),
        base_full[samples].sum(axis=0).to_numpy(), method="pearson",
    )
    obs_p = _corr(
        (run_full[samples] > 0).sum(axis=0).to_numpy(),
        (base_full[samples] > 0).sum(axis=0).to_numpy(), method="pearson",
    )

    return {
        "shared_asv_count": len(shared),
        "frac_baseline_reads_shared": round(frac_shared, 6),
        "abundance_pearson": ab_p,
        "abundance_spearman": ab_s,
        "depth_pearson": depth_p,
        "observed_asv_pearson": obs_p,
    }


def summarize_sweep(runs: list[SweepRun]) -> pd.DataFrame:
    baseline = next((r for r in runs if r.point.is_baseline), None)
    if baseline is None:
        raise MicrobiomeSuiteError("Sweep has no baseline run.")
    ordered = [baseline] + [r for r in runs if not r.point.is_baseline]
    rows: list[dict] = []
    for run in ordered:
        row: dict = {"name": run.point.name, "is_baseline": run.point.is_baseline, "status": run.status}
        row.update({f"param_{k}": v for k, v in run.point.params.items()})
        if run.status == "ok":
            row.update(run_metrics(run.table, run.stats))
            if baseline.status == "ok":
                row.update(compare_to_baseline(
                    table_path=run.table, rep_seqs_path=run.rep_seqs,
                    baseline_table_path=baseline.table, baseline_rep_seqs_path=baseline.rep_seqs,
                ))
        rows.append(row)
    return pd.DataFrame(rows)


def write_sweep_summary(summary: pd.DataFrame, out_path: Path) -> Path:
    summary.to_csv(out_path, sep="\t", index=False)
    return out_path
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_dada2_sweep_metrics.py -v`
Expected: PASS (7). If `_write_table`'s round-trip header differs (pandas writes the index name row), the `index_col=0` read still yields ASV ids as index and samples as columns — the tests assert values, so confirm they pass; do not change the on-disk contract (`index_col=0`).

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/methods/dada2_sweep.py tests/test_dada2_sweep_metrics.py
git commit -m "feat(dada2): sweep metrics, baseline comparison, summary"
```

---

### Task 2: grid builder (config + CLI axes)

**Files:**
- Modify: `src/microsuite/methods/dada2_sweep.py`
- Test: `tests/test_dada2_sweep_grid.py`

**Interfaces:**
- Consumes: `GridPoint` (Task 1).
- Produces: `build_grid(*, config: Path | None = None, axes: dict[str, list] | None = None) -> list[GridPoint]`; `SWEEP_AXES` tuple.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_dada2_sweep_grid.py
from __future__ import annotations

import json
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.dada2_sweep import build_grid


def test_grid_from_config(tmp_path: Path) -> None:
    cfg = tmp_path / "grid.json"
    cfg.write_text(json.dumps([
        {"name": "baseline", "baseline": True, "params": {"max_ee_f": 2, "max_ee_r": 2}},
        {"name": "relaxed", "params": {"max_ee_f": 3, "max_ee_r": 5}},
    ]), encoding="utf-8")
    grid = build_grid(config=cfg)
    assert [p.name for p in grid] == ["baseline", "relaxed"]
    assert sum(p.is_baseline for p in grid) == 1
    assert grid[0].params == {"max_ee_f": 2, "max_ee_r": 2}


def test_grid_config_requires_exactly_one_baseline(tmp_path: Path) -> None:
    cfg = tmp_path / "grid.json"
    cfg.write_text(json.dumps([{"name": "a", "params": {}}, {"name": "b", "params": {}}]), encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match="baseline"):
        build_grid(config=cfg)


def test_grid_from_axes_cartesian(tmp_path: Path) -> None:
    grid = build_grid(axes={"max_ee_f": [2, 3], "trunc_len_f": [0, 220]})
    assert len(grid) == 4  # 2 x 2
    assert sum(p.is_baseline for p in grid) == 1
    baseline = next(p for p in grid if p.is_baseline)
    assert baseline.params == {"max_ee_f": 2, "trunc_len_f": 0}  # first value of each axis


def test_grid_both_or_neither_errors(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        build_grid()
    cfg = tmp_path / "g.json"
    cfg.write_text("[]", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError):
        build_grid(config=cfg, axes={"max_ee_f": [2]})
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_dada2_sweep_grid.py -v`
Expected: FAIL (`cannot import name 'build_grid'`).

- [ ] **Step 3: Add `build_grid` to `methods/dada2_sweep.py`**

Add `import json` and `from itertools import product` at the top. Then:

```python
SWEEP_AXES = (
    "max_ee_f", "max_ee_r", "trunc_len_f", "trunc_len_r",
    "trunc_q", "min_overlap", "max_ee", "trunc_len",
)


def _grid_from_config(config: Path) -> list[GridPoint]:
    try:
        entries = json.loads(config.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        raise MicrobiomeSuiteError(f"Could not read grid config {config}: {exc}") from exc
    if not isinstance(entries, list) or not entries:
        raise MicrobiomeSuiteError("Grid config must be a non-empty JSON list of param sets.")
    points: list[GridPoint] = []
    names: set[str] = set()
    for entry in entries:
        name = str(entry.get("name", "")).strip()
        if not name:
            raise MicrobiomeSuiteError("Each grid entry needs a non-empty 'name'.")
        if name in names:
            raise MicrobiomeSuiteError(f"Duplicate grid point name: {name}.")
        names.add(name)
        points.append(GridPoint(
            name=name, params=dict(entry.get("params", {})),
            is_baseline=bool(entry.get("baseline", False)),
        ))
    baselines = sum(p.is_baseline for p in points)
    if baselines != 1:
        raise MicrobiomeSuiteError(
            f"Grid config must flag exactly one baseline (found {baselines})."
        )
    return points


def _axis_point_name(params: dict) -> str:
    return "_".join(f"{k}{v}" for k, v in params.items()) or "point"


def _grid_from_axes(axes: dict[str, list]) -> list[GridPoint]:
    keys = [k for k in SWEEP_AXES if k in axes]
    unknown = [k for k in axes if k not in SWEEP_AXES]
    if unknown:
        raise MicrobiomeSuiteError(f"Unknown sweep axes: {unknown}. Allowed: {list(SWEEP_AXES)}.")
    if not keys:
        raise MicrobiomeSuiteError("No sweep axes provided.")
    combos = list(product(*(axes[k] for k in keys)))
    points: list[GridPoint] = []
    for index, combo in enumerate(combos):
        params = dict(zip(keys, combo))
        is_baseline = index == 0
        name = "baseline" if is_baseline else _axis_point_name(params)
        points.append(GridPoint(name=name, params=params, is_baseline=is_baseline))
    return points


def build_grid(
    *, config: Path | None = None, axes: dict[str, list] | None = None
) -> list[GridPoint]:
    has_config = config is not None
    has_axes = bool(axes)
    if has_config == has_axes:
        raise MicrobiomeSuiteError(
            "Provide exactly one grid source: a --grid-config file or axis values."
        )
    return _grid_from_config(config) if has_config else _grid_from_axes(axes)
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_dada2_sweep_grid.py -v`
Expected: PASS (4).

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/methods/dada2_sweep.py tests/test_dada2_sweep_grid.py
git commit -m "feat(dada2): sweep grid builder (json config + cli axes)"
```

---

### Task 3: orchestration + CLI

**Files:**
- Modify: `src/microsuite/methods/dada2_sweep.py`
- Modify: `src/microsuite/cli/method_features_cmd.py`
- Test: `tests/test_dada2_sweep_cli.py`

**Interfaces:**
- Consumes: `build_grid`, `GridPoint`, `SweepRun`, `summarize_sweep`, `write_sweep_summary`; `denoise` (from `methods/denoise.py`).
- Produces: `run_dada2_sweep(*, input_dir, mode, output_dir, grid, runtime="local", dada2_image=None, threads=1, force=False, timeout=None) -> Path`; CLI `denoise-sweep`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_dada2_sweep_cli.py
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.methods.dada2_sweep import GridPoint, run_dada2_sweep

STATS = (
    "\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\n"
    "s1\t1000\t900\t880\t870\t800\t700\n"
)


def _fake_denoise_factory(fail_names=()):
    def fake_denoise(*, backend, demux, output_table, output_rep_seqs, output_stats, **kw):
        name = Path(output_table).parent.name
        if name in fail_names:
            raise MicrobiomeSuiteError(f"denoise failed for {name}")
        Path(output_table).write_text("\ts1\nASV1\t5\n", encoding="utf-8")
        Path(output_rep_seqs).write_text(">ASV1\nAAA\n", encoding="utf-8")
        Path(output_stats).write_text(STATS, encoding="utf-8")
    return fake_denoise


def test_run_dada2_sweep_writes_summary(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr("microsuite.methods.dada2_sweep.denoise", _fake_denoise_factory())
    grid = [
        GridPoint(name="baseline", params={"max_ee_f": 2}, is_baseline=True),
        GridPoint(name="relaxed", params={"max_ee_f": 3}, is_baseline=False),
    ]
    out = run_dada2_sweep(
        input_dir=tmp_path / "reads", mode="paired", output_dir=tmp_path / "out", grid=grid
    )
    assert out.exists()
    df = pd.read_csv(out, sep="\t")
    assert list(df["name"]) == ["baseline", "relaxed"]
    assert set(df["status"]) == {"ok"}


def test_run_dada2_sweep_failed_point_recorded(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(
        "microsuite.methods.dada2_sweep.denoise", _fake_denoise_factory(fail_names={"relaxed"})
    )
    grid = [
        GridPoint(name="baseline", params={"max_ee_f": 2}, is_baseline=True),
        GridPoint(name="relaxed", params={"max_ee_f": 1}, is_baseline=False),
    ]
    with pytest.warns(UserWarning, match="relaxed"):
        out = run_dada2_sweep(
            input_dir=tmp_path / "reads", mode="paired", output_dir=tmp_path / "out", grid=grid
        )
    df = pd.read_csv(out, sep="\t")
    assert df[df["name"] == "relaxed"].iloc[0]["status"] == "failed"


def test_run_dada2_sweep_baseline_failure_raises(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(
        "microsuite.methods.dada2_sweep.denoise", _fake_denoise_factory(fail_names={"baseline"})
    )
    grid = [GridPoint(name="baseline", params={}, is_baseline=True)]
    with pytest.raises(MicrobiomeSuiteError, match="[Bb]aseline"):
        run_dada2_sweep(
            input_dir=tmp_path / "reads", mode="paired", output_dir=tmp_path / "out", grid=grid
        )


def test_cli_denoise_sweep_axes(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr("microsuite.methods.dada2_sweep.denoise", _fake_denoise_factory())
    (tmp_path / "reads").mkdir()
    r = CliRunner().invoke(app, [
        "denoise-sweep", "--input-dir", str(tmp_path / "reads"), "--mode", "paired",
        "--output-dir", str(tmp_path / "out"), "--max-ee-f", "2,3",
    ])
    assert r.exit_code == 0, r.stdout
    assert (tmp_path / "out" / "dada2_sweep_summary.tsv").exists()


def test_cli_denoise_sweep_both_sources_errors(tmp_path) -> None:
    cfg = tmp_path / "g.json"
    cfg.write_text("[]", encoding="utf-8")
    r = CliRunner().invoke(app, [
        "denoise-sweep", "--input-dir", str(tmp_path / "reads"), "--mode", "paired",
        "--output-dir", str(tmp_path / "out"), "--grid-config", str(cfg), "--max-ee-f", "2",
    ])
    assert r.exit_code != 0
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_dada2_sweep_cli.py -v`
Expected: FAIL (`cannot import name 'run_dada2_sweep'`).

- [ ] **Step 3: Add `run_dada2_sweep` to `methods/dada2_sweep.py`**

Add `import warnings` at the top and `from microsuite.methods.denoise import denoise` at module level (so tests can monkeypatch `microsuite.methods.dada2_sweep.denoise`).

```python
def run_dada2_sweep(
    *, input_dir: Path, mode: str, output_dir: Path, grid: list[GridPoint],
    runtime: str = "local", dada2_image: str | None = None, threads: int = 1,
    force: bool = False, timeout: float | None = None,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    ordered = sorted(grid, key=lambda p: not p.is_baseline)  # baseline first
    runs: list[SweepRun] = []
    for point in ordered:
        point_dir = output_dir / point.name
        point_dir.mkdir(parents=True, exist_ok=True)
        table = point_dir / "table.tsv"
        rep_seqs = point_dir / "rep_seqs.fasta"
        stats = point_dir / "stats.tsv"
        status = "ok"
        try:
            denoise(
                backend="dada2-r", demux=input_dir, mode=mode,
                output_table=table, output_rep_seqs=rep_seqs, output_stats=stats,
                runtime=runtime, dada2_image=dada2_image, threads=threads,
                force=force, timeout=timeout, **point.params,
            )
        except MicrobiomeSuiteError as exc:
            if point.is_baseline:
                raise MicrobiomeSuiteError(
                    f"Baseline sweep run '{point.name}' failed: {exc}"
                ) from exc
            warnings.warn(f"Sweep run '{point.name}' failed: {exc}", stacklevel=2)
            status = "failed"
        runs.append(SweepRun(point=point, table=table, rep_seqs=rep_seqs, stats=stats, status=status))
    summary = summarize_sweep(runs)
    return write_sweep_summary(summary, output_dir / "dada2_sweep_summary.tsv")
```

- [ ] **Step 4: Add the `denoise-sweep` CLI command**

In `cli/method_features_cmd.py`, inside `register(app)` (beside `denoise_cmd`), add — mirroring `denoise_cmd`'s option style. Add a module-level helper to parse a comma list, and build the axes dict from only the provided options:

```python
    @app.command("denoise-sweep")
    def denoise_sweep_cmd(
        input_dir: Annotated[Path, typer.Option("--input-dir", help="Directory of demultiplexed FASTQs.")],
        mode: Annotated[str, typer.Option("--mode", help="DADA2 mode: single or paired.")],
        output_dir: Annotated[Path, typer.Option("--output-dir", help="Sweep output directory.")],
        grid_config: Annotated[Path | None, typer.Option("--grid-config", help="JSON grid config.")] = None,
        max_ee_f: Annotated[str | None, typer.Option("--max-ee-f", help="Comma list, e.g. 2,3,5.")] = None,
        max_ee_r: Annotated[str | None, typer.Option("--max-ee-r")] = None,
        trunc_len_f: Annotated[str | None, typer.Option("--trunc-len-f")] = None,
        trunc_len_r: Annotated[str | None, typer.Option("--trunc-len-r")] = None,
        trunc_q: Annotated[str | None, typer.Option("--trunc-q")] = None,
        min_overlap: Annotated[str | None, typer.Option("--min-overlap")] = None,
        max_ee: Annotated[str | None, typer.Option("--max-ee")] = None,
        trunc_len: Annotated[str | None, typer.Option("--trunc-len")] = None,
        runtime: Annotated[str, typer.Option("--runtime", help="local or docker.")] = "local",
        image: Annotated[str | None, typer.Option("--image", help="Container image (dada2-r docker).")] = None,
        threads: Annotated[int, typer.Option("--threads")] = 1,
        force: Annotated[bool, typer.Option("--force")] = False,
    ) -> None:
        from microsuite.methods.dada2_sweep import build_grid, run_dada2_sweep

        float_axes = {"max_ee_f": max_ee_f, "max_ee_r": max_ee_r, "max_ee": max_ee}
        int_axes = {
            "trunc_len_f": trunc_len_f, "trunc_len_r": trunc_len_r,
            "trunc_q": trunc_q, "min_overlap": min_overlap, "trunc_len": trunc_len,
        }
        axes: dict[str, list] = {}
        for key, raw in float_axes.items():
            if raw is not None:
                axes[key] = [float(v) for v in raw.split(",")]
        for key, raw in int_axes.items():
            if raw is not None:
                axes[key] = [int(v) for v in raw.split(",")]

        grid = build_grid(config=grid_config, axes=axes or None)
        run_dada2_sweep(
            input_dir=input_dir, mode=mode, output_dir=output_dir, grid=grid,
            runtime=runtime, dada2_image=image, threads=threads, force=force,
        )
```

(`Annotated`, `typer`, `Path` are already imported in `method_features_cmd.py`.)

- [ ] **Step 5: Run to verify pass + full suite + lint**

Run: `uv run pytest tests/test_dada2_sweep_cli.py -v` (5 pass), then `uv run pytest -q` (all green), then `uv run ruff check .` and `uv run ruff format --check .` (both clean; `uv run ruff format .` on the two changed files if needed, then re-check). Sanity: `uv run python -c "from typer.testing import CliRunner; from microsuite.cli.app import app; print(CliRunner().invoke(app, ['denoise-sweep','--help']).stdout)"` shows `--grid-config` and `--max-ee-f`.

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/methods/dada2_sweep.py src/microsuite/cli/method_features_cmd.py tests/test_dada2_sweep_cli.py
git commit -m "feat(dada2): denoise-sweep orchestration + CLI"
```

---

## Self-Review

**Spec coverage:**
- Grid via JSON config OR CLI axes, one baseline (axes → point 0) → Task 2 `build_grid`. ✓
- Per-run retention/ASV/depth/observed/chimera → Task 1 `run_metrics`. ✓
- Baseline similarity: shared-ASV (by sequence) abundance corr + per-sample corr → Task 1 `compare_to_baseline`. ✓
- Baseline-first `dada2_sweep_summary.tsv`, failed non-baseline recorded, baseline failure fatal → Task 1 `summarize_sweep` + Task 3 `run_dada2_sweep`. ✓
- `denoise-sweep` CLI, both/neither source errors → Task 3. ✓
- Metrics/grid offline-tested; orchestration with stubbed `denoise` → Tasks 1-3 tests. ✓
- Both CI gates → Task 3 Step 5. ✓

**Placeholder scan:** none — full module code across the three tasks, exact CLI additions, and complete tests with concrete fixtures and asserted numbers.

**Consistency:** `GridPoint`/`SweepRun`/`build_grid`/`run_metrics`/`compare_to_baseline`/`summarize_sweep`/`run_dada2_sweep` names and signatures match across the module, its tests, and the CLI; the ASV table on-disk contract (`index_col=0`, features×samples) is used identically by the writer fixtures and `_read_asv_table`; `denoise` is imported at module level so the tests' `monkeypatch.setattr("microsuite.methods.dada2_sweep.denoise", …)` targets the same name the orchestration calls; the CLI axis names match `SWEEP_AXES`.
