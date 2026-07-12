from __future__ import annotations

import json
from dataclasses import dataclass, field
from itertools import product
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


SWEEP_AXES = (
    "max_ee_f",
    "max_ee_r",
    "trunc_len_f",
    "trunc_len_r",
    "trunc_q",
    "min_overlap",
    "max_ee",
    "trunc_len",
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
        if not isinstance(entry, dict):
            raise MicrobiomeSuiteError("Each grid entry must be a JSON object.")
        name = str(entry.get("name", "")).strip()
        if not name:
            raise MicrobiomeSuiteError("Each grid entry needs a non-empty 'name'.")
        if name in names:
            raise MicrobiomeSuiteError(f"Duplicate grid point name: {name}.")
        names.add(name)
        points.append(
            GridPoint(
                name=name,
                params=dict(entry.get("params", {})),
                is_baseline=bool(entry.get("baseline", False)),
            )
        )
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
        params = dict(zip(keys, combo, strict=True))
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
    *,
    table_path: Path,
    rep_seqs_path: Path,
    baseline_table_path: Path,
    baseline_rep_seqs_path: Path,
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
        base_full[samples].sum(axis=0).to_numpy(),
        method="pearson",
    )
    obs_p = _corr(
        (run_full[samples] > 0).sum(axis=0).to_numpy(),
        (base_full[samples] > 0).sum(axis=0).to_numpy(),
        method="pearson",
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
        row: dict = {
            "name": run.point.name,
            "is_baseline": run.point.is_baseline,
            "status": run.status,
        }
        row.update({f"param_{k}": v for k, v in run.point.params.items()})
        if run.status == "ok":
            row.update(run_metrics(run.table, run.stats))
            if baseline.status == "ok":
                row.update(
                    compare_to_baseline(
                        table_path=run.table,
                        rep_seqs_path=run.rep_seqs,
                        baseline_table_path=baseline.table,
                        baseline_rep_seqs_path=baseline.rep_seqs,
                    )
                )
        rows.append(row)
    return pd.DataFrame(rows)


def write_sweep_summary(summary: pd.DataFrame, out_path: Path) -> Path:
    summary.to_csv(out_path, sep="\t", index=False)
    return out_path
