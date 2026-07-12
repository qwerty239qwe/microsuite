from __future__ import annotations

import math
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
        pooled["assigned_features"].to_numpy(dtype=float),
        total,
        out=np.zeros_like(total),
        where=total > 0,
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
    import pandas as pd

    summary = summarize_assignment(adata)
    per_sample = summary[summary["sample"] != POOLED_LABEL]
    pivot = per_sample.pivot(index="sample", columns="rank", values="assigned_read_frac")
    ranks = [r for r in LEVELS if r in pivot.columns]
    pivot = pivot[ranks]
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
