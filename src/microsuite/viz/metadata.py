from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
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
