from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.abundance import abundance_native
from microsuite.methods.normalize import normalize_native

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


_CLR_STYLES = ("boxplot", "heatmap", "violin")


def plot_clr_by_group(
    adata: ad.AnnData,
    *,
    level: str,
    group_by: str,
    output: Path,
    top_n: int = 10,
    style: str = "boxplot",
) -> None:
    if style not in _CLR_STYLES:
        raise MicrobiomeSuiteError(
            f"Unknown style '{style}'. Choose one of: {', '.join(_CLR_STYLES)}."
        )
    group = _require_obs_column(adata, group_by).astype(str)
    # samples x taxa (raises on level)
    counts = abundance_native(adata, level=level, relative=False)
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
                    series,
                    positions=positions,
                    widths=slot * 0.9,
                    patch_artist=True,
                    manage_ticks=False,
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
