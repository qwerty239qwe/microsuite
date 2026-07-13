from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity.beta import beta_diversity
from microsuite.methods.abundance import abundance_native
from microsuite.methods.normalize import normalize_native
from microsuite.ordination.pcoa import pcoa

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def _require_obs_column(adata: ad.AnnData, column: str) -> pd.Series:
    if column not in adata.obs.columns:
        available = ", ".join(str(c) for c in adata.obs.columns) or "(none)"
        raise MicrobiomeSuiteError(
            f"Metadata column '{column}' not found in obs; available: {available}."
        )
    return pd.Series(adata.obs[column])


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


def plot_taxa_by_group(
    adata: ad.AnnData,
    *,
    level: str,
    group_by: str,
    output: Path,
    top_n: int = 15,
    group_order: list[str] | None = None,
) -> None:
    group = _require_obs_column(adata, group_by)
    frame = abundance_native(adata, level=level, relative=True)  # samples x taxa (raises on level)
    grouped = frame.groupby(group.astype(str).to_numpy()).mean()  # group x taxa
    grouped = grouped.reindex(_resolve_group_order(grouped.index, group_order))
    order = grouped.mean(axis=0).sort_values(ascending=False)
    top = list(order.head(top_n).index)
    plot_df = grouped[top].copy()
    other = grouped.drop(columns=top).sum(axis=1)
    if (other > 0).any():
        plot_df["Other"] = other
    taxa_cols = [c for c in plot_df.columns if c != "Other"]
    color_map = dict(zip(taxa_cols, _qualitative_colors(len(taxa_cols)), strict=True))
    color_map["Other"] = "0.7"
    colors = [color_map[c] for c in plot_df.columns]
    width = max(6.0, len(plot_df.index) * 0.9)
    height = max(4.5, min(12.0, 2.5 + plot_df.shape[1] * 0.25))
    ax = plot_df.plot(kind="bar", stacked=True, figsize=(width, height), width=0.85, color=colors)
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
    group_order: list[str] | None = None,
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
    categories = _resolve_group_order(pd.unique(groups), group_order)

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
        colors = _qualitative_colors(len(categories))
        handles = []
        for gi, cat in enumerate(categories):
            series = [clr_top.loc[groups == cat, taxon].to_numpy() for taxon in top]
            series = [d if len(d) else np.array([np.nan]) for d in series]
            positions = np.arange(len(top)) + (gi - (n_groups - 1) / 2) * slot
            color = colors[gi]
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
                for body in vp["bodies"]:  # ty: ignore[not-iterable]
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


_ORDINATION_STYLES = ("scatter", "trajectory", "facet")


def _pc_coords(adata: ad.AnnData) -> pd.DataFrame:
    dist = beta_diversity(adata, "bray-curtis")
    coords = pcoa(dist, dimensions=2).set_index("sample_id")
    return coords.loc[[str(s) for s in adata.obs_names]]


def plot_braycurtis_ordination(
    adata: ad.AnnData,
    *,
    color_by: str,
    output: Path,
    subject: str | None = None,
    style: str | None = None,
    order_by: str | None = None,
    order: list[str] | None = None,
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

    order_column = _require_obs_column(adata, order_by) if order_by is not None else color
    order_vals = np.asarray([str(v) for v in order_column.to_numpy()])
    order_ranking = _resolve_group_order(pd.unique(order_vals), order)
    order_rank = {value: index for index, value in enumerate(order_ranking)}

    coords = _pc_coords(adata)
    x = coords["PC1"].to_numpy()
    y = coords["PC2"].to_numpy()
    xlab = f"PC1 ({coords['PC1_variance'].iloc[0] * 100:.1f}%)"
    ylab = f"PC2 ({coords['PC2_variance'].iloc[0] * 100:.1f}%)"
    color_vals = color.to_numpy()
    is_numeric = pd.api.types.is_numeric_dtype(color)

    if effective == "facet":
        assert subj is not None  # guaranteed above for facet/trajectory styles
        subjects = list(pd.unique(subj.to_numpy()))
        ncols = min(3, len(subjects))
        nrows = int(np.ceil(len(subjects) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.5 * nrows), squeeze=False)
        subj_vals = subj.to_numpy()

        if is_numeric:
            values = color_vals.astype(float)
            vmin = float(np.nanmin(values))
            vmax = float(np.nanmax(values))
            if vmin == vmax:
                vmax = vmin + 1
            cmap = plt.get_cmap("viridis")
            mappable = None
            for idx, sub in enumerate(subjects):
                ax = axes[idx // ncols][idx % ncols]
                mask = subj_vals == sub
                mappable = ax.scatter(
                    x[mask], y[mask], c=values[mask], cmap=cmap, vmin=vmin, vmax=vmax, s=40
                )
                ax.set_title(str(sub))
                ax.set_xlabel(xlab)
                ax.set_ylabel(ylab)
            for j in range(len(subjects), nrows * ncols):
                axes[j // ncols][j % ncols].axis("off")
            fig.tight_layout()
            if mappable is not None:
                fig.colorbar(mappable, ax=axes, label=color_by)
        else:
            categories = list(pd.unique(color_vals))
            cmap = plt.get_cmap("tab10")
            cat_color = {c: cmap(i % 10) for i, c in enumerate(categories)}
            for idx, sub in enumerate(subjects):
                ax = axes[idx // ncols][idx % ncols]
                mask = subj_vals == sub
                colors = [cat_color[v] for v in color_vals[mask]]
                ax.scatter(x[mask], y[mask], color=colors, s=40)
                ax.set_title(str(sub))
                ax.set_xlabel(xlab)
                ax.set_ylabel(ylab)
            for j in range(len(subjects), nrows * ncols):
                axes[j // ncols][j % ncols].axis("off")
            fig.tight_layout()
            handles = [
                plt.Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor=cat_color[cat],
                    label=str(cat),
                    markersize=8,
                )
                for cat in categories
            ]
            fig.legend(handles=handles, title=color_by, bbox_to_anchor=(1.02, 1), loc="upper left")

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
        assert subj is not None  # guaranteed above for facet/trajectory styles
        subj_vals = subj.to_numpy()
        for sub in pd.unique(subj_vals):
            idx = np.where(subj_vals == sub)[0]
            idx = sorted(idx, key=lambda i: order_rank.get(order_vals[i], len(order_rank)))
            ax.plot(x[idx], y[idx], color="gray", alpha=0.6, lw=1.0, zorder=0)

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(f"Bray-Curtis PCoA by {color_by}")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)
