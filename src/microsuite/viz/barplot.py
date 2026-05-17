from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity._matrix import dense_counts

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def taxonomy_barplot(adata: ad.AnnData, *, level: str, output: Path, top_n: int = 20) -> None:
    if level not in adata.var.columns:
        raise MicrobiomeSuiteError(f"Taxonomy level not found in table: {level}")

    counts = dense_counts(adata)
    taxa = adata.var[level].astype(str).replace("", "Unclassified")
    frame = pd.DataFrame(counts, index=adata.obs_names.astype(str), columns=taxa)
    grouped = frame.T.groupby(level=0).sum().T

    totals = grouped.sum(axis=1)
    relative = grouped.div(totals.replace(0, np.nan), axis=0).fillna(0.0)
    top_taxa = relative.sum(axis=0).sort_values(ascending=False).head(top_n).index
    plot_data = relative.loc[:, top_taxa].copy()
    other = relative.drop(columns=top_taxa).sum(axis=1)
    if (other > 0).any():
        plot_data["Other"] = other

    width = max(8.0, len(plot_data.index) * 0.45)
    height = max(4.5, min(12.0, 2.5 + plot_data.shape[1] * 0.25))
    ax = plot_data.plot(kind="bar", stacked=True, figsize=(width, height), width=0.9)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Relative abundance")
    ax.set_ylim(0, 1)
    ax.legend(title=level, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(output, dpi=160)
    plt.close()
