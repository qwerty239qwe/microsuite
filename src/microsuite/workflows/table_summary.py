from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from microsuite import __version__
from microsuite._paths import ensure_input, prepare_output
from microsuite.diversity.alpha import alpha_diversity
from microsuite.diversity.beta import beta_diversity
from microsuite.io.h5ad import write_h5ad
from microsuite.io.qza import read_qza
from microsuite.io.tsv import read_tsv
from microsuite.ordination.pcoa import pcoa
from microsuite.viz.barplot import taxonomy_barplot


def run_table_summary(
    *,
    output: Path,
    table: Path,
    metadata: Path,
    taxonomy: Path | None = None,
    taxonomy_artifact: Path | None = None,
    input_format: str = "tsv",
    alpha_metric: str = "shannon",
    beta_metric: str = "bray-curtis",
    barplot_level: str = "genus",
    force: bool = False,
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    if input_format == "tsv":
        adata = read_tsv(
            ensure_input(table),
            ensure_input(metadata),
            ensure_input(taxonomy) if taxonomy else None,
        )
    elif input_format == "qza":
        adata = read_qza(
            ensure_input(table),
            ensure_input(metadata),
            ensure_input(taxonomy_artifact) if taxonomy_artifact else None,
        )
    else:
        raise ValueError(f"Unsupported workflow input format: {input_format}")

    table_h5ad = prepare_output(output / "table.h5ad", force=force)
    write_h5ad(adata, table_h5ad)

    alpha_path = prepare_output(output / f"alpha-{alpha_metric}.tsv", force=force)
    alpha_diversity(adata, alpha_metric).to_csv(alpha_path, sep="\t", index=False)

    beta_path = prepare_output(output / f"beta-{beta_metric}.tsv", force=force)
    beta = beta_diversity(adata, beta_metric)
    beta.to_csv(beta_path, sep="\t")

    pcoa_path = prepare_output(output / "pcoa.tsv", force=force)
    pcoa(beta).to_csv(pcoa_path, sep="\t", index=False)

    barplot_path = prepare_output(output / f"barplot-{barplot_level}.png", force=force)
    taxonomy_barplot(adata, level=barplot_level, output=barplot_path)

    run = {
        "toolbox": "microsuite",
        "version": __version__,
        "workflow": "table-summary",
        "created_at": datetime.now(UTC).isoformat(),
        "inputs": {
            "format": input_format,
            "table": str(table),
            "metadata": str(metadata),
            "taxonomy": str(taxonomy) if taxonomy else None,
            "taxonomy_artifact": str(taxonomy_artifact) if taxonomy_artifact else None,
        },
        "parameters": {
            "alpha_metric": alpha_metric,
            "beta_metric": beta_metric,
            "barplot_level": barplot_level,
        },
        "outputs": {
            "table_h5ad": str(table_h5ad),
            "alpha": str(alpha_path),
            "beta": str(beta_path),
            "pcoa": str(pcoa_path),
            "barplot": str(barplot_path),
        },
    }
    pd.Series(run["outputs"]).to_json(output / "outputs.json", indent=2)
    (output / "run.json").write_text(json.dumps(run, indent=2), encoding="utf-8")
