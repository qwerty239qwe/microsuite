from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.io.h5ad import read_h5ad
from microsuite.methods._dispatch import require_backend
from microsuite.methods.abundance import abundance_native

SUPPORTED_BACKENDS = ("native",)


def shared_taxa(
    *,
    backend: str,
    table: Path,
    output: Path,
    level: str,
    group: str,
    force: bool = False,
) -> None:
    backend = require_backend(backend, SUPPORTED_BACKENDS, "shared_taxa")
    result = shared_taxa_native(read_h5ad(ensure_input(table)), level=level, group=group)
    result.to_csv(prepare_output(output, force=force), sep="\t", index=False)


def shared_taxa_native(adata: ad.AnnData, *, level: str, group: str) -> pd.DataFrame:
    obs = pd.DataFrame(adata.obs)
    if group not in obs.columns:
        raise MicrobiomeSuiteError(f"Sample metadata group not found: {group}")

    counts = abundance_native(adata, level=level, relative=False)
    presence = counts > 0
    groups = obs[group].astype(str)
    group_presence = presence.groupby(groups).any()
    rows = []
    for taxon in group_presence.columns:
        present_groups = group_presence.index[group_presence[taxon]].astype(str).tolist()
        rows.append(
            {
                "taxon": taxon,
                "n_groups": len(present_groups),
                "groups": ",".join(present_groups),
            }
        )
    return pd.DataFrame(rows).sort_values(["n_groups", "taxon"], ascending=[False, True])
