from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity._matrix import dense_counts

ANCOMBC_SCRIPT = Path(__file__).resolve().parents[3] / "scripts" / "r" / "ancombc.R"


def run_ancombc(adata: ad.AnnData, *, group: str, output: Path) -> None:
    if group not in adata.obs.columns:
        raise MicrobiomeSuiteError(f"Group column not found in sample metadata: {group}")
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise MicrobiomeSuiteError(
            "ANCOM-BC requires external Rscript and the R package 'ANCOMBC'. "
            "Install R, install ANCOMBC, then rerun this command."
        )

    with TemporaryDirectory() as temp_dir:
        temp = Path(temp_dir)
        counts_path = temp / "counts.tsv"
        metadata_path = temp / "metadata.tsv"

        pd.DataFrame(dense_counts(adata).T, index=adata.var_names, columns=adata.obs_names).to_csv(
            counts_path, sep="\t"
        )
        pd.DataFrame(adata.obs).to_csv(metadata_path, sep="\t")

        result = subprocess.run(
            [
                rscript,
                str(ANCOMBC_SCRIPT),
                str(counts_path),
                str(metadata_path),
                group,
                str(output),
            ],
            check=False,
            text=True,
            capture_output=True,
        )
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip() or "ANCOM-BC failed."
        raise MicrobiomeSuiteError(message)
