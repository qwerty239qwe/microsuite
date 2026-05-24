from __future__ import annotations

import shutil
from importlib.resources import files
from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity._matrix import dense_counts
from microsuite.runtime.runner import CommandLog, run_command

R_BACKEND_PACKAGES = {
    "aldex2": "ALDEx2",
    "maaslin2": "Maaslin2",
    "lefse": "lefser",
}


def run_r_diffab_backend(
    adata: ad.AnnData,
    *,
    backend: str,
    group: str,
    output: Path,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    if backend not in R_BACKEND_PACKAGES:
        raise MicrobiomeSuiteError(f"Unsupported R differential abundance backend: {backend}")
    if group not in adata.obs.columns:
        raise MicrobiomeSuiteError(f"Group column not found in sample metadata: {group}")
    rscript = shutil.which("Rscript")
    if rscript is None:
        package = R_BACKEND_PACKAGES[backend]
        raise MicrobiomeSuiteError(
            f"{backend} requires external Rscript and the R package '{package}'. "
            f"Install R, install {package}, then rerun this command."
        )

    script = files("microsuite.diffab.r").joinpath(f"{backend}.R")
    with TemporaryDirectory() as temp_dir:
        temp = Path(temp_dir)
        counts_path = temp / "counts.tsv"
        metadata_path = temp / "metadata.tsv"

        pd.DataFrame(dense_counts(adata).T, index=adata.var_names, columns=adata.obs_names).to_csv(
            counts_path, sep="\t"
        )
        pd.DataFrame(adata.obs).to_csv(metadata_path, sep="\t")

        run_command(
            [
                rscript,
                str(script),
                str(counts_path),
                str(metadata_path),
                group,
                str(output),
            ],
            failure_message=f"{backend} failed.",
            run_dir=run_dir,
            log=CommandLog(
                task="diff_abundance",
                backend=backend,
                inputs={"group": group},
                outputs={"output": str(output)},
            ),
            timeout=timeout,
        )
