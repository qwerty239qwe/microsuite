from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diffab._runner import invoke_r_backend
from microsuite.diversity._matrix import dense_counts
from microsuite.runtime.runner import CommandLog

R_BACKEND_PACKAGES = {
    "aldex2": "ALDEx2",
    "maaslin2": "Maaslin2",
}


def run_r_diffab_backend(
    adata: ad.AnnData,
    *,
    backend: str,
    group: str,
    output: Path,
    run_dir: Path | None = None,
    timeout: float | None = None,
    runtime: str = "local",
    image: str | None = None,
    engine: str = "docker",
) -> None:
    if backend not in R_BACKEND_PACKAGES:
        raise MicrobiomeSuiteError(f"Unsupported R differential abundance backend: {backend}")
    if group not in adata.obs.columns:
        raise MicrobiomeSuiteError(f"Group column not found in sample metadata: {group}")

    package = R_BACKEND_PACKAGES[backend]
    with TemporaryDirectory() as temp_dir:
        temp = Path(temp_dir)
        counts_path = temp / "counts.tsv"
        metadata_path = temp / "metadata.tsv"

        pd.DataFrame(dense_counts(adata).T, index=adata.var_names, columns=adata.obs_names).to_csv(
            counts_path, sep="\t"
        )
        pd.DataFrame(adata.obs).to_csv(metadata_path, sep="\t")

        invoke_r_backend(
            backend=backend,
            positional=[counts_path, metadata_path, group, output],
            runtime=runtime,
            image=image,
            engine=engine,
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(
                task="diff_abundance",
                backend=backend,
                inputs={"group": group},
                outputs={"output": str(output)},
            ),
            local_missing_message=(
                f"{backend} requires external Rscript and the R package '{package}'. "
                f"Install R, install {package}, then rerun this command, or use "
                f"--runtime docker with the r-diffab-{backend} image."
            ),
        )
