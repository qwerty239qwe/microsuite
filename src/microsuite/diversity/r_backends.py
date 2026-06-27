from __future__ import annotations

import shutil
from collections.abc import Iterable
from importlib.resources import files
from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity._matrix import dense_counts
from microsuite.runtime.runner import CommandLog, run_command

R_ALPHA_BACKENDS = ("breakaway", "inext")
R_ALPHA_PACKAGES = {"breakaway": "breakaway", "inext": "iNEXT"}


def r_alpha_diversity(
    adata: ad.AnnData,
    *,
    backend: str,
    q: Iterable[float] = (0.0, 1.0, 2.0),
    datatype: str = "abundance",
    knots: int = 40,
    se: bool = True,
    conf: float = 0.95,
    nboot: int = 50,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> pd.DataFrame:
    if backend not in R_ALPHA_BACKENDS:
        raise MicrobiomeSuiteError(f"Unsupported R alpha-diversity backend: {backend}")
    rscript = shutil.which("Rscript")
    if rscript is None:
        package = R_ALPHA_PACKAGES[backend]
        raise MicrobiomeSuiteError(
            f"{backend} alpha diversity requires external Rscript and the R package "
            f"'{package}'. Install R, install {package}, then rerun this command."
        )
    if datatype not in {"abundance", "incidence_freq", "incidence_raw"}:
        raise MicrobiomeSuiteError(
            "iNEXT datatype must be abundance, incidence_freq, or incidence_raw."
        )
    if knots < 1:
        raise MicrobiomeSuiteError("iNEXT knots must be at least 1.")
    if not 0 < conf < 1:
        raise MicrobiomeSuiteError("iNEXT conf must be greater than 0 and less than 1.")
    if nboot < 0:
        raise MicrobiomeSuiteError("iNEXT nboot must be non-negative.")

    with TemporaryDirectory(prefix=f"microsuite-{backend}-") as temp_dir:
        temp = Path(temp_dir)
        counts_path = temp / "counts.tsv"
        output_path = temp / "alpha.tsv"
        pd.DataFrame(dense_counts(adata).T, index=adata.var_names, columns=adata.obs_names).to_csv(
            counts_path, sep="\t"
        )

        script = files("microsuite.diversity.r").joinpath(f"{backend}_alpha.R")
        command = [rscript, str(script), str(counts_path), str(output_path)]
        if backend == "inext":
            command.extend(
                [
                    _format_q(q),
                    datatype,
                    str(knots),
                    str(se).lower(),
                    str(conf),
                    str(nboot),
                ]
            )

        run_command(
            command,
            f"{backend} alpha diversity failed.",
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(
                task="alpha_diversity",
                backend=backend,
                inputs={"table_shape": list(adata.shape)},
                outputs={"output": str(output_path)},
                params={
                    "q": list(q),
                    "datatype": datatype,
                    "knots": knots,
                    "se": se,
                    "conf": conf,
                    "nboot": nboot,
                }
                if backend == "inext"
                else None,
            ),
        )
        return pd.read_csv(output_path, sep="\t")


def _format_q(q: Iterable[float]) -> str:
    values = [float(value) for value in q]
    if not values:
        raise MicrobiomeSuiteError("iNEXT q must contain at least one diversity order.")
    return ",".join(f"{value:g}" for value in values)
