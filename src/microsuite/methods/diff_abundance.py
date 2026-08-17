from __future__ import annotations

from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.batch.value_type import require_value_types
from microsuite.diffab.ancombc import run_ancombc
from microsuite.diffab.maaslin3 import run_maaslin3
from microsuite.diffab.r_backends import run_r_diffab_backend
from microsuite.io.h5ad import read_h5ad
from microsuite.methods._dispatch import require_backend
from microsuite.methods._qiime import require_qiime, run_qiime

SUPPORTED_BACKENDS = ("ancombc", "qiime2-ancombc", "aldex2", "maaslin2", "maaslin3", "lefse")
R_BACKENDS = ("aldex2", "maaslin2", "lefse")
COUNT_REQUIRING_BACKENDS = ("ancombc", "aldex2")


def diff_abundance(
    *,
    backend: str,
    table: Path,
    group: str | None = None,
    output: Path,
    metadata: Path | None = None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
    runtime: str = "local",
    image: str | None = None,
    engine: str = "docker",
    formula: str | None = None,
    fix_formula: str | None = None,
    rand_formula: str | None = None,
    normalization: str = "TSS",
    transform: str = "LOG",
    min_prevalence: float = 0.0,
    min_abundance: float = 0.0,
) -> None:
    backend = require_backend(backend, SUPPORTED_BACKENDS, "differential abundance")
    if backend == "maaslin3":
        adata = read_h5ad(ensure_input(table))
        run_maaslin3(
            adata,
            output=output,
            group=group,
            formula=formula,
            fix_formula=fix_formula,
            rand_formula=rand_formula,
            normalization=normalization,
            transform=transform,
            min_prevalence=min_prevalence,
            min_abundance=min_abundance,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
            runtime=runtime,
            image=image,
            engine=engine,
        )
        return
    if any(value is not None for value in (formula, fix_formula, rand_formula)):
        raise MicrobiomeSuiteError(
            "--formula, --fix-formula, and --rand-formula are currently supported only "
            "by --backend maaslin3 on the unified diff_abundance command."
        )
    if (normalization, transform, min_prevalence, min_abundance) != ("TSS", "LOG", 0.0, 0.0):
        raise MicrobiomeSuiteError(
            "--normalization, --transform, --min-prevalence, and --min-abundance are "
            "currently supported only by --backend maaslin3."
        )
    if group is None:
        raise MicrobiomeSuiteError(f"--group is required for --backend {backend}.")
    if backend == "qiime2-ancombc":
        if metadata is None:
            raise MicrobiomeSuiteError("--metadata is required for --backend qiime2-ancombc.")
        qiime = require_qiime("QIIME 2 composition ANCOM-BC")
        ensure_input(table)
        ensure_input(metadata)
        prepare_output(output, force=force)
        run_qiime(
            [
                qiime,
                "composition",
                "ancombc",
                "--i-table",
                str(table),
                "--m-metadata-file",
                str(metadata),
                "--p-formula",
                group,
                "--o-differentials",
                str(output),
            ],
            "QIIME 2 ANCOM-BC failed.",
            run_dir=run_dir,
            timeout=timeout,
            task="diff_abundance",
            backend=backend,
        )
        return
    if backend in R_BACKENDS:
        adata = read_h5ad(ensure_input(table))
        if backend in COUNT_REQUIRING_BACKENDS:
            require_value_types(adata, ("counts",), operation=f"diff_abundance --backend {backend}")
        run_r_diffab_backend(
            adata,
            backend=backend,
            group=group,
            output=prepare_output(output, force=force),
            run_dir=run_dir,
            timeout=timeout,
            runtime=runtime,
            image=image,
            engine=engine,
        )
        return

    adata = read_h5ad(ensure_input(table))
    require_value_types(adata, ("counts",), operation=f"diff_abundance --backend {backend}")
    run_ancombc(
        adata,
        group=group,
        output=prepare_output(output, force=force),
        run_dir=run_dir,
        timeout=timeout,
        runtime=runtime,
        image=image,
        engine=engine,
    )
