from __future__ import annotations

from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.diffab.ancombc import run_ancombc
from microsuite.io.h5ad import read_h5ad
from microsuite.methods._qiime import require_qiime, run_qiime

SUPPORTED_BACKENDS = ("ancombc", "qiime2-ancombc", "aldex2", "maaslin2", "lefse")
PLANNED_BACKENDS = ("aldex2", "maaslin2", "lefse")


def diff_abundance(
    *,
    backend: str,
    table: Path,
    group: str,
    output: Path,
    metadata: Path | None = None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = backend.lower()
    if backend in PLANNED_BACKENDS:
        raise MicrobiomeSuiteError(
            f"Differential abundance backend '{backend}' is registered but not implemented yet. "
            "Use --backend ancombc for now."
        )
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
    if backend != "ancombc":
        raise MicrobiomeSuiteError(
            f"Unsupported differential abundance backend '{backend}'. "
            f"Choose one of: {', '.join(SUPPORTED_BACKENDS)}"
        )

    adata = read_h5ad(ensure_input(table))
    run_ancombc(
        adata,
        group=group,
        output=prepare_output(output, force=force),
        run_dir=run_dir,
        timeout=timeout,
    )
