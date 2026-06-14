from __future__ import annotations

from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods._dispatch import require_backend
from microsuite.methods._qiime import ensure_inputs, prepare_outputs, require_qiime, run_qiime

SUPPORTED_BACKENDS = ("qiime2-decontam",)


def decontam(
    *,
    backend: str,
    table: Path,
    metadata: Path,
    output: Path,
    method: str = "prevalence",
    prev_control_column: str | None = None,
    prev_control_indicator: str | None = None,
    freq_concentration_column: str | None = None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = require_backend(backend, SUPPORTED_BACKENDS, "decontam")
    decontam_qiime2(
        table=table,
        metadata=metadata,
        output=output,
        method=method,
        prev_control_column=prev_control_column,
        prev_control_indicator=prev_control_indicator,
        freq_concentration_column=freq_concentration_column,
        force=force,
        run_dir=run_dir,
        timeout=timeout,
    )


def decontam_qiime2(
    *,
    table: Path,
    metadata: Path,
    output: Path,
    method: str,
    prev_control_column: str | None,
    prev_control_indicator: str | None,
    freq_concentration_column: str | None,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if method in {"prevalence", "combined"} and (
        prev_control_column is None or prev_control_indicator is None
    ):
        raise MicrobiomeSuiteError(
            "--prev-control-column and --prev-control-indicator are required for "
            "prevalence or combined decontam methods."
        )
    if method in {"frequency", "combined"} and freq_concentration_column is None:
        raise MicrobiomeSuiteError(
            "--freq-concentration-column is required for frequency or combined decontam methods."
        )

    qiime = require_qiime("QIIME 2 quality-control decontam-identify")
    ensure_inputs(table, metadata)
    prepare_outputs(output, force=force)

    command = [
        qiime,
        "quality-control",
        "decontam-identify",
        "--i-table",
        str(table),
        "--m-metadata-file",
        str(metadata),
        "--p-method",
        method,
    ]
    if prev_control_column is not None:
        command.extend(["--p-prev-control-column", prev_control_column])
    if prev_control_indicator is not None:
        command.extend(["--p-prev-control-indicator", prev_control_indicator])
    if freq_concentration_column is not None:
        command.extend(["--p-freq-concentration-column", freq_concentration_column])
    command.extend(["--o-decontam-scores", str(output)])
    run_qiime(
        command,
        "QIIME 2 quality-control decontam-identify failed.",
        run_dir=run_dir,
        timeout=timeout,
        task="decontam",
        backend="qiime2-decontam",
    )
