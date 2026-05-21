from __future__ import annotations

from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods._qiime import ensure_inputs, prepare_outputs, require_qiime, run_qiime

SUPPORTED_BACKENDS = ("qiime2-taxonomy",)


def evaluate(
    *,
    backend: str,
    expected_taxa: Path,
    observed_taxa: Path,
    output: Path,
    feature_table: Path | None = None,
    depth: int = 7,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = backend.lower()
    if backend != "qiime2-taxonomy":
        backends = ", ".join(SUPPORTED_BACKENDS)
        raise MicrobiomeSuiteError(
            f"Unsupported evaluate backend '{backend}'. Choose one of: {backends}"
        )
    evaluate_qiime2_taxonomy(
        expected_taxa=expected_taxa,
        observed_taxa=observed_taxa,
        feature_table=feature_table,
        output=output,
        depth=depth,
        force=force,
        run_dir=run_dir,
        timeout=timeout,
    )


def evaluate_qiime2_taxonomy(
    *,
    expected_taxa: Path,
    observed_taxa: Path,
    feature_table: Path | None,
    output: Path,
    depth: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    qiime = require_qiime("QIIME 2 quality-control evaluate-taxonomy")
    ensure_inputs(expected_taxa, observed_taxa, feature_table)
    prepare_outputs(output, force=force)

    command = [
        qiime,
        "quality-control",
        "evaluate-taxonomy",
        "--i-expected-taxa",
        str(expected_taxa),
        "--i-observed-taxa",
        str(observed_taxa),
    ]
    if feature_table is not None:
        command.extend(["--i-feature-table", str(feature_table)])
    command.extend(
        [
            "--p-depth",
            str(depth),
            "--o-visualization",
            str(output),
        ]
    )
    run_qiime(
        command,
        "QIIME 2 quality-control evaluate-taxonomy failed.",
        run_dir=run_dir,
        timeout=timeout,
        task="evaluate",
        backend="qiime2-taxonomy",
    )
