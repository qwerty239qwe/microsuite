from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output


@dataclass(frozen=True)
class QIIME2DiversityMetric:
    action: str
    output_option: str
    thread_option: str | None
    requires_phylogeny: bool


QIIME2_DIVERSITY_METRICS = {
    "observed_features": QIIME2DiversityMetric("observed-features", "--o-vector", None, False),
    "observed-features": QIIME2DiversityMetric("observed-features", "--o-vector", None, False),
    "shannon": QIIME2DiversityMetric("shannon-entropy", "--o-vector", None, False),
    "shannon_entropy": QIIME2DiversityMetric("shannon-entropy", "--o-vector", None, False),
    "shannon-entropy": QIIME2DiversityMetric("shannon-entropy", "--o-vector", None, False),
    "pielou": QIIME2DiversityMetric("pielou-evenness", "--o-vector", None, False),
    "pielou_evenness": QIIME2DiversityMetric("pielou-evenness", "--o-vector", None, False),
    "pielou-evenness": QIIME2DiversityMetric("pielou-evenness", "--o-vector", None, False),
    "faith_pd": QIIME2DiversityMetric("faith-pd", "--o-vector", "--p-threads", True),
    "faith-pd": QIIME2DiversityMetric("faith-pd", "--o-vector", "--p-threads", True),
    "bray_curtis": QIIME2DiversityMetric("bray-curtis", "--o-distance-matrix", "--p-n-jobs", False),
    "bray-curtis": QIIME2DiversityMetric("bray-curtis", "--o-distance-matrix", "--p-n-jobs", False),
    "jaccard": QIIME2DiversityMetric("jaccard", "--o-distance-matrix", "--p-n-jobs", False),
    "unweighted_unifrac": QIIME2DiversityMetric(
        "unweighted-unifrac", "--o-distance-matrix", "--p-threads", True
    ),
    "unweighted-unifrac": QIIME2DiversityMetric(
        "unweighted-unifrac", "--o-distance-matrix", "--p-threads", True
    ),
    "weighted_unifrac": QIIME2DiversityMetric(
        "weighted-unifrac", "--o-distance-matrix", "--p-threads", True
    ),
    "weighted-unifrac": QIIME2DiversityMetric(
        "weighted-unifrac", "--o-distance-matrix", "--p-threads", True
    ),
}

SUPPORTED_METHODS = ("qiime2",)


def diversity_calc(
    *,
    backend: str,
    metric: str,
    table: Path,
    output: Path,
    phylogeny: Path | None = None,
    threads: str = "1",
    force: bool = False,
) -> None:
    backend = backend.lower()
    if backend != "qiime2":
        raise MicrobiomeSuiteError(
            f"Unsupported diversity calculation backend '{backend}'. "
            f"Choose one of: {', '.join(SUPPORTED_METHODS)}"
        )
    diversity_calc_qiime2(
        metric=metric,
        table=table,
        phylogeny=phylogeny,
        output=output,
        threads=threads,
        force=force,
    )


def diversity_calc_qiime2(
    *,
    metric: str,
    table: Path,
    output: Path,
    phylogeny: Path | None,
    threads: str,
    force: bool,
) -> None:
    metric_spec = QIIME2_DIVERSITY_METRICS.get(metric.lower())
    if metric_spec is None:
        raise MicrobiomeSuiteError(
            f"Unsupported QIIME 2 diversity metric '{metric}'. "
            f"Choose one of: {', '.join(sorted(QIIME2_DIVERSITY_METRICS))}"
        )
    if metric_spec.requires_phylogeny and phylogeny is None:
        raise MicrobiomeSuiteError(f"--phylogeny is required for metric '{metric}'.")

    qiime = shutil.which("qiime")
    if qiime is None:
        raise MicrobiomeSuiteError(
            "QIIME 2 diversity calculation requires the external 'qiime' command. "
            "Activate a QIIME 2 environment with the diversity-lib plugin and rerun this command."
        )

    ensure_input(table)
    if phylogeny is not None:
        ensure_input(phylogeny)
    prepare_output(output, force=force)

    command = [
        qiime,
        "diversity-lib",
        metric_spec.action,
        "--i-table",
        str(table),
        metric_spec.output_option,
        str(output),
    ]
    if phylogeny is not None:
        command.extend(["--i-phylogeny", str(phylogeny)])
    if metric_spec.thread_option is not None:
        command.extend([metric_spec.thread_option, threads])

    result = subprocess.run(command, check=False, text=True, capture_output=True)
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip() or "QIIME 2 diversity failed."
        raise MicrobiomeSuiteError(message)
