from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
NEXTFLOW = ROOT / "workflows" / "nextflow"

pytestmark = pytest.mark.skipif(
    os.environ.get("MICROSUITE_RUN_EXTERNAL_INTEGRATION") != "1",
    reason="set MICROSUITE_RUN_EXTERNAL_INTEGRATION=1 to run external-tool integration tests",
)


@pytest.mark.parametrize("layout", ["single-end", "paired-end"])
def test_nextflow_stub_runs_amplicon_qiime2_with_trim(layout: str, tmp_path: Path) -> None:
    if shutil.which("nextflow") is None:
        pytest.skip("nextflow is not installed on PATH")

    r1 = tmp_path / "s1_R1.fastq.gz"
    r1.write_bytes(gzip.compress(b"@r\nACGT\n+\nIIII\n"))
    if layout == "paired-end":
        r2 = tmp_path / "s1_R2.fastq.gz"
        r2.write_bytes(gzip.compress(b"@r\nACGT\n+\nIIII\n"))
        row = f"s1\t{r1}\t{r2}\n"
    else:
        # single-end: empty read2 -> FASTP emits ONE trimmed file -> FASTP.out.trimmed
        # is a scalar Path (not a list). This is the exact case that regressed.
        row = f"s1\t{r1}\t\n"

    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(f"sample_id\tread1\tread2\n{row}", encoding="utf-8")
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("sample-id\tgroup\ns1\ta\n", encoding="utf-8")
    classifier = tmp_path / "classifier.qza"
    classifier.write_text("stub", encoding="utf-8")
    outdir = tmp_path / "results"

    cmd = [
        "nextflow",
        "run",
        str(NEXTFLOW / "main.nf"),
        "-stub-run",
        "-profile",
        "local",
        "--workflow",
        "amplicon_qiime2",
        "--manifest",
        str(manifest),
        "--metadata",
        str(metadata),
        "--classifier",
        str(classifier),
        "--outdir",
        str(outdir),
        "--trim",
        "true",
    ]
    result = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True, timeout=900)
    assert result.returncode == 0, result.stderr
    assert (outdir / "trim" / "fastp").exists(), "fastp outputs were not published"
