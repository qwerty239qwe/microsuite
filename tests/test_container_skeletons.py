from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CONTAINERS = ROOT / "containers"


def test_container_dockerfiles_exist_with_expected_tools() -> None:
    expected = {
        "microsuite": ["microsuite", "uv"],
        "fastqc": ["fastqc", "openjdk"],
        "qiime2-amplicon": ["qiime", "QIIME 2"],
        "r-diffab": ["Rscript", "ANCOMBC"],
        "kraken2": ["kraken2", "Bracken support is planned"],
    }

    for name, tokens in expected.items():
        dockerfile = CONTAINERS / name / "Dockerfile"
        assert dockerfile.exists(), name
        text = dockerfile.read_text(encoding="utf-8")
        assert "org.opencontainers.image.title" in text
        assert "org.opencontainers.image.description" in text
        assert "# Expected commands:" in text
        for token in tokens:
            assert token in text


def test_dockerignore_excludes_local_artifacts() -> None:
    text = (ROOT / ".dockerignore").read_text(encoding="utf-8")

    for pattern in [".venv", ".git", "dist", "/runs", "/data", "/results", ".tokensave"]:
        assert pattern in text
    assert "src/microsuite/data" not in text


def test_readme_method_surface_links_backends_to_environments() -> None:
    text = (ROOT / "README.md").read_text(encoding="utf-8")

    assert "### Quality Reports" in text
    assert "### Quality Filtering" in text
    assert "### Differential Abundance" in text
    assert "| Subtopic |" not in text
    assert "Python function name" not in text
    assert "microsuite.api." not in text
    assert "microsuite.methods." not in text
    assert "microsuite.viz." not in text
    assert (
        "| Backend | Version | Status | CLI command | Python invocation | "
        "Image / environment | Operational tradeoff | Purpose |"
    ) in text
    assert "[microsuite Python](containers/microsuite/Dockerfile)" in text
    assert "[FastQC](containers/fastqc/Dockerfile)" in text
    assert "[QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile)" in text
    assert "[R diffab](containers/r-diffab/Dockerfile)" in text
    assert "[Kraken2](containers/kraken2/Dockerfile)" in text
    assert (
        "| `fastqc` | FastQC 0.12.1 | ready | "
        "`microsuite qc --backend fastqc` | "
        '`qc(backend="fastqc", inputs=[...], output_dir=...)` |'
    ) in text
    assert (
        "| `qiime2-demux-summarize` | QIIME 2 2024.10 | partial | "
        "`microsuite qc --backend qiime2-demux` | "
        '`qc(backend="qiime2-demux", demux=..., output=...)` |'
    ) in text
    assert (
        "| `qiime2-exclude-seqs` | QIIME 2 2026.4 API | partial | "
        "`microsuite qc_filter --backend qiime2-exclude-seqs` | "
        '`qc_filter(backend="qiime2-exclude-seqs", query_sequences=..., '
        "reference_sequences=...)` |"
    ) in text
    assert (
        "| `qiime2-bowtie2-build` | QIIME 2 2026.4 API | partial | "
        "`microsuite qc_filter --backend qiime2-bowtie2-build` | "
        '`qc_filter(backend="qiime2-bowtie2-build", sequences=..., output=...)` |'
    ) in text
    assert (
        "| `qiime2-decontam` | QIIME 2 2026.4 API | partial | "
        "`microsuite decontam --backend qiime2-decontam` | "
        '`decontam(backend="qiime2-decontam", table=..., metadata=..., output=...)` |'
    ) in text
    assert (
        "| `fastp` | User env | partial | "
        "`microsuite trim --backend fastp` | "
        '`trim(backend="fastp", read1=..., output1=...)` |'
    ) in text
    assert (
        "| `cutadapt` | Cutadapt >=4.x user env | partial | "
        "`microsuite trim --backend cutadapt` | "
        '`trim(backend="cutadapt", read1=..., output1=..., adapter=...)` |'
    ) in text
    assert (
        "| `trimmomatic` | Trimmomatic >=0.39 user env | partial | "
        "`microsuite trim --backend trimmomatic` | "
        '`trim(backend="trimmomatic", read1=..., output1=..., trimmomatic_steps=[...])` |'
    ) in text
    assert (
        "| `trim-galore` | Trim Galore >=0.6 user env | partial | "
        "`microsuite trim --backend trim-galore` | "
        '`trim(backend="trim-galore", read1=..., output1=..., adapter=...)` |'
    ) in text
    assert (
        "| `dada2-r` | DADA2 R user env | partial | "
        "`microsuite denoise --backend dada2-r` | "
        '`denoise(backend="dada2-r", demux=reads_dir, output_table=...)` |'
    ) in text
    assert (
        "| `ancombc` | R 4.4.0 image; ANCOMBC via "
        "Bioconductor | partial | `microsuite diff_abundance --backend ancombc` | "
        '`diff_abundance(backend="ancombc", table=..., group=..., output=...)` |'
    ) in text


@pytest.mark.skipif(
    os.environ.get("MICROSUITE_RUN_DOCKER_TESTS") != "1",
    reason="set MICROSUITE_RUN_DOCKER_TESTS=1 to run Docker-backed integration tests",
)
def test_fastqc_container_runs_tiny_fastq_end_to_end(tmp_path: Path) -> None:
    if shutil.which("docker") is None:
        pytest.skip("docker is not installed")

    image = "microsuite/fastqc:ci"
    image_check = subprocess.run(
        ["docker", "image", "inspect", image],
        check=False,
        text=True,
        capture_output=True,
        timeout=30,
    )
    if image_check.returncode != 0:
        pytest.skip(f"Docker image is not available locally: {image}")

    fixture = ROOT / "tests" / "fixtures" / "fastq" / "tiny.fastq"
    output_dir = tmp_path / "fastqc"
    output_dir.mkdir()

    command = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{fixture.parent.resolve()}:/input:ro",
        "-v",
        f"{output_dir.resolve()}:/output",
        image,
        "--outdir",
        "/output",
        "/input/tiny.fastq",
    ]
    try:
        result = subprocess.run(
            command,
            check=False,
            text=True,
            capture_output=True,
            timeout=90,
        )
    except subprocess.TimeoutExpired as exc:
        raise AssertionError(f"FastQC container timed out: {' '.join(command)}") from exc

    assert result.returncode == 0, result.stderr or result.stdout
    assert (output_dir / "tiny_fastqc.html").exists()
    assert (output_dir / "tiny_fastqc.zip").exists()
