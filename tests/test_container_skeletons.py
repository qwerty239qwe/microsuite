from __future__ import annotations

from pathlib import Path

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
    assert (
        "| Backend | Version | Status | CLI command | Python function name | "
        "Image / environment | Operational tradeoff | Purpose |"
    ) in text
    assert "[microsuite Python](containers/microsuite/Dockerfile)" in text
    assert "[FastQC](containers/fastqc/Dockerfile)" in text
    assert "[QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile)" in text
    assert "[R diffab](containers/r-diffab/Dockerfile)" in text
    assert "[Kraken2](containers/kraken2/Dockerfile)" in text
    assert (
        "| `fastqc` | FastQC 0.12.1 | ready | "
        "`microsuite qc --backend fastqc` | `microsuite.methods.qc.qc` |"
    ) in text
    assert (
        "| `qiime2-demux-summarize` | QIIME 2 2024.10 | partial | "
        "`microsuite qc --backend qiime2-demux` | `microsuite.methods.qc.qc` |"
    ) in text
    assert "| `qiime2-quality-control-exclude-seqs` | QIIME 2 2024.10 | planned |" in text
    assert (
        "| `fastp` | User env | partial | "
        "`microsuite trim --backend fastp` | `microsuite.methods.trim.trim` |"
    ) in text
    assert (
        "| `ancombc` | R 4.4.0 image; ANCOMBC via "
        "Bioconductor | partial | `microsuite diff_abundance --backend ancombc` | "
        "`microsuite.methods.diff_abundance.diff_abundance` |"
    ) in text
