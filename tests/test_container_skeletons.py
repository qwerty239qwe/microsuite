from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTAINERS = ROOT / "containers"


def test_container_dockerfiles_exist_with_expected_tools() -> None:
    expected = {
        "microsuite": ["microsuite", "uv"],
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

    assert "| Task | Backend | Status | Image / environment | Purpose |" in text
    assert "[microsuite Python](containers/microsuite/Dockerfile)" in text
    assert "[QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile)" in text
    assert "[R diffab](containers/r-diffab/Dockerfile)" in text
    assert "[Kraken2](containers/kraken2/Dockerfile)" in text
    assert "| `qc` | `fastqc` |" in text
    assert "| `trim` | `fastp` |" in text
    assert "| `diff_abundance` | `ancombc` |" in text
