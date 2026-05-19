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
