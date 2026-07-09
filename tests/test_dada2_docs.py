from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_dada2_doc_has_container_section() -> None:
    text = (ROOT / "docs" / "dada2.md").read_text(encoding="utf-8")
    assert "Running DADA2 in a container" in text
    assert "--runtime docker" in text
    assert "docker pull" in text
    assert "r-dada2" in text
    assert "MICROSUITE_R_DADA2_IMAGE" in text


def test_installation_links_to_dada2_doc() -> None:
    text = (ROOT / "docs" / "installation.md").read_text(encoding="utf-8")
    assert "dada2.md" in text
