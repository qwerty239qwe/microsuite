from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_multisample_doc_exists_with_sections() -> None:
    doc = ROOT / "docs" / "multisample.md"
    assert doc.exists(), doc
    text = doc.read_text(encoding="utf-8")
    # manifest section
    assert "sample_id" in text
    # concurrency guidance
    assert ("jobs × threads" in text) or ("jobs x threads" in text)
    # decision guide references all three paths
    assert "run_fastp_multiqc.sh" in text
    assert "amplicon_qiime2" in text


def test_readme_links_to_multisample_doc() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "multisample.md" in readme
