from __future__ import annotations

from pathlib import Path

DOCS = Path(__file__).resolve().parents[1] / "docs"


def test_three_api_docs_exist_and_name_api_boundaries() -> None:
    expected = [
        "api-cli.md",
        "api-python.md",
        "api-nextflow.md",
        "containers.md",
        "three-api-roadmap.md",
    ]
    for name in expected:
        assert (DOCS / name).exists(), name

    roadmap = (DOCS / "three-api-roadmap.md").read_text(encoding="utf-8")
    assert "Nextflow API" in roadmap
    assert "CLI API" in roadmap
    assert "Python SDK" in roadmap


def test_existing_docs_reflect_three_api_model() -> None:
    readme = (DOCS.parent / "README.md").read_text(encoding="utf-8")
    toolbox = (DOCS / "toolbox.md").read_text(encoding="utf-8")

    for text in [readme, toolbox]:
        assert "Nextflow" in text
        assert "CLI" in text
        assert "Python SDK" in text

    assert "The CLI is one public surface" in toolbox
    assert "The public surface is the `microbiome` CLI" not in toolbox

    nextflow = (DOCS / "api-nextflow.md").read_text(encoding="utf-8")
    assert "Current status" in nextflow
    assert "module files are placeholders" in nextflow

    python_api = (DOCS / "api-python.md").read_text(encoding="utf-8")
    assert "read_table" in python_api
    assert ".h5ad" in python_api
