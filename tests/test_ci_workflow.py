from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_github_actions_ci_runs_project_quality_gate() -> None:
    workflow = ROOT / ".github" / "workflows" / "ci.yml"
    assert workflow.exists()

    text = workflow.read_text(encoding="utf-8")
    assert "astral-sh/setup-uv" in text
    assert 'python-version: ["3.11", "3.12"]' in text
    assert "UV_PYTHON: ${{ matrix.python-version }}" in text
    assert "enable-cache: true" in text
    assert "cache-dependency-glob: uv.lock" in text
    assert "uv sync --all-extras --locked" in text
    assert "--extra dev --all-extras" not in text
    assert "uv run pytest" in text
    assert "uv run ruff check ." in text
    assert "uv run ruff format --check ." in text
    assert "uv run ty check" in text
    assert "uv build" in text
    assert "nextflow-smoke:" in text
    assert "nf-core/setup-nextflow@v2" in text
    assert "nextflow-io/setup-nextflow" not in text
    assert 'version: "25.10.0"' in text
    assert "nextflow run workflows/nextflow/main.nf" in text
    assert "--outdir results/nextflow-smoke" in text


def test_source_data_package_is_not_ignored() -> None:
    assert (ROOT / "src" / "microsuite" / "data" / "__init__.py").exists()
    assert (ROOT / "src" / "microsuite" / "data" / "moving_pictures.py").exists()
    assert (
        ROOT / "src" / "microsuite" / "data" / "fixtures" / "moving_pictures_small" / "table.tsv"
    ).exists()


def test_demo_data_attribution_is_documented() -> None:
    attribution = (ROOT / "docs" / "data-attribution.md").read_text(encoding="utf-8")
    fixture_readme = (
        ROOT / "src" / "microsuite" / "data" / "fixtures" / "moving_pictures_small" / "README.md"
    ).read_text(encoding="utf-8")

    for text in (attribution, fixture_readme):
        assert "not owned by this project" in text
        assert "10.1186/gb-2011-12-5-r50" in text
        assert "QIIME 2 Moving Pictures tutorial" in text
