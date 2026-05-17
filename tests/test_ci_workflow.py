from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_github_actions_ci_runs_project_quality_gate() -> None:
    workflow = ROOT / ".github" / "workflows" / "ci.yml"
    assert workflow.exists()

    text = workflow.read_text(encoding="utf-8")
    assert "astral-sh/setup-uv" in text
    assert 'python-version: ["3.11", "3.12"]' in text
    assert "enable-cache: true" in text
    assert "cache-dependency-glob: uv.lock" in text
    assert "uv sync --extra dev --all-extras" in text
    assert "uv run pytest" in text
    assert "uv run ruff check ." in text
    assert "uv run ruff format --check ." in text
    assert "uv run ty check" in text
    assert "uv build" in text
