from __future__ import annotations

from pathlib import Path

from microsuite.batch.backends import BATCH_BACKENDS
from microsuite.runtime.container import resolve_batch_image

_REPO_ROOT = Path(__file__).resolve().parents[1]


def test_default_image_is_per_backend() -> None:
    assert resolve_batch_image("mmuphin", None, image=BATCH_BACKENDS["mmuphin"].image) == (
        "ghcr.io/qwerty239qwe/microsuite/r-batch-mmuphin:latest"
    )


def test_explicit_override_wins_over_env(monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_R_BATCH_MMUPHIN_IMAGE", "from-env:1")
    assert resolve_batch_image("mmuphin", "explicit:2") == "explicit:2"


def test_env_override_is_uppercased_and_underscored(monkeypatch) -> None:
    # The backend is 'combat-seq'; the env var cannot contain a hyphen.
    monkeypatch.setenv("MICROSUITE_R_BATCH_COMBAT_SEQ_IMAGE", "from-env:1")
    assert resolve_batch_image("combat-seq", None) == "from-env:1"


def test_every_backend_default_image_matches_container_dir_and_ci() -> None:
    """The capability table, the container tree, and CI must agree on names.

    This is the invariant Fix 1 exists to enforce: `resolve_batch_image`'s
    default (via each record's declared `image`) must name a real
    `containers/<basename>/` directory and a real `image:` entry in the
    Docker CI matrix -- not just the backend name with hyphens stripped.
    """
    workflow_text = (_REPO_ROOT / ".github/workflows/docker.yml").read_text(encoding="utf-8")

    for backend, record in BATCH_BACKENDS.items():
        resolved = resolve_batch_image(backend, None, image=record.image)
        basename = resolved.rsplit("/", 1)[-1].split(":", 1)[0]
        assert basename == record.image, backend
        assert (_REPO_ROOT / "containers" / basename).is_dir(), (
            f"{backend}: no containers/{basename}/ directory"
        )
        assert f"image: {basename}" in workflow_text, (
            f"{backend}: {basename} not built by .github/workflows/docker.yml"
        )
