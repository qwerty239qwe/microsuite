from __future__ import annotations

from microsuite.runtime.container import resolve_batch_image


def test_default_image_is_per_backend() -> None:
    assert resolve_batch_image("mmuphin", None) == (
        "ghcr.io/qwerty239qwe/microsuite/r-batch-mmuphin:latest"
    )


def test_explicit_override_wins_over_env(monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_R_BATCH_MMUPHIN_IMAGE", "from-env:1")
    assert resolve_batch_image("mmuphin", "explicit:2") == "explicit:2"


def test_env_override_is_uppercased_and_underscored(monkeypatch) -> None:
    # The backend is 'combat-seq'; the env var cannot contain a hyphen.
    monkeypatch.setenv("MICROSUITE_R_BATCH_COMBAT_SEQ_IMAGE", "from-env:1")
    assert resolve_batch_image("combat-seq", None) == "from-env:1"
