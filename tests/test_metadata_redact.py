from __future__ import annotations

from enum import Enum
from pathlib import Path

from microsuite.metadata.redact import (
    MASK,
    json_safe,
    redact_command,
    redact_params,
    redact_text,
)


class _Mode(Enum):
    FAST = "fast"


def test_redact_params_masks_sensitive_and_returns_secrets() -> None:
    masked, secrets = redact_params(
        {"trunc_len": 240, "auth_token": "supersecret", "nested": {"api_key": "abc123"}}
    )
    assert masked["trunc_len"] == 240
    assert masked["auth_token"] == MASK
    assert masked["nested"]["api_key"] == MASK
    assert secrets == {"supersecret", "abc123"}


def test_redact_params_json_safe_coercion() -> None:
    masked, _ = redact_params(
        {"path": Path("/a/b"), "mode": _Mode.FAST, "items": (1, 2), "s": {3, 3}}
    )
    assert masked["path"] == "/a/b"
    assert masked["mode"] == "fast"
    assert masked["items"] == [1, 2]
    assert masked["s"] == [3]


def test_redact_params_ignores_empty_secret() -> None:
    _, secrets = redact_params({"token": "   "})
    assert secrets == set()


def test_redact_command_flag_space_and_equals_and_bare() -> None:
    masked, discovered = redact_command(
        ["tool", "--token", "abc", "--api-key=xyz", "leak-abc"], {"leak-abc"}
    )
    assert masked == ["tool", "--token", MASK, "--api-key=" + MASK, MASK]
    assert discovered == {"abc", "xyz"}


def test_redact_text_longest_first_and_short_boundary() -> None:
    # command-only secret "abc" (found via --token) scrubbed from a message
    text = "Authentication rejected for abc; token abcdef also bad"
    out = redact_text(text, {"abc", "abcdef"})
    assert "abcdef" not in out
    assert "for " + MASK in out


def test_redact_text_short_secret_boundary_aware() -> None:
    # a <4-char secret only masked at word boundaries, not inside another word
    out = redact_text("pw is xy but xylophone stays", {"xy"})
    assert "xylophone" in out
    assert "is " + MASK in out


def test_json_safe_passthrough_bool() -> None:
    assert json_safe(True) is True
    assert json_safe(None) is None
