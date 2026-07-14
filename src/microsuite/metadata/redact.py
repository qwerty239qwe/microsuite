from __future__ import annotations

"""Shared secret redaction for the ``stage-result.v1`` envelope.

One mechanism scrubs ``params``, every subprocess ``command``, and the free-text
``error.message``. Values under sensitive keys/flags become ``"***"``; the raw
values are collected so a secret that appears only on a command line (never in
``params``) is still removed from the error message. All parameter values are
coerced to JSON-safe types so a valid payload never emits tokens Microboard's
``JSON.parse`` rejects.
"""

import re
from enum import Enum
from pathlib import Path
from typing import Any

MASK = "***"
SENSITIVE_KEY_RE = re.compile(r"(?i)(token|secret|password|passwd|api[-_]?key|credential|auth)")


def _is_secret_value(value: Any) -> bool:
    return isinstance(value, str) and value.strip() != ""


def json_safe(value: Any) -> Any:
    """Coerce a value into JSON-serialisable form (no NaN/inf handling here)."""
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, Enum):
        return json_safe(value.value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [json_safe(v) for v in value]
    return str(value)


def redact_params(mapping: Any) -> tuple[Any, set[str]]:
    """Return a JSON-safe deep copy with sensitive keys masked + the secret values."""
    secrets: set[str] = set()

    def _walk(value: Any, key_is_sensitive: bool) -> Any:
        if isinstance(value, dict):
            out: dict[str, Any] = {}
            for k, v in value.items():
                sensitive = bool(SENSITIVE_KEY_RE.search(str(k)))
                if sensitive and not isinstance(v, (dict, list, tuple, set)):
                    if _is_secret_value(v):
                        secrets.add(v)
                    out[str(k)] = MASK
                else:
                    out[str(k)] = _walk(v, sensitive)
            return out
        if isinstance(value, (list, tuple, set)):
            return [_walk(v, key_is_sensitive) for v in value]
        return json_safe(value)

    masked = _walk(mapping, False)
    return masked, secrets


def _flag_name(arg: str) -> str:
    return arg.lstrip("-").split("=", 1)[0]


def redact_command(argv: list[str], secrets: set[str]) -> tuple[list[str], set[str]]:
    """Mask sensitive flag values + known bare secrets; return masked argv + discovered."""
    masked: list[str] = []
    discovered: set[str] = set()
    skip_next = False
    for i, arg in enumerate(argv):
        if skip_next:
            skip_next = False
            continue
        if arg.startswith("-") and "=" in arg and SENSITIVE_KEY_RE.search(_flag_name(arg)):
            flag, value = arg.split("=", 1)
            if _is_secret_value(value):
                discovered.add(value)
            masked.append(f"{flag}={MASK}")
            continue
        if arg.startswith("-") and SENSITIVE_KEY_RE.search(_flag_name(arg)):
            nxt = argv[i + 1] if i + 1 < len(argv) else None
            if nxt is not None:
                if _is_secret_value(nxt):
                    discovered.add(nxt)
                masked.append(arg)
                masked.append(MASK)
                skip_next = True
                continue
            masked.append(arg)
            continue
        if arg in secrets:
            masked.append(MASK)
            continue
        masked.append(arg)
    return masked, discovered


def redact_text(text: str, secrets: set[str]) -> str:
    """Scrub every non-empty secret from free text, longest-first."""
    out = text
    for secret in sorted((s for s in secrets if s), key=len, reverse=True):
        if len(secret) < 4:
            out = re.sub(rf"(?<!\w){re.escape(secret)}(?!\w)", MASK, out)
        else:
            out = out.replace(secret, MASK)
    return out
