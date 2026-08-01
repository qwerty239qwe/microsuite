"""Shared backend-dispatch helpers for method modules.

``require_backend`` centralizes the "normalize and validate the requested
backend" step that every method function repeats, so the rejection message is
identical everywhere and the per-method dispatch no longer needs a trailing
``else: raise`` branch.
"""

from __future__ import annotations

from collections.abc import Iterable

from microsuite._errors import MicrobiomeSuiteError


def require_backend(backend: str, supported: Iterable[str], label: str) -> str:
    """Return the lower-cased backend, or raise if it is not supported."""
    options = tuple(supported)
    normalized = backend.lower()
    if normalized not in options:
        choices = ", ".join(options)
        raise MicrobiomeSuiteError(
            f"Unsupported {label} backend '{normalized}'. Choose one of: {choices}"
        )
    return normalized


def reject_options(backend: str, options: dict[str, object | None]) -> None:
    """Raise if any option that the chosen backend does not support was supplied.

    Silently ignoring an unsupported option is worse than failing: the command
    succeeds and returns a result computed without it.
    """
    rejected = [
        option
        for option, value in options.items()
        if value is not None and value is not False and value != []
    ]
    if rejected:
        raise MicrobiomeSuiteError(f"{', '.join(rejected)} not supported by --backend {backend}.")
