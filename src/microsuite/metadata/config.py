"""Resolved-configuration snapshot writer for the end-to-end runner.

No single microsuite command owns the whole analysis configuration, so the
orchestrator aggregates defaults + overrides and calls ``write_resolved_config``
to persist a redacted, versioned snapshot next to the run. Per-stage resolved
parameters already live in each ``stage-result.v1`` envelope's ``params``.
"""

from __future__ import annotations

import time
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from microsuite import __version__
from microsuite._errors import MicrobiomeSuiteError
from microsuite.metadata.redact import redact_params
from microsuite.metadata.schemas import RESOLVED_CONFIG_VERSION
from microsuite.metadata.stage import _atomic_write
from microsuite.metadata.validate import validate

DEFAULT_FILENAME = "resolved_config.json"


def _utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def write_resolved_config(
    run_dir: Path, config: Mapping[str, Any], *, name: str = DEFAULT_FILENAME
) -> Path:
    """Write a redacted ``resolved-config.v1`` snapshot (defaults in, secrets out)."""
    masked, _ = redact_params(config)
    payload = {
        "schema_version": RESOLVED_CONFIG_VERSION,
        "generated_at": _utc_now(),
        "producer": {"name": "microsuite", "version": __version__},
        "config": masked,
    }
    errors = validate(payload, RESOLVED_CONFIG_VERSION)
    if errors:
        raise MicrobiomeSuiteError(f"Refusing to write an invalid resolved-config: {errors}")
    run_dir.mkdir(parents=True, exist_ok=True)
    target = run_dir / name
    _atomic_write(target, payload)
    return target
