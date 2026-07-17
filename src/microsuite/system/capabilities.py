from __future__ import annotations

from collections.abc import Iterable
from typing import Any

from microsuite.system.version import package_version

# Capability IDs are stable integration contracts. Increment ``api`` when a
# consumer must adapt to a behavior or output-schema change.
_CAPABILITIES: dict[str, dict[str, Any]] = {
    "diversity.beta_significance.permdisp.native": {"available": True, "api": 1},
    "metadata.resolved_config.v1": {"available": True, "api": 1},
    "metadata.stage_result.v1": {"available": True, "api": 1},
    "methods.dada2.454": {"available": True, "api": 1},
    "methods.dada2.parameter_sweep": {"available": True, "api": 1},
    "system.capabilities": {"available": True, "api": 1},
    "system.doctor": {"available": True, "api": 1},
    "system.version": {"available": True, "api": 1},
    "table.normalize.clr": {"available": True, "api": 1},
    "viz.metadata_grouping": {"available": True, "api": 1},
}


def capability_payload() -> dict[str, Any]:
    return {
        "schema_version": "microsuite-capabilities.v1",
        "producer": {"name": "microsuite", "version": package_version()},
        "capabilities": {key: dict(value) for key, value in sorted(_CAPABILITIES.items())},
    }


def require_capabilities(required: Iterable[str]) -> tuple[list[str], list[str]]:
    available: list[str] = []
    missing: list[str] = []
    for capability in required:
        entry = _CAPABILITIES.get(capability)
        if entry is not None and entry.get("available") is True:
            available.append(capability)
        else:
            missing.append(capability)
    return available, missing
