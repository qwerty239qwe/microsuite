from __future__ import annotations

import os
from pathlib import Path

VALID_BUILD_TARGETS: tuple[str, ...] = ("vsearch", "blast", "qiime2")


def refdb_cache_dir() -> Path:
    env = os.environ.get("MICROSUITE_REFDB_DIR")
    root = Path(env) if env else Path.home() / ".cache" / "microsuite" / "refdb"
    root.mkdir(parents=True, exist_ok=True)
    return root
