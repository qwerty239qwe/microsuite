from __future__ import annotations

import gzip
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError


def validate_output_file(path: Path, *, allow_empty: bool = False) -> None:
    if not path.exists():
        raise MicrobiomeSuiteError(f"Expected output was not created: {path}")
    if not allow_empty and path.stat().st_size == 0:
        raise MicrobiomeSuiteError(
            f"Output is empty: {path}. This can mean an incomplete run or an "
            "unsynced cloud-storage placeholder file."
        )
    if path.name.endswith(".gz"):
        try:
            with gzip.open(path, "rb") as handle:
                handle.read(65536)
        except (gzip.BadGzipFile, EOFError, OSError) as exc:
            raise MicrobiomeSuiteError(f"Output is not a valid gzip file: {path}.") from exc


def validate_outputs(outputs: dict[str, str], *, allow_empty: bool = False) -> None:
    for path_str in outputs.values():
        validate_output_file(Path(path_str), allow_empty=allow_empty)
