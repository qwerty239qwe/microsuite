from __future__ import annotations

import shutil
import urllib.request
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError

FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "moving_pictures_small"
REMOTE_FILES = {
    "sample-metadata.tsv": (
        "https://data.qiime2.org/2024.10/tutorials/moving-pictures/sample_metadata.tsv"
    ),
    "table.qza": "https://docs.qiime2.org/2024.10/data/tutorials/moving-pictures/table.qza",
    "taxonomy.qza": "https://docs.qiime2.org/2024.10/data/tutorials/moving-pictures/taxonomy.qza",
}


def copy_small_fixture(output: Path, *, force: bool = False) -> None:
    if not FIXTURE_DIR.exists():
        raise MicrobiomeSuiteError(f"Bundled fixture directory is missing: {FIXTURE_DIR}")
    output.mkdir(parents=True, exist_ok=True)
    for source in FIXTURE_DIR.iterdir():
        target = output / source.name
        if target.exists() and not force:
            raise MicrobiomeSuiteError(f"Output exists, pass --force to overwrite: {target}")
        shutil.copy2(source, target)


def fetch_moving_pictures(output: Path, *, force: bool = False) -> None:
    output.mkdir(parents=True, exist_ok=True)
    for filename, url in REMOTE_FILES.items():
        target = output / filename
        if target.exists() and not force:
            raise MicrobiomeSuiteError(f"Output exists, pass --force to overwrite: {target}")
        try:
            urllib.request.urlretrieve(url, target)
        except OSError as exc:
            raise MicrobiomeSuiteError(f"Failed to download {url}: {exc}") from exc
