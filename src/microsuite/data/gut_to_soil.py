from __future__ import annotations

import urllib.request
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError

BASE_URL = "https://gut-to-soil-tutorial.readthedocs.io/en/2026.4/data/gut-to-soil"

REMOTE_FILES = {
    "sample-metadata.tsv": f"{BASE_URL}/sample-metadata.tsv",
    "asv-table-ms2.qza": f"{BASE_URL}/asv-table-ms2.qza",
    "taxonomy.qza": f"{BASE_URL}/taxonomy.qza",
}


def fetch_gut_to_soil(output: Path, *, force: bool = False) -> None:
    output.mkdir(parents=True, exist_ok=True)
    for filename, url in REMOTE_FILES.items():
        target = output / filename
        if target.exists() and not force:
            continue
        try:
            _download(url, target)
        except OSError as exc:
            raise MicrobiomeSuiteError(f"Failed to download {url}: {exc}") from exc


def _download(url: str, target: Path) -> None:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "microsuite-real-data-integration/0.1"},
    )
    with urllib.request.urlopen(request) as response:
        target.write_bytes(response.read())
