from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import register_provider
from microsuite.refdb.providers._base import RefDbProvider
from microsuite.refdb.spec import RawRefDb, RefDbSpec


def _load_biodbs_fetch() -> Callable[[str, str, str], tuple[str, str]]:
    from biodbs.amplicon import fetch_reference  # type: ignore[import-not-found]

    return fetch_reference


class BiodbsProvider(RefDbProvider):
    name = "biodbs"

    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            fetch = _load_biodbs_fetch()
        except ImportError as exc:
            raise MicrobiomeSuiteError(
                "The default 'biodbs' provider requires the biodbs package with "
                "amplicon-reference support. Install/upgrade biodbs, or pass a raw "
                "--classifier path, or use --provider rescript."
            ) from exc
        seqs, tax = fetch(spec.name, spec.version, str(out_dir))
        return RawRefDb(sequences=Path(seqs), taxonomy=Path(tax))


register_provider(BiodbsProvider())
