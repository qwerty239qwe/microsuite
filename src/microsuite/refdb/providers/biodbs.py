from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import register_provider
from microsuite.refdb.providers._base import RefDbProvider
from microsuite.refdb.spec import RawRefDb, RefDbSpec


def _load_biodbs():
    try:
        import biodbs  # type: ignore[import-not-found]
    except ImportError as exc:
        raise MicrobiomeSuiteError(
            "The default 'biodbs' provider requires biodbs>=0.4.0. Install the "
            "refdb extra (e.g. `uv sync --extra refdb`), pass a raw --classifier "
            "path, or use --provider rescript."
        ) from exc
    return biodbs


# HOMD serves 16S RefSeq under a `current/` symlink dir; the wrapper
# homd_download_16s_refseq() is broken (probe TL;DR #2/#5), so download the
# fasta and its sibling QIIME-style taxonomy via the generic homd_download_file.
_HOMD_BASE = "ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/current"
_HOMD_STEM = "HOMD_16S_rRNA_RefSeq_V16.03"


def _homd_adapter(bd, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
    seqs = Path(bd.homd_download_file(f"{_HOMD_BASE}/{_HOMD_STEM}.fasta", str(out_dir)))
    # The .qiime.taxonomy file is already an id-first `seqID<TAB>lineage` TSV.
    tax = Path(bd.homd_download_file(f"{_HOMD_BASE}/{_HOMD_STEM}.qiime.taxonomy", str(out_dir)))
    return RawRefDb(sequences=seqs, taxonomy=tax)


_DB_ADAPTERS: dict[str, Callable[[object, RefDbSpec, Path], RawRefDb]] = {
    "homd": _homd_adapter,
}


class BiodbsProvider(RefDbProvider):
    name = "biodbs"

    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
        out_dir.mkdir(parents=True, exist_ok=True)
        adapter = _DB_ADAPTERS.get(spec.name.lower())
        if adapter is None:
            supported = ", ".join(sorted(_DB_ADAPTERS))
            raise MicrobiomeSuiteError(
                f"biodbs provider has no adapter for DB '{spec.name}'. Supported: {supported}."
            )
        bd = _load_biodbs()
        return adapter(bd, spec, out_dir)


register_provider(BiodbsProvider())
