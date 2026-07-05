from __future__ import annotations

import gzip
import shutil
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


def _gunzip(path: Path) -> Path:
    if path.suffix != ".gz":
        return path
    out = path.with_suffix("")
    with gzip.open(path, "rb") as src, out.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    return out


# biodbs's SILVA_Fetcher joins relative paths against a stale, broken base URL
# (probe TL;DR #4), but an absolute URL bypasses that (broken) urljoin and
# downloads fine. The SSURef NR99 FASTA embeds taxonomy in each header as
# `>{id} {lineage}` (verified live 2026-07-06 against the real file, e.g.
# `>AY846379.1.1791 Eukaryota;Archaeplastida;...`), so one file yields both
# sequences and taxonomy.
_SILVA_HOST = "https://ftp.arb-silva.de/current/Exports"


def _silva_adapter(bd, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
    version = spec.version or "138.2"
    url = f"{_SILVA_HOST}/SILVA_{version}_SSURef_NR99_tax_silva.fasta.gz"
    raw_fa = _gunzip(Path(bd.silva_download_file(url, str(out_dir))))
    seqs = out_dir / "silva_seqs.fasta"
    tax = out_dir / "silva.tax.tsv"
    with raw_fa.open() as src, seqs.open("w") as sfh, tax.open("w") as tfh:
        for line in src:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                sid, _, lineage = header.partition(" ")
                sfh.write(f">{sid}\n")
                tfh.write(f"{sid}\t{lineage}\n")
            else:
                sfh.write(line)
    return RawRefDb(sequences=seqs, taxonomy=tax)


# GTDB's SSU (16S) reps live under the release directory, e.g.
# `latest/genomic_files_reps/bac120_ssu_reps.fna.gz` — the release segment is
# required: gtdb_download_file() joins a relative path against GTDB_Fetcher's
# base_url ("https://data.gtdb.ecogenomic.org/releases/") with a plain
# urljoin, so a path lacking the leading "{release}/" 404s (verified live
# 2026-07-06: dropping "latest/" 404s, the full path returns 200). Taxonomy is
# a headerless `{genome_accession}\t{lineage}` TSV. Verified live 2026-07-06
# that the real SSU fasta header's first token IS the bare genome accession
# (e.g. `>RS_GCF_031457235.1 d__Bacteria;p__...`) with no `~`/contig suffix —
# unlike the brief's assumed `{genome_acc}~{contig}` shape. We still split on
# "~" defensively: str.split("~", 1)[0] is a no-op when "~" is absent, so this
# same code path handles both the live no-tilde format and any historical/
# other-release file that does encode a contig suffix after "~".
def _gtdb_adapter(bd, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
    domain = "bac120"
    release = spec.version or "latest"
    ssu = _gunzip(
        Path(
            bd.gtdb_download_file(
                f"{release}/genomic_files_reps/{domain}_ssu_reps.fna.gz", str(out_dir)
            )
        )
    )
    tax_gz = Path(
        bd.gtdb_download_taxonomy(
            domain=domain, dest=str(out_dir), release=release, compressed=True
        )
    )
    tax_file = _gunzip(tax_gz)
    acc_to_lineage: dict[str, str] = {}
    for row in tax_file.read_text().splitlines():
        if not row.strip():
            continue
        acc, _, lineage = row.partition("\t")
        acc_to_lineage[acc] = lineage
    seqs = out_dir / "gtdb_seqs.fasta"
    tax = out_dir / "gtdb.tax.tsv"
    with ssu.open() as src, seqs.open("w") as sfh, tax.open("w") as tfh:
        for line in src:
            if line.startswith(">"):
                rec_id = line[1:].split()[0]  # e.g. RS_GCF_031457235.1 or RS_GCF_1~ctg1
                genome_acc = rec_id.split("~", 1)[0]
                lineage = acc_to_lineage.get(genome_acc, "")
                sfh.write(f">{rec_id}\n")
                tfh.write(f"{rec_id}\t{lineage}\n")
            else:
                sfh.write(line)
    return RawRefDb(sequences=seqs, taxonomy=tax)


_DB_ADAPTERS: dict[str, Callable[[object, RefDbSpec, Path], RawRefDb]] = {
    "homd": _homd_adapter,
    "silva": _silva_adapter,
    "gtdb": _gtdb_adapter,
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
