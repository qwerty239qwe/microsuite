from __future__ import annotations

import gzip
import shutil
import tarfile
import zipfile
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


def _extract_zip_member(zip_path: Path, suffix: str, out_path: Path) -> Path:
    """Extract the single member of `zip_path` whose name ends with `suffix`
    (e.g. a QIIME2 .qza artifact's "data/dna-sequences.fasta") to `out_path`."""
    with zipfile.ZipFile(zip_path) as zf:
        candidates = [n for n in zf.namelist() if n.endswith(suffix)]
        if len(candidates) != 1:
            raise MicrobiomeSuiteError(
                f"Expected exactly one member ending with '{suffix}' in {zip_path}, "
                f"found {len(candidates)}: {candidates}."
            )
        with zf.open(candidates[0]) as src, out_path.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    return out_path


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


# GreenGenes ships QIIME2 .qza artifacts (which are plain ZIP files) under a
# release directory, e.g. "2022.7-rc1/2022.7.backbone.v4.fna.qza" — like GTDB,
# greengenes_download_file() needs the release segment in the path (verified
# live 2026-07-06: "2022.7.backbone.tax.qza" with no release prefix 404s;
# "2022.7-rc1/2022.7.backbone.tax.qza" returns 200). Verified live that the
# sequence .qza's data/dna-sequences.fasta member already has id-only headers
# (e.g. ">MJ006-1-barcode39-umi49105bins-ubs-7"), and the taxonomy .qza's
# data/taxonomy.tsv has a QIIME2 header row ("Feature ID\tTaxon\n") followed by
# id-first rows whose ids match the fasta headers exactly — so we strip the
# header and keep columns 0 (id) and 1 (lineage).
_GREENGENES_RELEASE = "2022.7-rc1"
_GREENGENES_SEQ_FILE = "2022.7.backbone.v4.fna.qza"
_GREENGENES_TAX_FILE = "2022.7.backbone.tax.qza"


def _greengenes_adapter(bd, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
    release = spec.version or _GREENGENES_RELEASE
    seq_qza = Path(bd.greengenes_download_file(f"{release}/{_GREENGENES_SEQ_FILE}", str(out_dir)))
    tax_qza = Path(bd.greengenes_download_file(f"{release}/{_GREENGENES_TAX_FILE}", str(out_dir)))
    seqs = _extract_zip_member(
        seq_qza, "data/dna-sequences.fasta", out_dir / "greengenes_seqs.fasta"
    )
    raw_tax = _extract_zip_member(tax_qza, "data/taxonomy.tsv", out_dir / "greengenes_raw_tax.tsv")
    tax = out_dir / "greengenes.tax.tsv"
    with raw_tax.open() as src, tax.open("w") as dst:
        rows = iter(src)
        next(rows, None)  # drop the QIIME2 header row ("Feature ID\tTaxon[\tConfidence]")
        for row in rows:
            if not row.strip():
                continue
            fields = row.rstrip("\n").split("\t")
            dst.write(f"{fields[0]}\t{fields[1]}\n")
    return RawRefDb(sequences=seqs, taxonomy=tax)


# UNITE ships ONE .tgz containing several RESCRIPt-style fasta+taxonomy pairs
# clustered at different thresholds (e.g. "dynamic", "97", "99"), each
# duplicated again under a developer/ subdirectory (verified live 2026-07-06
# against the 2020-02-20 fungi release: 6 top-level *.fasta candidates and 6
# *.txt candidates, not the single pair the initial brief assumed). The
# "dynamic" clustering is UNITE's recommended general-purpose release, so we
# select the member containing "dynamic" in its name, outside any developer/
# path, and require the selection to be unambiguous. Fasta headers are
# id-only (">SH1546528.08FU_JF832665_refs") and the matching taxonomy .txt is
# already an id-first "accession\tlineage" TSV with no header row.
def _select_unite_member(names: list[str], *, suffix: str, keyword: str = "dynamic") -> str:
    candidates = [n for n in names if n.endswith(suffix) and "developer" not in Path(n).parts]
    preferred = [n for n in candidates if keyword in Path(n).name]
    pool = preferred or candidates
    if len(pool) != 1:
        raise MicrobiomeSuiteError(
            f"Expected exactly one non-developer '{keyword}' member ending with "
            f"'{suffix}' in the UNITE archive, found {len(pool)}: {pool or candidates}."
        )
    return pool[0]


def _unite_adapter(bd, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
    version = spec.version or "2020-02-20"
    archive = Path(bd.unite_download(version, str(out_dir), taxon_group="fungi", singletons=False))
    seqs = out_dir / "unite_seqs.fasta"
    tax = out_dir / "unite.tax.tsv"
    with tarfile.open(archive, "r:*") as tf:
        names = tf.getnames()
        fasta_name = _select_unite_member(names, suffix=".fasta")
        tax_name = _select_unite_member(names, suffix=".txt")
        fasta_src = tf.extractfile(fasta_name)
        if fasta_src is None:
            raise MicrobiomeSuiteError(f"UNITE archive member '{fasta_name}' is not a file.")
        with seqs.open("wb") as dst:
            shutil.copyfileobj(fasta_src, dst)
        tax_src = tf.extractfile(tax_name)
        if tax_src is None:
            raise MicrobiomeSuiteError(f"UNITE archive member '{tax_name}' is not a file.")
        with tax.open("wb") as dst:
            shutil.copyfileobj(tax_src, dst)
    return RawRefDb(sequences=seqs, taxonomy=tax)


# PR2 ships sequences and taxonomy as two separate gzipped GitHub-release
# assets in "mothur" format: "pr2_version_{v}_SSU_mothur.fasta.gz" (id-only
# headers, e.g. ">AB353770.1.1740_U") and "pr2_version_{v}_SSU_mothur.tax.gz"
# (id-first "accession\tlineage;" TSV, headerless, verified live 2026-07-06).
# The lineage has a trailing ";" that we strip so it matches the SILVA/GTDB
# adapters' lineage shape.
_PR2_FASTA_TMPL = "pr2_version_{v}_SSU_mothur.fasta.gz"
_PR2_TAX_TMPL = "pr2_version_{v}_SSU_mothur.tax.gz"


def _pr2_adapter(bd, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
    version = spec.version or "5.1.1"
    seq_gz = Path(bd.pr2_download_asset(_PR2_FASTA_TMPL.format(v=version), str(out_dir)))
    tax_gz = Path(bd.pr2_download_asset(_PR2_TAX_TMPL.format(v=version), str(out_dir)))
    seqs = _gunzip(seq_gz)
    raw_tax = _gunzip(tax_gz)
    tax = out_dir / "pr2.tax.tsv"
    with raw_tax.open() as src, tax.open("w") as dst:
        for row in src:
            if not row.strip():
                continue
            acc, _, lineage = row.rstrip("\n").partition("\t")
            dst.write(f"{acc}\t{lineage.rstrip(';')}\n")
    return RawRefDb(sequences=seqs, taxonomy=tax)


_DB_ADAPTERS: dict[str, Callable[[object, RefDbSpec, Path], RawRefDb]] = {
    "homd": _homd_adapter,
    "silva": _silva_adapter,
    "gtdb": _gtdb_adapter,
    "greengenes": _greengenes_adapter,
    "unite": _unite_adapter,
    "pr2": _pr2_adapter,
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
