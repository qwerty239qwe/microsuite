# biodbs v0.4.0 API probe (Task R1)

Captured 2026-07-05 against the real, live HOMD / SILVA / GTDB / GreenGenes /
UNITE / PR2 endpoints, using biodbs installed from the git tag (PyPI does not
yet have 0.4.0):

```
uv pip install "git+https://github.com/qwerty239qwe/biodbs.git@biodbs_v0.4.0"
```

This pinned to commit `96cab22908ef9603d3f7bcddcecfe68465756553` (tag
`biodbs_v0.4.0`), `biodbs.__version__ == "0.4.0"`.

This document is the source of truth for R2–R4's per-DB provider code and for
R5's live-integration test fixtures. Findings below come from real network
calls (a handful of small listing/HEAD/GET requests, no multi-hundred-MB
downloads).

## TL;DR — the four things that will bite R2

1. **`TableData.to_csv(path)` always raises `NotImplementedError`** for every
   DB's table class (`HOMDTableData`, `GTDBTableData` — SILVA/GreenGenes/
   UNITE/PR2 don't even have a table class). Only the base
   `BaseFetchedData.to_csv` exists and it is never overridden anywhere in
   0.4.0. **R2's `_write_taxonomy_tsv` must always go through
   `.as_dataframe(engine="pandas").to_csv(path, sep="\t", index=False)`**
   (or hand-roll from `.as_dict()`/`.rows`), never call `.to_csv()` on the
   biodbs object itself.
2. **`homd_get_hmt_lineage()` returns a mis-parsed table.** The live file has
   a junk title line before the real tab-separated header, so
   `csv.DictReader` treats the title as a 1-column header and stuffs the
   real header's fields into a list under the `None` restkey. `.as_dataframe()`
   "succeeds" (no exception) but the DataFrame is garbage (columns
   `['HOMD.org Taxon Data::Taxonomic Lineage', None]`). **Do not use this
   function for HOMD taxonomy.** Use `homd_get_gtdb_taxonomy()` instead (see
   HOMD table below) — it parses cleanly with a proper header and a
   `GENOME_ID` first column, or fetch the raw text via `homd_get_text(url)`
   and parse manually (skip line 1) if the HMT lineage's Domain..Species
   columns are specifically required.
3. **`homd_get_taxon_table()` is broken.** Its keyword-based link matching
   (`_get_table_by_keywords("taxon")`) currently resolves to an HTML page,
   not a data file — `.as_dataframe()` "succeeds" but returns a nonsense
   single column (`<style>`) full of CSS text. Do not use it.
4. **SILVA's base URL is stale.** `SILVA_Fetcher` hits
   `https://www.arb-silva.de/current-release/`, which used to be a plain
   Apache autoindex but is now a TYPO3 CMS page — it returns real HTML (200
   OK) but contains none of the `<a href>` autoindex links biodbs's parser
   expects. Result: `silva_list_current_files()` returns an **empty**
   `SILVAFileListData` (len 0), and `silva_get_version()` /
   `silva_get_readme()` / `silva_get_citation()` return the CMS page's HTML
   instead of the real text file. The real files still exist, just on a
   different host: `https://ftp.arb-silva.de/current/` (still a working
   nginx autoindex). **R2/R5 must not rely on `silva_list_current_files()`
   or `silva_get_version()`**; either hardcode the known-good
   `ftp.arb-silva.de` paths, or pass a full absolute URL to
   `silva_download_file()`/`silva_get_...()`-style calls (absolute URLs
   bypass the (broken) base_url join and download fine). Treat SILVA listing
   support in biodbs 0.4.0 as non-functional until upstream fixes it.

5. (Bonus / lower severity) **`homd_download_16s_refseq()`'s zero-arg default
   is also broken** against the live site: it auto-discovers a filename by
   listing `ftp/16S_rRNA_refseq/` directly, but the real fasta files now live
   one level deeper, under
   `ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/current/`. Auto-discovery raises
   `APIValidationError: No matching file found in '16S_rRNA_refseq'`, and
   even passing `filename="HOMD_16S_rRNA_RefSeq_V16.03.fasta"` explicitly
   fails (404) because the function's path template omits the subdirectory.
   **Workaround:** use the lower-level `homd_download_file()` (or
   `homd_list_ftp()` to discover the URL, then `homd_download_file(url, dest)`)
   with the full path
   `ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/current/HOMD_16S_rRNA_RefSeq_V16.03.fasta`.

---

## Per-DB findings

### HOMD

| Aspect | Finding |
|---|---|
| Sequence fetch fn | `homd_download_file(path_or_url, dest)` (generic) — `homd_download_16s_refseq(dest, filename="", overwrite=False)` exists but its zero-arg auto-discovery **and** its path template are broken against the live site (see bullet 5 above); call `homd_download_file` directly with a full discovered path/URL instead. |
| Taxonomy fetch fn | `homd_get_hmt_lineage()` (BROKEN parse, see above), `homd_get_taxon_table()` (BROKEN, resolves to HTML), `homd_get_gtdb_taxonomy()` (WORKS — recommended) |
| Return type (listing) | `homd_list_16s_refseq()` → `HOMDFileListData` (3 subdirectory entries, not files — must descend via `homd_list_ftp("16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/current")`) |
| Return type (taxonomy) | `HOMDTableData` for all `homd_get_*` table functions |
| `.as_dataframe(engine="pandas")` | Works (no exception) for all three table getters, but only `homd_get_gtdb_taxonomy()` produces a *usable* result. |
| Columns — `homd_get_hmt_lineage()` | `['HOMD.org Taxon Data::Taxonomic Lineage', None]` — **garbage**, do not use. |
| Columns — `homd_get_taxon_table()` | `['<style>']` — **garbage** (matched an HTML page), do not use. |
| Columns — `homd_get_gtdb_taxonomy()` | `GENOME_ID, HMT_ID, strain, hmt_naming_status, hmt_cultivation_status, hmt_primary_body_site_w_abundance, organism, contigs, combined_size, MAG, GC, url, GTDB_taxonomy, bioproject, taxid, biosample, assembly_name, assembly_level, assembly_method, submission_date, geo_loc_name, isolation_source, seqtech, submitter, coverage, ANI, checkM_completeness, checkM_contamination, checkM2_completeness, checkM2_contamination, refseq_assembly, WGS, prokka_CDS, prokka_gene, prokka_mRNA, prokka_misc_RNA, prokka_rRNA, prokka_tRNA, prokka_tmRNA, pangenome` — first column `GENOME_ID` is a genome accession (not a 16S sequence id); `GTDB_taxonomy` column holds the full `d__;p__;c__;...` lineage string. Usable for genome-level taxonomy but is NOT a per-16S-record id — if a per-sequence taxonomy TSV keyed by HMT-ID is required, prefer manually parsing `homd_get_text(<hmt_lineage_url>)` after skipping line 1. |
| `.to_csv(path)` | Raises `NotImplementedError` for all `HOMDTableData` instances (base class not overridden). Must serialize via `.as_dataframe().to_csv(path, sep="\t", index=False)`. |
| Smallest real reference file (for R5) | `https://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/current/HOMD_16S_rRNA_RefSeq_V16.03.fasta` (~10.2 MB) paired with `.../current/HOMD_16S_rRNA_RefSeq_V16.03.qiime.taxonomy` (~0.84 MB). Call: `homd_download_file("ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/current/HOMD_16S_rRNA_RefSeq_V16.03.fasta", dest=...)` and the sibling `...qiime.taxonomy` file the same way. (Discovered via `homd_list_ftp("16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/current")`.) |

### SILVA

| Aspect | Finding |
|---|---|
| Sequence fetch fn | `silva_download_file(path, dest)` / `silva_download_classifier(kind, filename, dest)` |
| Taxonomy fetch fn | **None** — SILVA has no `*TableData` class in 0.4.0 at all (only `SILVAFileListData`, `SILVAReleaseListData`, `SILVATextData`). Taxonomy ships bundled inside the downloaded reference archives (e.g. `tax_slv_ssu_*.txt.gz`, or the `_tax_silva.fasta.gz` headers) — microsuite must parse it itself after download; there is no biodbs-side table object to call `.as_dataframe()`/`.to_csv()` on. |
| Return type (listing) | `silva_list_current_files()` → `SILVAFileListData`, but **returns 0 items** live (base URL is stale — see TL;DR #4). `silva_get_version()` → `SILVATextData`, but its `.text` is the CMS HTML page, not `VERSION.txt` content. |
| `.as_dataframe()` / `.to_csv()` | N/A (no table class); `SILVAFileListData.as_dataframe()` works structurally (`name,url,is_dir` columns) but is empty against the live site. |
| Working base URL | `https://ftp.arb-silva.de/current/` (autoindex still works) vs. biodbs's `https://www.arb-silva.de/current-release/` (CMS page, no autoindex). `VERSION.txt` at the real host is a real 6-byte text file (`138.2`). |
| Smallest real reference file (for R5) | Taxonomy: `https://ftp.arb-silva.de/current/Exports/taxonomy/tax_slv_ssu_138.2.txt.gz` (~211 KB). A larger sequence+taxonomy classifier set lives under `https://ftp.arb-silva.de/current/QIIME2/2025.7/SSU/...` (not sized further — nested one level deeper by region). Because biodbs's listing is broken, R5 will need to pass the **full absolute URL** to `silva_download_file()` (absolute URLs bypass the broken base_url join and download correctly) rather than a relative path discovered via `silva_list_current_files()`. |

### GTDB

| Aspect | Finding |
|---|---|
| Sequence/taxonomy fetch fn | `gtdb_download_taxonomy(domain="bac120"|"ar53", dest, release="latest", compressed=True)`; in-memory: `gtdb_get_taxonomy(domain, release)` |
| Return type (listing) | `gtdb_list_releases()` → `GTDBFileListData` (12 entries: `latest`, `release202`...`release232`, `release80`...`release86`); `gtdb_list_release_files("latest")` → `GTDBFileListData` (20 entries incl. `ar53_taxonomy.tsv[.gz]`, `bac120_taxonomy.tsv[.gz]`) |
| Return type (taxonomy) | `GTDBTableData` |
| `.as_dataframe(engine="pandas")` | **Works cleanly.** `gtdb_get_taxonomy(domain="ar53")` → shape `(22343, 2)`. |
| Columns | `accession, classification` (fixed via the `fieldnames=["accession","classification"]` passed internally by `GTDB_Fetcher.get_taxonomy` — the raw GTDB `.tsv` has no header row, so this hardcoded fieldnames list is required and correct). `classification` is the full `d__Archaea;p__...;s__...` string; first column `accession` is the genome accession (e.g. `RS_GCF_964307495.1`) — this is the "record id" column R2 should treat as column 0 in the emitted TSV. |
| `.to_csv(path)` | Raises `NotImplementedError` (base class not overridden). Must use `.as_dataframe().to_csv(path, sep="\t", index=False)`. |
| Current GTDB release | `latest` resolves to `v232` (per `gtdb_get_version()`). |
| Smallest real reference file (for R5) | `ar53_taxonomy.tsv.gz` (~312 KB compressed; ~3.18 MB decompressed, 22,343 rows) — much smaller than `bac120_taxonomy.tsv.gz` (~9.9 MB). Call: `gtdb_download_taxonomy(domain="ar53", dest=..., release="latest", compressed=True)`. |

**R3 correction (verified live 2026-07-06):** `gtdb_download_file()` is a thin
wrapper around `GTDB_Fetcher.download_file`, which resolves a relative
`path_or_url` via a plain `urljoin` against `base_url =
"https://data.gtdb.ecogenomic.org/releases/"` — it does **not** insert the
release segment for you (unlike `download_taxonomy`, which resolves the full
URL via `_find_release_file`/`list_release_files`). So
`gtdb_download_file("genomic_files_reps/bac120_ssu_reps.fna.gz", dest)`
404s; the release must be included explicitly:
`gtdb_download_file(f"{release}/genomic_files_reps/{domain}_ssu_reps.fna.gz", dest)`.
Confirmed live: dropping `latest/` 404s
(`https://data.gtdb.ecogenomic.org/releases/genomic_files_reps/bac120_ssu_reps.fna.gz`),
while the full path returns 200
(`https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/bac120_ssu_reps.fna.gz`,
~30.8 MB).

Also, R3 originally assumed (per an earlier design note) that the SSU FASTA
header's first token is `{genome_accession}~{contig}...`. **Live sampling of
the real `bac120_ssu_reps.fna.gz` (50 records) shows no `~` at all** — the
header is simply `>{genome_accession} {lineage}[bracketed locus/location
metadata]`, e.g. `>RS_GCF_031457235.1
d__Bacteria;p__Pseudomonadota;...;s__Hydrogenophaga laconesensis
[locus_tag=...] [location=...] [ssu_len=...] [contig_len=...]`, and this
accession matches the taxonomy TSV's first column (`{accession}\t{lineage}`,
headerless) exactly. The implemented adapter still does
`rec_id.split("~", 1)[0]` to derive the taxonomy join key — that is a no-op
when `~` is absent, so it produces correct output against the real live file
without any special-casing, while remaining compatible with any
tilde-suffixed id shape should one appear for a different release/domain.

### GreenGenes

| Aspect | Finding |
|---|---|
| Sequence/taxonomy fetch fn | `greengenes_download_file(path, dest)` — **generic file download only**; there is no `greengenes_get_*` table/text accessor and no `*TableData` class for GreenGenes in 0.4.0 (only `GreenGenesFileListData`/`GreenGenesReleaseListData`). Taxonomy ships as a downloadable file (`.tsv.gz` or bundled in a `.qza`); microsuite must download + parse it itself. |
| Return type (listing) | `greengenes_list_releases()` → `GreenGenesReleaseListData` (12 releases: `2022.7-rc1` .. `2022.10`, `2024.09`, `current`, `gg_12_8`, `gg_12_10`, `gg_13_5`, `gg_13_8_otus`, `unversioned`); `greengenes_list_files("2022.7-rc1")` → `GreenGenesFileListData` (23 files) |
| `.as_dataframe()` / `.to_csv()` | List types only (`name,url,is_dir`/`name,url,modified`); works fine, N/A for taxonomy since there's no table type. |
| Smallest real reference file (for R5) | Full per-ASV taxonomy files are huge (`2022.7.taxonomy.asv.tsv.gz` ≈ 187 MB, `2022.10.taxonomy.asv.tsv.gz` ≈ 188 MB) — avoid these. The much smaller **"backbone"** (representative/clustered) set is the right pick: `2022.7.backbone.v4.fna.qza` (~7.7 MB, V4-region seqs) + `2022.7.backbone.tax.qza` (~2.9 MB, matching taxonomy) under `2022.7-rc1`. Call: `greengenes_download_file("2022.7.backbone.v4.fna.qza", dest=...)` and `greengenes_download_file("2022.7.backbone.tax.qza", dest=...)`. Note: GreenGenes files are `.qza` (QIIME 2 zipped artifacts) — R2/R5 will need microsuite's existing `io/qza.py` reader to pull sequences/taxonomy back out. |

### UNITE

| Aspect | Finding |
|---|---|
| Sequence+taxonomy fetch fn | `unite_download(version, dest, taxon_group="fungi"|"eukaryotes", singletons=False)` — single call downloads one `.tgz`/`.gz` archive containing **both** sequences and taxonomy (UNITE ships them together per RESCRIPt convention); there is no separate taxonomy-table getter and no `*TableData` class for UNITE in 0.4.0 (module only exposes DOI resolution + download). |
| Helper fns | `unite_resolve_doi(version, taxon_group, singletons)` → `str` (DOI, no network call — pure dict lookup); `unite_get_download_url(version, taxon_group, singletons)` → `str` (one cheap PlutoF JSON API call, resolves DOI to the actual media URL) |
| Verified versions | `UNITE_DOIS` table hardcodes 6 releases: `2020-02-20`, `2021-05-10`, `2022-10-16`, `2023-07-18`, `2024-04-04`, `2025-02-19`. `unite_resolve_doi("2025-02-19", "fungi", False)` → `10.15156/BIO/3301241`; `unite_get_download_url(...)` → resolves to a signed S3 URL (`https://s3.hpc.ut.ee/plutof-public/original/<uuid>.tgz`). |
| `.as_dataframe()` / `.to_csv()` | N/A — no table object returned at any point; the archive must be downloaded and extracted, then taxonomy parsed from the extracted files (fasta headers / accompanying `.txt`) by microsuite itself. |
| Smallest real reference file (for R5) | Across the 3 releases checked, size scales with release age: `2020-02-20` fungi non-singleton ≈ 39.2 MB (smallest), `2021-05-10` ≈ 50.4 MB, `2025-02-19` ≈ 78.9 MB. Use the oldest release for the cheapest live test: `unite_download("2020-02-20", dest=..., taxon_group="fungi", singletons=False)`. Still not tiny (~39 MB) — there is no smaller UNITE archive available; consider marking any live UNITE integration test as slow/opt-in. |

### PR2

| Aspect | Finding |
|---|---|
| Sequence+taxonomy fetch fn | `pr2_download_asset(name, dest, tag=None)` — downloads one named GitHub-release asset; sequences and taxonomy are separate assets (paired by filename convention), no single combined table. No `*TableData` class for PR2 in 0.4.0 (only `PR2AssetListData`/`PR2ReleaseListData`). |
| Return type (listing) | `pr2_list_releases()` → `PR2ReleaseListData` (19 releases, newest first: `v5.1.1`, `v5.1.0.0`, `v5.0.0`, ...); `pr2_list_assets()` → `PR2AssetListData` (10 assets for latest release `v5.1.1`) |
| `.as_dataframe()` | Works for both list types. `pr2_list_assets()` columns: `name, url, size` (size in bytes, from GitHub API — no HTTP HEAD needed to size assets, a real cost-saver for R5). |
| `.to_csv()` | Raises `NotImplementedError` (base class). |
| Asset inventory (v5.1.1) | `pr2_version_5.1.1_chimera.xlsx` (1.18 MB, not a taxonomy/seq file), `pr2_version_5.1.1_SSU_mothur.tax.gz` (4.69 MB, taxonomy), `pr2_version_5.1.1_SSU_mothur.fasta.gz` (51.3 MB, sequences — smallest seq file), `pr2_version_5.1.1_SSU_dada2.fasta.gz` (54.8 MB), `pr2_version_5.1.1_SSU_UTAX.fasta.gz` (58.4 MB), `pr2_version_5.1.1_SSU_taxo_long.fasta.gz` (59.6 MB), `pr2_version_5.1.1_taxonomy.xlsx` (3.95 MB, taxonomy only, no matching seq file), `pr2_version_5.1.1_merged.xlsx` (123 MB), `pr2_version_5.1.1_unassigned.xlsx` (128 MB), `pr2_version_5.1.1_emu.zip` (75.2 MB). |
| Smallest real reference file (for R5) | Smallest **matched seq+tax pair**: mothur format — `pr2_version_5.1.1_SSU_mothur.fasta.gz` (~51.3 MB) + `pr2_version_5.1.1_SSU_mothur.tax.gz` (~4.69 MB). Calls: `pr2_download_asset("pr2_version_5.1.1_SSU_mothur.fasta.gz", dest=...)` and `pr2_download_asset("pr2_version_5.1.1_SSU_mothur.tax.gz", dest=...)`. (`.xlsx`/`.zip` assets are not fasta/plain-text taxonomy and are out of scope for a refdb provider.) |

---

## Serialization rule for R2's `_write_taxonomy_tsv`

Because **no** `*TableData` subclass in biodbs 0.4.0 overrides `to_csv`, R2 must
never call `<biodbs_table_obj>.to_csv(...)`. The only reliable path is:

```python
df = biodbs_table_obj.as_dataframe(engine="pandas")
df.to_csv(dest, sep="\t", index=False)
```

And even `.as_dataframe()` is only trustworthy for `homd_get_gtdb_taxonomy()`
and `gtdb_get_taxonomy()` among the taxonomy-table-shaped calls probed here —
HOMD's other two table getters (`homd_get_hmt_lineage`, `homd_get_taxon_table`)
return structurally valid but semantically garbage DataFrames and must not be
used as-is. SILVA, GreenGenes, UNITE, and PR2 have **no** taxonomy table
object at all in 0.4.0 — for those four, taxonomy always arrives as a
downloaded file (fasta headers, `.txt`/`.tsv` archive members, or an `.xlsx`)
that microsuite's own provider code must parse after `*_download_*`/
`*_download_asset` completes.
