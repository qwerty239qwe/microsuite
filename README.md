# microsuite

`microsuite` is a multi-environment microbiome toolbox with three APIs:

- Nextflow API for reproducible full workflows.
- CLI API for ergonomic one-step task commands.
- Python SDK for programmatic table/statistics functions.

- run built-in workflows from feature tables or QIIME 2 artifacts
- run method-oriented tasks such as taxonomy classification
- import TSV, BIOM, or QIIME 2-compatible `.qza` feature tables into AnnData
- compute alpha and beta diversity
- run PCoA
- draw taxonomy barplots
- fetch and run demo datasets

See [docs/three-api-roadmap.md](docs/three-api-roadmap.md) for the architecture.
Demo data attribution and citation details are in
[docs/data-attribution.md](docs/data-attribution.md).

## Method Surface

### Quality Reports

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `fastqc` | FastQC 0.12.1 | ready | `microsuite qc --backend fastqc` | `qc(backend="fastqc", inputs=[...], output_dir=...)` | [FastQC](containers/fastqc/Dockerfile) or external `fastqc` | Ready as a CLI wrapper and standalone container; Nextflow raw-read wiring remains planned. | Raw-read quality reports. |
| `multiqc` | User env | partial | `microsuite qc --backend multiqc` | `qc(backend="multiqc", input_dir=..., output_dir=...)` | External `multiqc`; container planned | Good aggregation layer; depends on upstream report files. | Aggregate QC reports. |
| `qiime2-demux-summarize` | QIIME 2 2024.10 | partial | `microsuite qc --backend qiime2-demux` | `qc(backend="qiime2-demux", demux=..., output=...)` | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Best for QIIME artifacts; current CLI backend name should be renamed or aliased. | Demultiplexed-read quality visualization. |

### Quality Filtering

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `qiime2-exclude-seqs` | QIIME 2 2026.4 API | partial | `microsuite qc_filter --backend qiime2-exclude-seqs` | `qc_filter(backend="qiime2-exclude-seqs", query_sequences=..., reference_sequences=...)` | External QIIME 2 with `q2-quality-control`; [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) pin update pending | Strong contaminant/non-target sequence filtering; needs reference sequences and threshold choices. | Exclude or retain feature sequences by alignment to reference sequences. |
| `qiime2-filter-reads` | QIIME 2 2026.4 API | partial | `microsuite qc_filter --backend qiime2-filter-reads` | `qc_filter(backend="qiime2-filter-reads", demux=..., database=..., output=...)` | External QIIME 2 with `q2-quality-control`; [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) pin update pending | Useful for host/contaminant read removal; requires a Bowtie2 index. | Filter demultiplexed reads by alignment to a reference database. |
| `qiime2-bowtie2-build` | QIIME 2 2026.4 API | partial | `microsuite qc_filter --backend qiime2-bowtie2-build` | `qc_filter(backend="qiime2-bowtie2-build", sequences=..., output=...)` | External QIIME 2 with `q2-quality-control`; [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) pin update pending | Completes the filter-reads setup path; still requires suitable reference sequences. | Build a Bowtie2 index artifact for read filtering. |
| `qiime2-decontam` | QIIME 2 2026.4 API | partial | `microsuite decontam --backend qiime2-decontam` | `decontam(backend="qiime2-decontam", table=..., metadata=..., output=...)` | External QIIME 2 with `q2-quality-control`; [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) pin update pending | Useful contamination screening; requires negative controls or concentration metadata. | Identify likely contaminant features with decontam. |
| `qiime2-quality-filter-q-score` | QIIME 2 2024.10 | planned | planned | planned | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Native QIIME quality-score filtering; mainly useful before downstream QIIME artifact workflows. | Filter demultiplexed reads by quality scores. |

### Trimming

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `fastp` | User env | partial | `microsuite trim --backend fastp` | `trim(backend="fastp", read1=..., output1=...)` | External `fastp`; container planned | Fast all-in-one preprocessing; primer-specific trimming is less explicit than Cutadapt. | Adapter trimming, quality filtering, HTML/JSON reports. |
| `cutadapt` | Cutadapt >=4.x user env | partial | `microsuite trim --backend cutadapt` | `trim(backend="cutadapt", read1=..., output1=..., adapter=...)` | External `cutadapt` on `PATH`; container planned | Precise primer/adaptor trimming with explicit adapter control; requires users to choose primer/adapter sequences. | Adapter/primer trimming and read filtering. |
| `trimmomatic` | Trimmomatic >=0.39 user env | partial | `microsuite trim --backend trimmomatic` | `trim(backend="trimmomatic", read1=..., output1=..., trimmomatic_steps=[...])` | External `trimmomatic` on `PATH`; container planned | Mature Java trimmer with explicit step pipeline; paired mode requires unpaired output files. | Sliding-window, length, quality, and adapter trimming. |
| `trim-galore` | Trim Galore >=0.6 user env | partial | `microsuite trim --backend trim-galore` | `trim(backend="trim-galore", read1=..., output1=..., adapter=...)` | External `trim_galore` on `PATH`; container planned | Convenient Cutadapt/FastQC wrapper; output names are mostly tool-controlled. | Adapter/quality trimming with integrated QC conventions. |
| `qiime2-cutadapt` | QIIME 2 2024.10 | planned | `microsuite trim --backend qiime2-cutadapt` | planned | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Fits QIIME artifact workflows; less convenient for raw FASTQ-only runs. | QIIME 2 Cutadapt wrapper. |

### Denoising And Clustering

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `qiime2-dada2` | QIIME 2 2024.10 | partial | `microsuite denoise --backend qiime2-dada2` | `denoise(backend="qiime2-dada2", demux=..., output_table=...)` | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Strong amplicon default; needs careful truncation choices. | DADA2 ASV inference from demultiplexed reads. |
| `qiime2-deblur` | QIIME 2 2024.10 | partial | `microsuite denoise --backend qiime2-deblur` | `denoise(backend="qiime2-deblur", demux=..., output_table=...)` | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Reproducible fixed-error model; mainly 16S-oriented. | Deblur ASV inference from demultiplexed reads. |
| `dada2-r` | DADA2 R user env | partial | `microsuite denoise --backend dada2-r` | `denoise(backend="dada2-r", demux=reads_dir, output_table=...)` | External `Rscript` with R package `dada2`; dedicated DADA2 image later | Direct R ecosystem access; currently expects a FASTQ directory and writes TSV/FASTA outputs. | R/DADA2 ASV inference from raw or trimmed FASTQ files. |
| `vsearch` | QIIME 2 2024.10 | partial | `microsuite cluster --backend vsearch` | `cluster(backend="vsearch", table=..., rep_seqs=...)` | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Useful for OTU workflows; less ASV-centric. | QIIME 2 VSEARCH de novo feature clustering. |

### Taxonomy And Phylogeny

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `qiime2` | QIIME 2 2024.10 | partial | `microsuite tax_classify --backend qiime2` | `tax_classify(backend="qiime2", rep_seqs=..., classifier=..., output=...)` | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Strong classifier ecosystem; requires trained classifier artifacts. | QIIME 2 taxonomy classification. |
| `qiime2-taxonomy` | QIIME 2 2026.4 API | partial | `microsuite evaluate --backend qiime2-taxonomy` | `evaluate(backend="qiime2-taxonomy", expected_taxa=..., observed_taxa=...)` | External QIIME 2 with `q2-quality-control`; [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) pin update pending | Best for mock-community or known-composition validation; needs expected taxonomy. | Evaluate observed taxonomy against expected assignments. |
| `kraken2` | Kraken2 2.1.3 | planned | `microsuite tax_classify --backend kraken2` | planned | [Kraken2](containers/kraken2/Dockerfile) | Fast profiling; requires large databases. | Taxonomic profiling/classification. |
| `bracken` | Planned | planned | `microsuite tax_classify --backend bracken` | planned | [Kraken2](containers/kraken2/Dockerfile), Bracken planned | Improves abundance estimates; depends on Kraken2 database setup. | Abundance re-estimation from Kraken2 output. |
| `metaphlan` | Planned | planned | `microsuite tax_classify --backend metaphlan` | planned | Image not added yet | Good marker-gene profiling; separate database/runtime needed. | Marker-gene taxonomic profiling. |
| `qiime2-phylogeny` | QIIME 2 2024.10 | partial | planned | planned | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Integrated with QIIME artifacts; heavier runtime. | Alignment, masking, tree construction, and rooting. |
| `mafft-fasttree` | Planned | planned | planned | planned | Image not added yet | Lightweight standalone path; more file-format handling needed. | Standalone MAFFT/FastTree phylogeny. |

### Table Transforms And Summaries

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `native-normalize` | microsuite 0.1.0 | ready | `microsuite normalize --backend native` | `normalize_table(adata, method="relative")` | [microsuite Python](containers/microsuite/Dockerfile) | Fast and portable; narrower than specialized compositional packages. | Relative abundance, CLR, and table transforms. |
| `native-abundance` | microsuite 0.1.0 | ready | `microsuite abundance --backend native` | `abundance_table(adata, level="genus")` | [microsuite Python](containers/microsuite/Dockerfile) | Simple summary output; depends on taxonomy quality. | Summarize abundance at taxonomy levels. |
| `native-shared-taxa` | microsuite 0.1.0 | ready | `microsuite shared_taxa --backend native` | `shared_taxa_table(adata, level="genus", group=...)` | [microsuite Python](containers/microsuite/Dockerfile) | Easy group comparison; descriptive rather than inferential. | Compare shared taxa across sample groups. |
| `native-rarefy` | microsuite 0.1.0 | ready | `microsuite rarefy --backend native` | `rarefy_table(adata, depth=...)` | [microsuite Python](containers/microsuite/Dockerfile) | Reproducible with seeds; discards reads by design. | Rarefy feature tables to a fixed depth. |

### Diversity And Ecological Statistics

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `native` | microsuite 0.1.0 | ready | lower-level `microsuite diversity ...` | `alpha_diversity(adata, metric=...)` or `beta_diversity(adata, metric=...)` | [microsuite Python](containers/microsuite/Dockerfile) | Lightweight and Windows-friendly; phylogenetic metrics are limited. | Native alpha/beta diversity. |
| `qiime2-diversity-lib` | QIIME 2 2024.10 | partial | `microsuite diversity_calc --backend qiime2` | `diversity_calc(backend="qiime2", metric=..., table=..., output=...)` | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Broader metric coverage; requires QIIME artifacts. | QIIME 2 diversity-lib metrics. |
| `qiime2-beta-significance` | QIIME 2 2024.10 | planned | planned | planned | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Established PERMANOVA workflow; artifact-based. | PERMANOVA and beta-diversity tests. |
| `native-beta-significance` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) | Easier SDK/web integration; needs statistical validation. | Native beta-diversity tests. |
| `mantel` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) or R image later | Useful distance association; sensitive to design assumptions. | Mantel association testing. |
| `rda` | Planned | planned | planned | planned | R image later | Interpretable constrained ordination; needs careful preprocessing. | Redundancy analysis. |
| `cca` | Planned | planned | planned | planned | R image later | Handles unimodal ecological gradients; requires ecological assumptions. | Canonical correspondence analysis. |
| `db-rda` | Planned | planned | planned | planned | R image later | Works with distance matrices; permutation design matters. | Distance-based redundancy analysis. |
| `native-gamma-diversity` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) | Useful group-level summary; definitions must be explicit. | Region/group-level diversity summaries. |
| `beta-turnover` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) | Good ecological interpretation; metric choices matter. | Community turnover analysis. |
| `taxa-turnover` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) | Taxon-centric interpretation; sensitive to filtering. | Taxa turnover analysis. |

### Differential Abundance

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `ancombc` | R 4.4.0 image; ANCOMBC via Bioconductor | partial | `microsuite diff_abundance --backend ancombc` | `diff_abundance(backend="ancombc", table=..., group=..., output=...)` | [R diffab](containers/r-diffab/Dockerfile) | Strong compositional method; R/Bioconductor install is heavier. | ANCOM-BC differential abundance. |
| `aldex2` | Planned | planned | `microsuite diff_abundance --backend aldex2` | planned | [R diffab](containers/r-diffab/Dockerfile) | Good compositional alternative; R wrapper pending. | ALDEx2 differential abundance. |
| `maaslin2` | Planned | planned | `microsuite diff_abundance --backend maaslin2` | planned | [R diffab](containers/r-diffab/Dockerfile) | Flexible covariate modeling; more complex formula interface. | MaAsLin2 multivariable association testing. |
| `lefse` | Planned | planned | `microsuite diff_abundance --backend lefse` | planned | Image not added yet | Familiar legacy workflow; weaker modern compositional assumptions. | LEfSe legacy differential abundance. |

### Networks, Function, Machine Learning, And Reporting

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `native-correlation` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) | Simple and transparent; compositional bias risk. | Correlation network analysis. |
| `sparcc` | Planned | planned | planned | planned | Image not added yet | Compositional-aware; runtime/database setup pending. | SparCC association network analysis. |
| `spieceasi` | Planned | planned | planned | planned | R image later | Strong ecological network method; R dependency and tuning burden. | SPIEC-EASI network inference. |
| `flashweave` | Planned | planned | planned | planned | Image not added yet | Handles heterogeneous metadata; separate runtime needed. | FlashWeave network inference. |
| `picrust2` | Planned | planned | planned | planned | Image not added yet | Popular marker-gene function prediction; reference-dependent. | Predict function from marker-gene profiles. |
| `tax4fun2` | Planned | planned | planned | planned | Image not added yet | Alternative functional prediction; separate runtime needed. | Tax4Fun2 function prediction. |
| `humann` | Planned | planned | planned | planned | Image not added yet | Strong metagenomic functional profiling; heavy databases. | Functional profiling from metagenomic data. |
| `randomforest` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) | Interpretable baseline ML; needs validation/splits. | Supervised sample classification. |
| `xgboost` | Planned | planned | planned | planned | Image not added yet | Strong predictive model; optional dependency and tuning needed. | Optional XGBoost sample classification. |
| `native-time-series` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) | Web-friendly outputs; design choices still open. | Longitudinal microbiome analysis. |
| `native-visualize` | microsuite 0.1.0 | ready | `microsuite viz ...` | `taxonomy_barplot(adata, level=..., output=...)` | [microsuite Python](containers/microsuite/Dockerfile) | Lightweight static figures; not an interactive dashboard. | Barplots, ordination plots, and heatmaps. |
| `native-report` | microsuite 0.1.0 | ready | `microsuite report --backend native` | `report(backend="native", run_dir=..., output=...)` | [microsuite Python](containers/microsuite/Dockerfile) | Good provenance summary; not full narrative reporting yet. | HTML provenance reports from run metadata. |

## Install

```bash
uv sync --extra dev
```

Optional compatibility extras:

```bash
uv sync --extra biom --extra dev
uv sync --extra qza --extra dev
uv sync --extra all --extra dev
```

## CLI

```bash
microsuite workflow list
microsuite workflow moving-pictures --out runs/moving-pictures --force
microsuite workflow table-summary \
  --format tsv \
  --table table.tsv \
  --metadata metadata.tsv \
  --taxonomy taxonomy.tsv \
  --out runs/table-summary
microsuite qc --backend fastqc --input sample_R1.fastq.gz --output-dir qc/fastqc
microsuite qc --backend multiqc --input-dir qc/fastqc --output-dir qc/multiqc
microsuite qc --backend qiime2-demux --demux demux.qza -o qc/demux.qzv
microsuite trim \
  --backend fastp \
  --read1 sample_R1.fastq.gz \
  --read2 sample_R2.fastq.gz \
  --output1 trimmed_R1.fastq.gz \
  --output2 trimmed_R2.fastq.gz \
  --html qc/fastp.html \
  --json-report qc/fastp.json
microsuite denoise \
  --backend qiime2-dada2 \
  --demux demux.qza \
  --output-table table.qza \
  --output-rep-seqs rep-seqs.qza \
  --output-stats stats.qza \
  --trunc-len 150
microsuite cluster \
  --backend vsearch \
  --table table.qza \
  --rep-seqs rep-seqs.qza \
  --output-table clustered-table.qza \
  --output-rep-seqs clustered-rep-seqs.qza \
  --identity 0.97
microsuite normalize --backend native --method relative --table table.h5ad -o relative.h5ad
microsuite abundance --backend native --table table.h5ad --level genus -o abundance.tsv
microsuite shared_taxa --backend native --table table.h5ad --level genus --group body_site -o shared.tsv
microsuite rarefy --backend native --table table.h5ad --depth 10000 -o rarefied.h5ad
microsuite methods
microsuite tax_classify \
  --backend qiime2 \
  --rep-seqs rep-seqs.qza \
  --classifier classifier.qza \
  -o taxonomy.qza
microsuite diversity_calc \
  --backend qiime2 \
  --metric bray-curtis \
  --table table.qza \
  -o bray-curtis.qza
microsuite diff_abundance --backend ancombc --table table.h5ad --group treatment -o diff.tsv
microsuite report --backend native --run-dir runs/table-summary -o report.html

microsuite import tsv table.tsv --metadata metadata.tsv --taxonomy taxonomy.tsv -o table.h5ad
microsuite diversity alpha table.h5ad --metric shannon -o alpha.tsv
microsuite diversity beta table.h5ad --metric bray-curtis -o beta.tsv
microsuite ordination pcoa beta.tsv -o pcoa.tsv
microsuite viz barplot table.h5ad --level genus -o barplot.png
microsuite data fetch moving-pictures -o data/moving-pictures-real --full
microsuite qiime inspect data/moving-pictures-real/table.qza
microsuite qiime extract data/moving-pictures-real/taxonomy.qza -o data/taxonomy-extract
```

Use method-oriented commands such as `tax_classify` when you know the task you
want to run. Use `workflow` commands for complete pipelines. Use `import`,
`diversity`, `ordination`, `viz`, and `qiime` as lower-level building blocks.
Commands overwrite outputs only when `--force` is supplied.

## Development

```bash
uv run ruff check .
uv run ruff format --check .
uv run ty check
uv run pytest
uv build
```

GitHub Actions runs the same quality gate on Python 3.11 and 3.12.
