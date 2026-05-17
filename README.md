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

| Task | Backends | Status | Purpose |
| --- | --- | --- | --- |
| `qc` | `fastqc`, `multiqc`, `qiime2-demux` | partial | Raw-read and demultiplexed-read quality reports. |
| `trim` | `fastp`; planned: `cutadapt`, `qiime2-cutadapt` | partial | Adapter trimming, quality filtering, and fastp reports. |
| `denoise` | `qiime2-dada2`, `qiime2-deblur`; planned: `dada2-r` | partial | Amplicon denoising from demultiplexed QIIME 2 artifacts. |
| `cluster` | `vsearch` | partial | QIIME 2 VSEARCH de novo feature clustering. |
| `tax_classify` | `qiime2`, `kraken2`, `bracken`, `metaphlan` | partial | Taxonomy assignment or taxonomic profiling. |
| `phylogeny` | `qiime2`, `mafft-fasttree` | partial | Alignment, masking, tree construction, and rooting. |
| `normalize` | `native` | ready | Relative abundance, CLR, and related table transforms. |
| `abundance` | `native` | ready | Summarize abundance at taxonomy levels. |
| `shared_taxa` | `native` | ready | Compare shared taxa across sample groups. |
| `rarefy` | `native` | ready | Rarefy feature tables to a fixed depth. |
| `diversity_calc` | `native`, `qiime2` | partial | Alpha/beta diversity calculation. |
| `beta_significance` | `qiime2`, `native` | planned | PERMANOVA and related beta-diversity tests. |
| `diff_abundance` | `ancombc`, `aldex2`, `maaslin2`, `lefse` | partial | Differential abundance testing. |
| `env_assoc` | `mantel`, `rda`, `cca`, `db-rda` | planned | Environmental association and constrained ordination. |
| `network` | `native-correlation`, `sparcc`, `spieceasi`, `flashweave` | planned | Taxa association network analysis. |
| `functional_predict` | `picrust2`, `tax4fun2` | planned | Predict function from marker-gene profiles. |
| `functional_profile` | `humann` | planned | Functional profiling from metagenomic data. |
| `classify_samples` | `randomforest`, `xgboost` | planned | Supervised sample classification. |
| `time_series` | `native` | planned | Longitudinal microbiome analysis. |
| `gamma_diversity` | `native` | planned | Region/group-level diversity summaries. |
| `turnover` | `beta-turnover`, `taxa-turnover` | planned | Community and taxa turnover analysis. |
| `visualize` | `native` | ready | Barplots, ordination plots, and heatmaps. |
| `report` | `native` | ready | HTML provenance reports from run metadata. |

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
