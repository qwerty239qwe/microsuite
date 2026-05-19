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

| Task | Backend | Status | Image / environment | Purpose |
| --- | --- | --- | --- | --- |
| `qc` | `fastqc` | partial | External `fastqc`; container planned | Raw-read quality reports. |
| `qc` | `multiqc` | partial | External `multiqc`; container planned | Aggregate QC reports. |
| `qc` | `qiime2-demux` | partial | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Demultiplexed-read quality visualization. |
| `trim` | `fastp` | partial | External `fastp`; container planned | Adapter trimming, quality filtering, HTML/JSON reports. |
| `trim` | `cutadapt` | planned | Image not added yet | Adapter/primer trimming. |
| `trim` | `qiime2-cutadapt` | planned | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | QIIME 2 Cutadapt wrapper. |
| `denoise` | `qiime2-dada2` | partial | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | DADA2 ASV inference from demultiplexed reads. |
| `denoise` | `qiime2-deblur` | partial | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Deblur ASV inference from demultiplexed reads. |
| `denoise` | `dada2-r` | planned | [R diffab](containers/r-diffab/Dockerfile) or dedicated DADA2 image later | R/DADA2 denoising. |
| `cluster` | `vsearch` | partial | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | QIIME 2 VSEARCH de novo feature clustering. |
| `tax_classify` | `qiime2` | partial | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | QIIME 2 taxonomy classification. |
| `tax_classify` | `kraken2` | planned | [Kraken2](containers/kraken2/Dockerfile) | Taxonomic profiling/classification. |
| `tax_classify` | `bracken` | planned | [Kraken2](containers/kraken2/Dockerfile), Bracken planned | Abundance re-estimation from Kraken2 output. |
| `tax_classify` | `metaphlan` | planned | Image not added yet | Marker-gene taxonomic profiling. |
| `phylogeny` | `qiime2` | partial | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Alignment, masking, tree construction, and rooting. |
| `phylogeny` | `mafft-fasttree` | planned | Image not added yet | Standalone MAFFT/FastTree phylogeny. |
| `normalize` | `native` | ready | [microsuite Python](containers/microsuite/Dockerfile) | Relative abundance, CLR, and table transforms. |
| `abundance` | `native` | ready | [microsuite Python](containers/microsuite/Dockerfile) | Summarize abundance at taxonomy levels. |
| `shared_taxa` | `native` | ready | [microsuite Python](containers/microsuite/Dockerfile) | Compare shared taxa across sample groups. |
| `rarefy` | `native` | ready | [microsuite Python](containers/microsuite/Dockerfile) | Rarefy feature tables to a fixed depth. |
| `diversity_calc` | `native` | ready | [microsuite Python](containers/microsuite/Dockerfile) | Native alpha/beta diversity. |
| `diversity_calc` | `qiime2` | partial | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | QIIME 2 diversity-lib metrics. |
| `beta_significance` | `qiime2` | planned | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | PERMANOVA and beta-diversity tests. |
| `beta_significance` | `native` | planned | [microsuite Python](containers/microsuite/Dockerfile) | Native beta-diversity tests. |
| `diff_abundance` | `ancombc` | partial | [R diffab](containers/r-diffab/Dockerfile) | ANCOM-BC differential abundance. |
| `diff_abundance` | `aldex2` | planned | [R diffab](containers/r-diffab/Dockerfile) | ALDEx2 differential abundance. |
| `diff_abundance` | `maaslin2` | planned | [R diffab](containers/r-diffab/Dockerfile) | MaAsLin2 multivariable association testing. |
| `diff_abundance` | `lefse` | planned | Image not added yet | LEfSe legacy differential abundance. |
| `env_assoc` | `mantel` | planned | [microsuite Python](containers/microsuite/Dockerfile) or R image later | Mantel association testing. |
| `env_assoc` | `rda` | planned | R image later | Redundancy analysis. |
| `env_assoc` | `cca` | planned | R image later | Canonical correspondence analysis. |
| `env_assoc` | `db-rda` | planned | R image later | Distance-based redundancy analysis. |
| `network` | `native-correlation` | planned | [microsuite Python](containers/microsuite/Dockerfile) | Correlation network analysis. |
| `network` | `sparcc` | planned | Image not added yet | SparCC association network analysis. |
| `network` | `spieceasi` | planned | R image later | SPIEC-EASI network inference. |
| `network` | `flashweave` | planned | Image not added yet | FlashWeave network inference. |
| `functional_predict` | `picrust2` | planned | Image not added yet | Predict function from marker-gene profiles. |
| `functional_predict` | `tax4fun2` | planned | Image not added yet | Tax4Fun2 function prediction. |
| `functional_profile` | `humann` | planned | Image not added yet | Functional profiling from metagenomic data. |
| `classify_samples` | `randomforest` | planned | [microsuite Python](containers/microsuite/Dockerfile) | Supervised sample classification. |
| `classify_samples` | `xgboost` | planned | Image not added yet | Optional XGBoost sample classification. |
| `time_series` | `native` | planned | [microsuite Python](containers/microsuite/Dockerfile) | Longitudinal microbiome analysis. |
| `gamma_diversity` | `native` | planned | [microsuite Python](containers/microsuite/Dockerfile) | Region/group-level diversity summaries. |
| `turnover` | `beta-turnover` | planned | [microsuite Python](containers/microsuite/Dockerfile) | Community turnover analysis. |
| `turnover` | `taxa-turnover` | planned | [microsuite Python](containers/microsuite/Dockerfile) | Taxa turnover analysis. |
| `visualize` | `native` | ready | [microsuite Python](containers/microsuite/Dockerfile) | Barplots, ordination plots, and heatmaps. |
| `report` | `native` | ready | [microsuite Python](containers/microsuite/Dockerfile) | HTML provenance reports from run metadata. |

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
