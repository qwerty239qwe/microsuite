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

| Subtopic | Backend | Version | Status | CLI command | Python function name | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `qc` | `fastqc` | User env | partial | `microsuite qc --backend fastqc` | `microsuite.methods.qc.qc` | External `fastqc`; container planned | Fast and familiar; requires separate install. | Raw-read quality reports. |
| `qc` | `multiqc` | User env | partial | `microsuite qc --backend multiqc` | `microsuite.methods.qc.qc` | External `multiqc`; container planned | Good aggregation layer; depends on upstream report files. | Aggregate QC reports. |
| `qc` | `qiime2-demux` | QIIME 2 2024.10 | partial | `microsuite qc --backend qiime2-demux` | `microsuite.methods.qc.qc` | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Best for QIIME artifacts; heavier runtime. | Demultiplexed-read quality visualization. |
| `trim` | `fastp` | User env | partial | `microsuite trim --backend fastp` | `microsuite.methods.trim.trim` | External `fastp`; container planned | Fast all-in-one preprocessing; primer-specific trimming is less explicit than Cutadapt. | Adapter trimming, quality filtering, HTML/JSON reports. |
| `trim` | `cutadapt` | Planned | planned | `microsuite trim --backend cutadapt` | planned | Image not added yet | Precise primer/adaptor trimming; more parameters to expose. | Adapter/primer trimming. |
| `trim` | `qiime2-cutadapt` | QIIME 2 2024.10 | planned | `microsuite trim --backend qiime2-cutadapt` | planned | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Fits QIIME artifact workflows; less convenient for raw FASTQ-only runs. | QIIME 2 Cutadapt wrapper. |
| `denoise` | `qiime2-dada2` | QIIME 2 2024.10 | partial | `microsuite denoise --backend qiime2-dada2` | `microsuite.methods.denoise.denoise` | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Strong amplicon default; needs careful truncation choices. | DADA2 ASV inference from demultiplexed reads. |
| `denoise` | `qiime2-deblur` | QIIME 2 2024.10 | partial | `microsuite denoise --backend qiime2-deblur` | `microsuite.methods.denoise.denoise` | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Reproducible fixed-error model; mainly 16S-oriented. | Deblur ASV inference from demultiplexed reads. |
| `denoise` | `dada2-r` | Planned | planned | `microsuite denoise --backend dada2-r` | planned | [R diffab](containers/r-diffab/Dockerfile) or dedicated DADA2 image later | Direct R ecosystem access; needs a dedicated DADA2 runtime. | R/DADA2 denoising. |
| `cluster` | `vsearch` | QIIME 2 2024.10 | partial | `microsuite cluster --backend vsearch` | `microsuite.methods.cluster.cluster` | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Useful for OTU workflows; less ASV-centric. | QIIME 2 VSEARCH de novo feature clustering. |
| `tax_classify` | `qiime2` | QIIME 2 2024.10 | partial | `microsuite tax_classify --backend qiime2` | `microsuite.methods.tax_classify.tax_classify` | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Strong classifier ecosystem; requires trained classifier artifacts. | QIIME 2 taxonomy classification. |
| `tax_classify` | `kraken2` | Kraken2 2.1.3 | planned | `microsuite tax_classify --backend kraken2` | planned | [Kraken2](containers/kraken2/Dockerfile) | Fast profiling; requires large databases. | Taxonomic profiling/classification. |
| `tax_classify` | `bracken` | Planned | planned | `microsuite tax_classify --backend bracken` | planned | [Kraken2](containers/kraken2/Dockerfile), Bracken planned | Improves abundance estimates; depends on Kraken2 database setup. | Abundance re-estimation from Kraken2 output. |
| `tax_classify` | `metaphlan` | Planned | planned | `microsuite tax_classify --backend metaphlan` | planned | Image not added yet | Good marker-gene profiling; separate database/runtime needed. | Marker-gene taxonomic profiling. |
| `phylogeny` | `qiime2` | QIIME 2 2024.10 | partial | planned | planned | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Integrated with QIIME artifacts; heavier runtime. | Alignment, masking, tree construction, and rooting. |
| `phylogeny` | `mafft-fasttree` | Planned | planned | planned | planned | Image not added yet | Lightweight standalone path; more file-format handling needed. | Standalone MAFFT/FastTree phylogeny. |
| `normalize` | `native` | microsuite 0.1.0 | ready | `microsuite normalize --backend native` | `microsuite.api.normalize_table` | [microsuite Python](containers/microsuite/Dockerfile) | Fast and portable; narrower than specialized compositional packages. | Relative abundance, CLR, and table transforms. |
| `abundance` | `native` | microsuite 0.1.0 | ready | `microsuite abundance --backend native` | `microsuite.api.abundance_table` | [microsuite Python](containers/microsuite/Dockerfile) | Simple summary output; depends on taxonomy quality. | Summarize abundance at taxonomy levels. |
| `shared_taxa` | `native` | microsuite 0.1.0 | ready | `microsuite shared_taxa --backend native` | `microsuite.api.shared_taxa_table` | [microsuite Python](containers/microsuite/Dockerfile) | Easy group comparison; descriptive rather than inferential. | Compare shared taxa across sample groups. |
| `rarefy` | `native` | microsuite 0.1.0 | ready | `microsuite rarefy --backend native` | `microsuite.api.rarefy_table` | [microsuite Python](containers/microsuite/Dockerfile) | Reproducible with seeds; discards reads by design. | Rarefy feature tables to a fixed depth. |
| `diversity_calc` | `native` | microsuite 0.1.0 | ready | lower-level `microsuite diversity ...` | `microsuite.api.alpha_diversity`, `microsuite.api.beta_diversity` | [microsuite Python](containers/microsuite/Dockerfile) | Lightweight and Windows-friendly; phylogenetic metrics are limited. | Native alpha/beta diversity. |
| `diversity_calc` | `qiime2` | QIIME 2 2024.10 | partial | `microsuite diversity_calc --backend qiime2` | `microsuite.methods.diversity_calc.diversity_calc` | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Broader metric coverage; requires QIIME artifacts. | QIIME 2 diversity-lib metrics. |
| `beta_significance` | `qiime2` | QIIME 2 2024.10 | planned | planned | planned | [QIIME 2 amplicon](containers/qiime2-amplicon/Dockerfile) | Established PERMANOVA workflow; artifact-based. | PERMANOVA and beta-diversity tests. |
| `beta_significance` | `native` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) | Easier SDK/web integration; needs statistical validation. | Native beta-diversity tests. |
| `diff_abundance` | `ancombc` | R 4.4.0 image; ANCOMBC via Bioconductor | partial | `microsuite diff_abundance --backend ancombc` | `microsuite.methods.diff_abundance.diff_abundance` | [R diffab](containers/r-diffab/Dockerfile) | Strong compositional method; R/Bioconductor install is heavier. | ANCOM-BC differential abundance. |
| `diff_abundance` | `aldex2` | Planned | planned | `microsuite diff_abundance --backend aldex2` | planned | [R diffab](containers/r-diffab/Dockerfile) | Good compositional alternative; R wrapper pending. | ALDEx2 differential abundance. |
| `diff_abundance` | `maaslin2` | Planned | planned | `microsuite diff_abundance --backend maaslin2` | planned | [R diffab](containers/r-diffab/Dockerfile) | Flexible covariate modeling; more complex formula interface. | MaAsLin2 multivariable association testing. |
| `diff_abundance` | `lefse` | Planned | planned | `microsuite diff_abundance --backend lefse` | planned | Image not added yet | Familiar legacy workflow; weaker modern compositional assumptions. | LEfSe legacy differential abundance. |
| `env_assoc` | `mantel` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) or R image later | Useful distance association; sensitive to design assumptions. | Mantel association testing. |
| `env_assoc` | `rda` | Planned | planned | planned | planned | R image later | Interpretable constrained ordination; needs careful preprocessing. | Redundancy analysis. |
| `env_assoc` | `cca` | Planned | planned | planned | planned | R image later | Handles unimodal ecological gradients; requires ecological assumptions. | Canonical correspondence analysis. |
| `env_assoc` | `db-rda` | Planned | planned | planned | planned | R image later | Works with distance matrices; permutation design matters. | Distance-based redundancy analysis. |
| `network` | `native-correlation` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) | Simple and transparent; compositional bias risk. | Correlation network analysis. |
| `network` | `sparcc` | Planned | planned | planned | planned | Image not added yet | Compositional-aware; runtime/database setup pending. | SparCC association network analysis. |
| `network` | `spieceasi` | Planned | planned | planned | planned | R image later | Strong ecological network method; R dependency and tuning burden. | SPIEC-EASI network inference. |
| `network` | `flashweave` | Planned | planned | planned | planned | Image not added yet | Handles heterogeneous metadata; separate runtime needed. | FlashWeave network inference. |
| `functional_predict` | `picrust2` | Planned | planned | planned | planned | Image not added yet | Popular marker-gene function prediction; reference-dependent. | Predict function from marker-gene profiles. |
| `functional_predict` | `tax4fun2` | Planned | planned | planned | planned | Image not added yet | Alternative functional prediction; separate runtime needed. | Tax4Fun2 function prediction. |
| `functional_profile` | `humann` | Planned | planned | planned | planned | Image not added yet | Strong metagenomic functional profiling; heavy databases. | Functional profiling from metagenomic data. |
| `classify_samples` | `randomforest` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) | Interpretable baseline ML; needs validation/splits. | Supervised sample classification. |
| `classify_samples` | `xgboost` | Planned | planned | planned | planned | Image not added yet | Strong predictive model; optional dependency and tuning needed. | Optional XGBoost sample classification. |
| `time_series` | `native` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) | Web-friendly outputs; design choices still open. | Longitudinal microbiome analysis. |
| `gamma_diversity` | `native` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) | Useful group-level summary; definitions must be explicit. | Region/group-level diversity summaries. |
| `turnover` | `beta-turnover` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) | Good ecological interpretation; metric choices matter. | Community turnover analysis. |
| `turnover` | `taxa-turnover` | Planned | planned | planned | planned | [microsuite Python](containers/microsuite/Dockerfile) | Taxon-centric interpretation; sensitive to filtering. | Taxa turnover analysis. |
| `visualize` | `native` | microsuite 0.1.0 | ready | `microsuite viz ...` | `microsuite.viz.barplot.taxonomy_barplot` | [microsuite Python](containers/microsuite/Dockerfile) | Lightweight static figures; not an interactive dashboard. | Barplots, ordination plots, and heatmaps. |
| `report` | `native` | microsuite 0.1.0 | ready | `microsuite report --backend native` | `microsuite.methods.report.report` | [microsuite Python](containers/microsuite/Dockerfile) | Good provenance summary; not full narrative reporting yet. | HTML provenance reports from run metadata. |

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
