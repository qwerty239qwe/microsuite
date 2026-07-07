# Method reference

Full backend catalog for `microsuite`. Every supported method, its CLI command,
Python invocation, container/environment, and operational tradeoff.

For a friendly overview and quick start, see the [README](../README.md).

## Backend Validation Status

This reference separates biological task support from runtime validation. The
`Status` column in the method tables describes API maturity: `ready`, `partial`,
or `planned`. Runtime validation is tracked separately:

| Level | Meaning |
| --- | --- |
| CI smoke-tested | The command or container is exercised in GitHub Actions with a lightweight smoke test. |
| Unit-tested wrapper | Command construction, validation, and error handling are covered by Python tests, but the external tool is not executed in CI. |
| Static only | Files, docs, or container skeletons are checked, but the backend is not runnable as part of the default test suite. |
| User environment | The backend requires tools, plugins, databases, R packages, or QIIME 2 environments supplied by the user. |
| Planned | Listed to reserve API shape, but not implemented for 0.1.0. |

| Backend family | API status | Validation level | Notes |
| --- | --- | --- | --- |
| Native table/statistics/report methods | ready | CI smoke-tested | Covered by unit and CLI workflow tests. |
| FastQC | ready | CI smoke-tested | CLI wrapper and container are smoke-tested. |
| MultiQC, fastp, Cutadapt, Trimmomatic, Trim Galore | ready | CI smoke-tested + unit-tested wrapper | fastp and MultiQC containers are smoke-tested; other wrappers are command-tested and user-supplied. |
| QIIME 2 method wrappers | ready | Unit-tested wrapper + user environment | Command construction is tested; QIIME 2/plugin version validation is user supplied. |
| QIIME 2 MOSHPIT | partial | Unit-tested wrapper + user environment | Initial `mosh` command wrappers cover MOSHPIT MEGAHIT assembly and MetaBAT2 contig binning; broader MOSHPIT actions remain planned. See https://moshpit.qiime2.org/en/stable/. |
| R differential-abundance methods | ready | Unit-tested wrapper + user environment | Python wrappers and runtime logs are tested; R/Bioconductor runtime is user supplied or containerized. |
| Nextflow workflows | ready | CI stub-tested + user environment | The process graph is exercised with Nextflow `-stub-run`; real QIIME 2 execution requires local tools or containers. |
| Kraken2, Bracken, MetaPhlAn, EMU, ALDEx2, MaAsLin2, LEfSe | ready | Kraken2 CI smoke-tested; MetaPhlAn heavy image manual-gated; EMU/R wrappers unit-tested | External runtimes and databases remain user supplied. |
| Metagenome assembly and binning wrappers | ready | Unit-tested wrapper + user environment | MEGAHIT, metaSPAdes, IDBA-UD, MetaBAT2, MaxBin2, and CONCOCT command construction is tested; real execution requires user-supplied tools and input matrices. |

## Method Surface

### Quality Reports

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `fastqc` | FastQC 0.12.1 | ready | `microsuite qc --backend fastqc` | `qc(backend="fastqc", inputs=[...], output_dir=...)` | [FastQC](../containers/fastqc/Dockerfile) or external `fastqc` | Ready as a CLI wrapper and standalone container; Nextflow raw-read wiring remains planned. | Raw-read quality reports. |
| `multiqc` | MultiQC user env | ready | `microsuite qc --backend multiqc` | `qc(backend="multiqc", input_dir=..., output_dir=...)` | [MultiQC](../containers/multiqc/Dockerfile) or external `multiqc` | Good aggregation layer; depends on upstream report files. | Aggregate QC reports. |
| `qiime2-demux-summarize` | QIIME 2 2024.10 | ready | `microsuite qc --backend qiime2-demux` | `qc(backend="qiime2-demux", demux=..., output=...)` | [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Best for QIIME artifacts; wrapper requires an activated QIIME 2 environment or container. | Demultiplexed-read quality visualization. |

### Quality Filtering

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `qiime2-exclude-seqs` | QIIME 2 user env | ready | `microsuite qc_filter --backend qiime2-exclude-seqs` | `qc_filter(backend="qiime2-exclude-seqs", query_sequences=..., reference_sequences=...)` | External QIIME 2 with `q2-quality-control`; [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Strong contaminant/non-target sequence filtering; needs reference sequences and threshold choices. | Exclude or retain feature sequences by alignment to reference sequences. |
| `qiime2-filter-reads` | QIIME 2 user env | ready | `microsuite qc_filter --backend qiime2-filter-reads` | `qc_filter(backend="qiime2-filter-reads", demux=..., database=..., output=...)` | External QIIME 2 with `q2-quality-control`; [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Useful for host/contaminant read removal; requires a Bowtie2 index. | Filter demultiplexed reads by alignment to a reference database. |
| `qiime2-bowtie2-build` | QIIME 2 user env | ready | `microsuite qc_filter --backend qiime2-bowtie2-build` | `qc_filter(backend="qiime2-bowtie2-build", sequences=..., output=...)` | External QIIME 2 with `q2-quality-control`; [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Completes the filter-reads setup path; still requires suitable reference sequences. | Build a Bowtie2 index artifact for read filtering. |
| `qiime2-decontam` | QIIME 2 user env | ready | `microsuite decontam --backend qiime2-decontam` | `decontam(backend="qiime2-decontam", table=..., metadata=..., output=...)` | External QIIME 2 with `q2-quality-control`; [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Useful contamination screening; requires negative controls or concentration metadata. | Identify likely contaminant features with decontam. |
| `qiime2-quality-filter-q-score` | QIIME 2 2024.10 | ready | `microsuite qc_filter --backend qiime2-quality-filter-q-score` | `qc_filter(backend="qiime2-quality-filter-q-score", demux=..., output=..., sequence_hits=...)` | [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Native QIIME quality-score filtering; mainly useful before downstream QIIME artifact workflows. | Filter demultiplexed reads by quality scores. |

### Trimming

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `fastp` | fastp user env | ready | `microsuite trim --backend fastp` | `trim(backend="fastp", read1=..., output1=...)` | [fastp](../containers/fastp/Dockerfile) or external `fastp` | Fast all-in-one preprocessing; supports shared quality, length, N, and adapter options; primer-specific trimming is less explicit than Cutadapt. | Adapter trimming, quality filtering, HTML/JSON reports. |
| `cutadapt` | Cutadapt >=4.x user env | ready | `microsuite trim --backend cutadapt` | `trim(backend="cutadapt", read1=..., output1=..., adapter=...)` | [Cutadapt](../containers/cutadapt/Dockerfile) or external `cutadapt` | Precise primer/adaptor trimming with explicit adapter control; requires users to choose primer/adapter sequences. | Adapter/primer trimming and read filtering. |
| `trimmomatic` | Trimmomatic >=0.39 user env | ready | `microsuite trim --backend trimmomatic` | `trim(backend="trimmomatic", read1=..., output1=..., trimmomatic_steps=[...])` | [Trimmomatic](../containers/trimmomatic/Dockerfile) or external `trimmomatic` | Mature Java trimmer with explicit step pipeline; paired mode requires unpaired output files. | Sliding-window, length, quality, and adapter trimming. |
| `trim-galore` | Trim Galore 0.6.x or v2.x user env | ready | `microsuite trim --backend trim-galore` | `trim(backend="trim-galore", read1=..., output1=..., trim_galore_version="auto")` | [Trim Galore](../containers/trim-galore/Dockerfile) or external `trim_galore` | Lets users keep tool-default behavior or explicitly select the v2 mode; output names are tool-controlled and validated. | Adapter/quality trimming with integrated QC conventions. |
| `qiime2-cutadapt` | QIIME 2 2024.10 | ready | `microsuite trim --backend qiime2-cutadapt --read1 demux.qza --output1 trimmed.qza` | `trim(backend="qiime2-cutadapt", read1=demux, output1=trimmed, adapter=...)` | [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Fits QIIME artifact workflows; less convenient for raw FASTQ-only runs. | QIIME 2 Cutadapt wrapper. |

> For batching many samples and choosing threads vs parallel jobs, see [multisample runs and concurrency](multisample.md).

### Feature Table Generation

ASV tables contain inferred exact sequence variants. OTU-style tables contain
features clustered by sequence identity, commonly 97%.

| Goal | Backend | Input | Output table type | Main outputs | Notes |
| --- | --- | --- | --- | --- | --- |
| Generate ASV table | `qiime2-dada2` | Demultiplexed QIIME 2 reads artifact | ASV feature table | `table.qza`, `rep-seqs.qza`, denoising stats, optional base-transition diagnostics | Main QIIME 2 ASV path, including single, paired, CCS, and pyro modes. |
| Generate ASV table | `qiime2-deblur` | Demultiplexed QIIME 2 reads artifact | ASV feature table | `table.qza`, `rep-seqs.qza`, Deblur stats | 16S-oriented Deblur path. |
| Generate ASV table | `dada2-r` | FASTQ directory | ASV count table | table TSV with `ASV*` feature IDs, representative-sequence FASTA, stats TSV | Direct R/DADA2 path; not a QIIME artifact workflow. |
| Generate OTU-style table | `vsearch` | FASTA sequences with sample IDs in labels | OTU count table | table TSV, centroid FASTA, sidecar `.uc` mapping | Standalone VSEARCH clustering. |
| Generate OTU-style table | `usearch` | FASTA sequences with sample IDs in labels | OTU count table | table TSV, centroid FASTA, sidecar `.uc` mapping | Standalone USEARCH clustering. |
| Generate QIIME-clustered table | `qiime2-vsearch` | QIIME 2 feature table and representative sequences | Clustered QIIME 2 feature table | clustered table `.qza`, clustered representative sequences `.qza` | QIIME VSEARCH de novo clustering. |

### Denoising And Clustering Backends

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `qiime2-dada2` | q2-dada2 2026.4.0 target | ready | `microsuite denoise --backend qiime2-dada2 --mode single` | `denoise(backend="qiime2-dada2", demux=..., output_table=..., mode="single")` | [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) or external QIIME 2 | Best for QIIME artifact workflows and q2-dada2 modes such as CCS and pyro; needs careful truncation choices. | DADA2 ASV inference from demultiplexed reads. |
| `qiime2-deblur` | QIIME 2 2024.10 | ready | `microsuite denoise --backend qiime2-deblur` | `denoise(backend="qiime2-deblur", demux=..., output_table=...)` | [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Reproducible fixed-error model; mainly 16S-oriented. | Deblur ASV inference from demultiplexed reads. |
| `dada2-r` | DADA2 1.40.0 / Bioconductor 3.23 target | ready | `microsuite denoise --backend dada2-r` | `denoise(backend="dada2-r", demux=reads_dir, output_table=...)` | [R DADA2](../containers/r-dada2/Dockerfile) or external `Rscript` with R package `dada2` | Best for direct FASTQ-to-TSV/FASTA workflows; expects a FASTQ directory and writes importable outputs with matching ASV IDs. | R/DADA2 ASV inference from raw or trimmed FASTQ files. |
| `vsearch` | VSEARCH user env | ready | `microsuite cluster --backend vsearch` | `cluster(backend="vsearch", rep_seqs=..., output_table=..., output_rep_seqs=...)` | [VSEARCH](../containers/vsearch/Dockerfile) or external `vsearch` | Standalone OTU-style clustering with TSV count table output; sample IDs are inferred from sequence labels. | VSEARCH `cluster_fast` sequence clustering. |
| `usearch` | USEARCH 12 | ready | `microsuite cluster --backend usearch` | `cluster(backend="usearch", rep_seqs=..., output_table=..., output_rep_seqs=...)` | [USEARCH 12](../containers/usearch/Dockerfile) or external `usearch` | Fast standalone OTU-style clustering with TSV count table output; sample IDs are inferred from sequence labels. | USEARCH `cluster_fast` sequence clustering. |
| `qiime2-vsearch` | QIIME 2 2024.10 | ready | `microsuite cluster --backend qiime2-vsearch` | `cluster(backend="qiime2-vsearch", table=..., rep_seqs=...)` | [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Useful for QIIME artifact workflows; requires QIIME 2 and q2-vsearch. | QIIME 2 VSEARCH de novo feature clustering. |

### Metagenome Assembly And Binning

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `megahit` | MEGAHIT user env | ready | `microsuite assemble --backend megahit --read1 R1.fq.gz --read2 R2.fq.gz --output-dir assembly` | `assemble(backend="megahit", read1=..., read2=..., output_dir=...)` | External `megahit` | Fast metagenome assembler with a compact command surface; paired or single reads are supported. | Assemble metagenomic reads into contigs. |
| `mosh-megahit` | MOSHPIT 2026.4 target | ready | `microsuite assemble --backend mosh-megahit --reads reads.qza --output-contigs contigs.qza` | `assemble(backend="mosh-megahit", reads=..., output_contigs=...)` | External MOSHPIT `mosh` CLI | QIIME 2 provenance-preserving MEGAHIT assembly over demultiplexed metagenome read artifacts. | Assemble MOSHPIT/QIIME metagenome reads into contigs artifacts. |
| `metaspades` | SPAdes/metaSPAdes user env | ready | `microsuite assemble --backend metaspades --read1 R1.fq.gz --read2 R2.fq.gz --output-dir assembly` | `assemble(backend="metaspades", read1=..., read2=..., output_dir=...)` | External `metaspades.py` | Strong metagenome assembly option; can be heavier than MEGAHIT. | Assemble metagenomic reads into contigs. |
| `idba-ud` | IDBA-UD user env | ready | `microsuite assemble --backend idba-ud --reads reads.fa --output-dir assembly` | `assemble(backend="idba-ud", reads=..., output_dir=...)` | External `idba_ud` | Expects FASTA input, so FASTQ conversion remains a preprocessing step. | Assemble metagenomic reads into contigs. |
| `metabat2` | MetaBAT2 user env | ready | `microsuite bin --backend metabat2 --contigs contigs.fa --depth depth.tsv --output-dir bins` | `bin_contigs(backend="metabat2", contigs=..., depth=..., output_dir=...)` | External `metabat2` | Requires a depth matrix, usually generated from read mapping. | Bin contigs into MAG candidates. |
| `mosh-metabat2` | MOSHPIT 2026.4 target | ready | `microsuite bin --backend mosh-metabat2 --contigs contigs.qza --alignment-maps reads-to-contigs-aln.qza --output-dir bins` | `bin_contigs(backend="mosh-metabat2", contigs=..., alignment_maps=..., output_dir=...)` | External MOSHPIT `mosh` CLI | QIIME 2 provenance-preserving MetaBAT2 binning using MOSHPIT contigs and read-alignment artifacts. | Bin MOSHPIT contigs into MAG artifacts. |
| `maxbin2` | MaxBin2 user env | ready | `microsuite bin --backend maxbin2 --contigs contigs.fa --abundance abundance.tsv --output-dir bins` | `bin_contigs(backend="maxbin2", contigs=..., abundance=..., output_dir=...)` | External `run_MaxBin.pl` | Requires an abundance table compatible with MaxBin2. | Bin contigs into MAG candidates. |
| `concoct` | CONCOCT user env | ready | `microsuite bin --backend concoct --contigs contigs.fa --coverage coverage.tsv --output-dir bins` | `bin_contigs(backend="concoct", contigs=..., coverage=..., output_dir=...)` | External `concoct` | Requires a coverage table; downstream bin FASTA extraction is not automated yet. | Bin contigs into MAG candidates. |

### Taxonomy And Phylogeny

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `qiime2` | QIIME 2 2024.10 | ready | `microsuite tax_classify --backend qiime2` | `tax_classify(backend="qiime2", rep_seqs=..., classifier=..., output=...)` | [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Strong classifier ecosystem; requires trained classifier artifacts. | QIIME 2 taxonomy classification. |
| `qiime2-taxonomy` | QIIME 2 user env | ready | `microsuite evaluate --backend qiime2-taxonomy` | `evaluate(backend="qiime2-taxonomy", expected_taxa=..., observed_taxa=...)` | External QIIME 2 with `q2-quality-control`; [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Best for mock-community or known-composition validation; needs expected taxonomy. | Evaluate observed taxonomy against expected assignments. |
| `kraken2` | Kraken2 2.1.3 | ready | `microsuite tax_classify --backend kraken2 --classifier DB` | `tax_classify(backend="kraken2", rep_seqs=..., classifier=database, output=...)` | [Kraken2](../containers/kraken2/Dockerfile) or external `kraken2` | Fast profiling; requires user-supplied databases. `--output` is the report and `.kraken` sidecar stores per-read assignments. | Taxonomic profiling/classification. |
| `bracken` | Bracken user env | ready | `microsuite tax_classify --backend bracken --classifier DB --level S --read-length 150` | `tax_classify(backend="bracken", rep_seqs=kraken_report, classifier=database, output=..., level="S", read_length=150)` | [Kraken2/Bracken](../containers/kraken2/Dockerfile) or external `bracken` | Re-estimates abundance from a Kraken report; requires a Bracken-prepared Kraken database and matched read length. | Abundance re-estimation from Kraken2 output. |
| `metaphlan` | MetaPhlAn user env | ready | `microsuite tax_classify --backend metaphlan --input-type fastq` | `tax_classify(backend="metaphlan", rep_seqs=..., output=..., input_type="fastq")` | [MetaPhlAn](../containers/metaphlan/Dockerfile) or external `metaphlan` | Marker-gene profiling; can use MetaPhlAn's configured default database or `--classifier` as a Bowtie2 database directory. Writes a `.bowtie2.bz2` sidecar. | Marker-gene taxonomic profiling. |
| `emu` | EMU user env | ready | `microsuite tax_classify --backend emu --input-type map-ont --classifier EMU_DB -o sample_rel-abundance.tsv` | `tax_classify(backend="emu", rep_seqs=long_reads, classifier=emu_db, output=..., input_type="map-ont")` | External `emu` with an EMU database | Species-level profiling for full-length long-read 16S; output path must end `_rel-abundance.tsv` to match EMU naming. | Long-read amplicon taxonomic profiling. |
| `qiime2-phylogeny` | QIIME 2 2024.10 | ready | `microsuite phylogeny --backend qiime2-mafft-fasttree` | `phylogeny(backend="qiime2-mafft-fasttree", rep_seqs=...)` | [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Integrated with QIIME artifacts; heavier runtime. | Alignment, masking, tree construction, and rooting. |
| `mafft-fasttree` | MAFFT + FastTree user env | ready | `microsuite phylogeny --backend mafft-fasttree` | `phylogeny(backend="mafft-fasttree", rep_seqs=..., output_aligned=..., output_tree=...)` | [MAFFT/FastTree](../containers/mafft-fasttree/Dockerfile) or external `mafft` and `FastTree` | Lightweight FASTA-to-Newick path; optional masked/rooted outputs are compatibility copies, not QIIME-style masking/rooting. | Standalone MAFFT/FastTree phylogeny. |

### Table Transforms And Summaries

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `native-normalize` | microsuite 0.1.0 | ready | `microsuite normalize --backend native` | `normalize_table(adata, method="relative")` | [microsuite Python](../containers/microsuite/Dockerfile) | Fast and portable; narrower than specialized compositional packages. | Relative abundance, CLR, and table transforms. |
| `native-abundance` | microsuite 0.1.0 | ready | `microsuite abundance --backend native` | `abundance_table(adata, level="genus")` | [microsuite Python](../containers/microsuite/Dockerfile) | Simple summary output; depends on taxonomy quality. | Summarize abundance at taxonomy levels. |
| `native-shared-taxa` | microsuite 0.1.0 | ready | `microsuite shared_taxa --backend native` | `shared_taxa_table(adata, level="genus", group=...)` | [microsuite Python](../containers/microsuite/Dockerfile) | Easy group comparison; descriptive rather than inferential. | Compare shared taxa across sample groups. |
| `native-rarefy` | microsuite 0.1.0 | ready | `microsuite rarefy --backend native` | `rarefy_table(adata, depth=...)` | [microsuite Python](../containers/microsuite/Dockerfile) | Reproducible with seeds; discards reads by design. | Rarefy feature tables to a fixed depth. |
| `qiime2-feature-table` | QIIME 2 2024.10 | ready | `microsuite feature_summarize --backend qiime2 --mode summarize` | `feature_summarize(backend="qiime2", mode="summarize", table=...)` | [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Keeps QIIME-native summaries and sequence tabulation in artifact form. | Feature-table summary and representative-sequence visualization. |
| `qiime2-taxa` | QIIME 2 2024.10 | ready | `microsuite tax_barplot --backend qiime2`; `microsuite tax_collapse --backend qiime2` | `tax_barplot(...)`; `tax_collapse(...)` | [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Artifact-native taxonomy visualization and collapse. | Taxa barplots and taxonomy-level table collapse. |

### Diversity And Ecological Statistics

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `native` | microsuite 0.1.0 | ready | lower-level `microsuite diversity ...` | `alpha_diversity(adata, metric=...)` or `beta_diversity(adata, metric=...)` | [microsuite Python](../containers/microsuite/Dockerfile) | Lightweight and Windows-friendly; phylogenetic metrics are limited. | Native alpha/beta diversity. |
| `qiime2-diversity-lib` | QIIME 2 2024.10 | ready | `microsuite diversity_calc --backend qiime2` | `diversity_calc(backend="qiime2", metric=..., table=..., output=...)` | [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Broader metric coverage; requires QIIME artifacts. | QIIME 2 diversity-lib metrics. |
| `qiime2-core-metrics-phylogenetic` | QIIME 2 2024.10 | ready | `microsuite diversity_core --backend qiime2-core-metrics-phylogenetic` | `diversity_core(...)` | [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Matches the Moving Pictures core diversity path. | Core phylogenetic diversity metrics. |
| `qiime2-alpha/beta-significance` | QIIME 2 2024.10 | ready | `microsuite diversity_test --backend qiime2-beta-group-significance` | `diversity_test(...)` | [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Established group-significance workflow; artifact-based. | Alpha group tests and PERMANOVA beta-diversity tests. |
| `native-beta-significance` | microsuite 0.1.0 | ready | `microsuite diversity beta-significance beta.tsv --metadata metadata.tsv --column body_site` | `beta_significance(distance_matrix, metadata, column=..., method="permanova")` | [microsuite Python](../containers/microsuite/Dockerfile) | Permutation PERMANOVA/ANOSIM without QIIME artifacts; users must choose valid grouping designs. | Native beta-diversity tests. |
| `mantel` | microsuite 0.1.0 | ready | `microsuite diversity mantel beta-a.tsv beta-b.tsv -o mantel.tsv` | `mantel_test(matrix_a, matrix_b, method="spearman")` | [microsuite Python](../containers/microsuite/Dockerfile) | Useful distance association; sensitive to spatial and repeated-measures design assumptions. | Mantel association testing. |
| `rda` | microsuite 0.1.0 | ready | `microsuite diversity constrained-ordination table.h5ad --constraint body_site --method rda -o rda.tsv` | `constrained_ordination(adata, constraints=[...], method="rda")` | [microsuite Python](../containers/microsuite/Dockerfile) | Dependency-light constrained ordination; advanced permutation designs remain external-tool territory. | Redundancy analysis. |
| `cca` | microsuite 0.1.0 | ready | `microsuite diversity constrained-ordination table.h5ad --constraint body_site --method cca -o cca.tsv` | `constrained_ordination(adata, constraints=[...], method="cca")` | [microsuite Python](../containers/microsuite/Dockerfile) | Chi-square standardized approximation; users still need ecological judgment. | Canonical correspondence analysis. |
| `db-rda` | microsuite 0.1.0 | ready | `microsuite diversity constrained-ordination table.h5ad --constraint body_site --method db-rda -o dbrda.tsv` | `constrained_ordination(adata, constraints=[...], method="db-rda")` | [microsuite Python](../containers/microsuite/Dockerfile) | Uses Bray-Curtis PCoA coordinates before constrained fitting; permutation design remains explicit user responsibility. | Distance-based redundancy analysis. |
| `native-gamma-diversity` | microsuite 0.1.0 | ready | `microsuite diversity gamma table.h5ad --group body_site -o gamma.tsv` | `gamma_diversity(adata, group=..., metric="observed_features")` | [microsuite Python](../containers/microsuite/Dockerfile) | Pooled group-level alpha metric; definition is explicit in output columns. | Region/group-level diversity summaries. |
| `beta-turnover` | microsuite 0.1.0 | ready | `microsuite diversity beta-turnover table.h5ad --level genus -o turnover.tsv` | `beta_turnover(adata, level="genus")` | [microsuite Python](../containers/microsuite/Dockerfile) | Presence/absence Sorensen decomposition; abundance-weighted turnover is separate. | Community turnover analysis. |
| `taxa-turnover` | microsuite 0.1.0 | ready | `microsuite diversity taxa-turnover table.h5ad --group body_site --level genus -o taxa-turnover.tsv` | `taxa_turnover(adata, group=..., level="genus")` | [microsuite Python](../containers/microsuite/Dockerfile) | Taxon-centric group comparison; sensitive to feature filtering and taxonomy level. | Taxa turnover analysis. |

### Differential Abundance

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `ancombc` | R image; ANCOMBC via Bioconductor | ready | `microsuite diff_abundance --backend ancombc` | `diff_abundance(backend="ancombc", table=..., group=..., output=...)` | [R diffab](../containers/r-diffab/Dockerfile) | Strong compositional method; R/Bioconductor image is heavier. | ANCOM-BC differential abundance. |
| `qiime2-ancombc` | QIIME 2 2024.10 | ready | `microsuite diff_abundance --backend qiime2-ancombc --metadata sample-metadata.tsv` | `diff_abundance(backend="qiime2-ancombc", table=..., metadata=..., group=...)` | [QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile) | Keeps composition results in QIIME artifact form. | QIIME composition ANCOM-BC. |
| `aldex2` | R image; ALDEx2 via Bioconductor | ready | `microsuite diff_abundance --backend aldex2` | `diff_abundance(backend="aldex2", table=..., group=..., output=...)` | [R diffab](../containers/r-diffab/Dockerfile) | Good compositional alternative; uses two-group t/Wilcoxon plus effect size or Kruskal-Wallis for multi-group designs. | ALDEx2 differential abundance. |
| `maaslin2` | R image; MaAsLin2 via Bioconductor | ready | `microsuite diff_abundance --backend maaslin2` | `diff_abundance(backend="maaslin2", table=..., group=..., output=...)` | [R diffab](../containers/r-diffab/Dockerfile) | Flexible covariate modeling; current wrapper exposes the selected group as a single fixed effect. | MaAsLin2 multivariable association testing. |
| `lefse` | R image; lefser via Bioconductor | ready | `microsuite diff_abundance --backend lefse` | `diff_abundance(backend="lefse", table=..., group=..., output=...)` | [R diffab](../containers/r-diffab/Dockerfile) | Familiar legacy workflow; current wrapper supports two-class comparisons through lefser. | LEfSe legacy differential abundance. |

### Networks

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `native-correlation` | microsuite 0.1.0 | ready | `microsuite network infer --backend native-correlation --table table.h5ad -o network.tsv` | `network(backend="native-correlation", table=..., output=...)` | [microsuite Python](../containers/microsuite/Dockerfile) | Simple and transparent; compositional bias risk. | Correlation network analysis. |
| `sparcc` | microsuite 0.1.0 | ready | `microsuite network infer --backend sparcc --table table.h5ad -o sparcc.tsv` | `network(backend="sparcc", table=..., output=...)` | [microsuite Python](../containers/microsuite/Dockerfile) | CLR correlation approximation for compositional data; full SparCC bootstrapping remains external. | SparCC-style association network analysis. |
| `spieceasi` | SpiecEasi R user env | ready | `microsuite network infer --backend spieceasi --table table.h5ad -o spieceasi.tsv` | `network(backend="spieceasi", table=..., output=...)` | External R/SpiecEasi environment | Strong ecological network method; R dependency and tuning burden. | SPIEC-EASI network inference. |
| `flashweave` | FlashWeave Julia user env | ready | `microsuite network infer --backend flashweave --table table.h5ad -o flashweave.edgelist` | `network(backend="flashweave", table=..., output=...)` | External Julia/FlashWeave environment | Handles heterogeneous metadata; separate Julia runtime needed. | FlashWeave network inference. |

### Functional Profiling

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `picrust2` | PICRUSt2 user env | ready | `microsuite functional_profile --backend picrust2 --table table.biom --rep-seqs rep-seqs.fasta --output-dir functions` | `functional_profile(backend="picrust2", table=..., rep_seqs=..., output_dir=...)` | External PICRUSt2 environment | Popular marker-gene function prediction; reference-dependent and expects PICRUSt2 input formats. | Predict function from marker-gene profiles. |
| `tax4fun2` | Tax4Fun2 R user env | ready | `microsuite functional_profile --backend tax4fun2 --table otu-table.tsv --rep-seqs otus.fasta --database Tax4Fun2_ReferenceData_v2 --output-dir functions` | `functional_profile(backend="tax4fun2", table=..., rep_seqs=..., database=..., output_dir=...)` | External R/Tax4Fun2 environment | Alternative functional prediction; requires Tax4Fun2 reference data and BLAST dependencies. | Tax4Fun2 function prediction. |
| `humann` | HUMAnN user env | ready | `microsuite functional_profile --backend humann --reads reads.fastq.gz --output-dir functions` | `functional_profile(backend="humann", reads=..., output_dir=...)` | External HUMAnN environment | Strong metagenomic functional profiling; heavy nucleotide/protein databases. | Functional profiling from metagenomic data. |

### Machine Learning And Longitudinal Analysis

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `randomforest` | microsuite 0.1.0 | ready | `microsuite ml classify --backend randomforest --table table.h5ad --target body_site -o predictions.tsv` | `ml_classify(backend="randomforest", table=..., target=..., output=...)` | [microsuite Python](../containers/microsuite/Dockerfile) | Interpretable baseline ML; uses scikit-learn when available and a deterministic fallback otherwise. | Supervised sample classification. |
| `xgboost` | XGBoost optional Python package | ready | `microsuite ml classify --backend xgboost --table table.h5ad --target body_site -o predictions.tsv` | `ml_classify(backend="xgboost", table=..., target=..., output=...)` | Optional external `xgboost` Python package | Strong predictive model; optional dependency and tuning burden. | Optional XGBoost sample classification. |
| `native-time-series` | microsuite 0.1.0 | ready | `microsuite ml longitudinal --backend native-time-series --table table.h5ad --subject subject --time day -o slopes.tsv` | `longitudinal(backend="native-time-series", table=..., subject=..., time=..., output=...)` | [microsuite Python](../containers/microsuite/Dockerfile) | Web-friendly per-feature slope summaries; mixed-effects modeling remains external-tool territory. | Longitudinal microbiome analysis. |

### Visualization And Reporting

| Backend | Version | Status | CLI command | Python invocation | Image / environment | Operational tradeoff | Purpose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `native-visualize` | microsuite 0.1.0 | ready | `microsuite viz ...` | `taxonomy_barplot(adata, level=..., output=...)` | [microsuite Python](../containers/microsuite/Dockerfile) | Lightweight static figures; not an interactive dashboard. | Barplots, ordination plots, and heatmaps. |
| `native-report` | microsuite 0.1.0 | ready | `microsuite report --backend native` | `report(backend="native", run_dir=..., output=...)` | [microsuite Python](../containers/microsuite/Dockerfile) | Good provenance summary; not full narrative reporting yet. | HTML provenance reports from run metadata. |
