# microsuite

**One microbiome toolbox that runs the same methods anywhere — your laptop, a
container, an HPC cluster, or a cloud VM.**

`microsuite` lets you go from raw reads to feature tables, taxonomy, diversity,
and differential abundance without rebuilding your pipeline for each compute
environment. It wraps the tools you already know (FastQC, Cutadapt, DADA2,
VSEARCH, QIIME 2, ANCOM-BC, and more) behind one consistent interface, and adds
fast native methods for table handling and statistics.

## Try it in 60 seconds

```bash
uv sync --extra dev          # install for local development
uv run microsuite --help     # see what's available
uv run microsuite workflow list
```

Run a complete demo pipeline on bundled example data:

```bash
uv run microsuite workflow moving-pictures --out runs/moving-pictures --force
```

No QIIME 2, R, or external databases needed for the native demo — those tools
are only required for backends that explicitly depend on them.

## Three ways to use it

Pick whichever fits how you work:

| If you want to... | Use | Example |
| --- | --- | --- |
| Run a full multi-step pipeline reproducibly | **Nextflow workflows** | `workflows/nextflow/main.nf` |
| Run one task at a time on local files | **CLI commands** | `microsuite trim --backend fastp ...` |
| Call methods from a notebook or script | **Python SDK** | `from microsuite.api import alpha_diversity` |

Processing many samples at once? See the [multisample & concurrency guide](docs/multisample.md).

The native data object is **AnnData**, but QIIME 2 artifacts (`.qza`), FASTQ
files, and external-tool formats are supported at the workflow boundaries.

## What can it do?

`microsuite` covers the common amplicon and metagenomics steps. Each task can
run through several backends, so you can swap tools without changing your
workflow shape:

| Task | CLI command | Example backends |
| --- | --- | --- |
| Quality reports | `qc` | FastQC, MultiQC, QIIME 2 demux |
| Quality filtering | `qc_filter`, `decontam` | QIIME 2 quality-control, decontam |
| Read trimming | `trim` | fastp, Cutadapt, Trimmomatic, Trim Galore |
| Primer preflight | `primer-check` | Native IUPAC-aware FASTQ inspection |
| Denoise to ASVs | `denoise` | QIIME 2 DADA2, Deblur, R DADA2 |
| Cluster to OTUs | `cluster` | VSEARCH, USEARCH, QIIME 2 VSEARCH |
| Assemble / bin metagenomes | `assemble`, `bin` | MEGAHIT, metaSPAdes, MetaBAT2, MaxBin2, CONCOCT |
| Taxonomy & phylogeny | `tax_classify`, `phylogeny` | QIIME 2, Kraken2, Bracken, MetaPhlAn, EMU, MAFFT/FastTree |
| Table transforms & summaries | `normalize`, `abundance`, `rarefy` | native, QIIME 2 feature-table |
| Batch effect correction | `batch correct` | MMUPHin, ComBat-seq, ConQuR, PLSDA-batch, MetaDICT |
| Diversity & ecology stats | `diversity`, `diversity_calc` | native, QIIME 2 diversity-lib |
| Differential abundance | `diff_abundance` | ANCOM-BC, ALDEx2, MaAsLin2, LEfSe |
| Networks | `network infer` | native correlation, SparCC, SPIEC-EASI, FlashWeave |
| Functional profiling | `functional_profile` | PICRUSt2, Tax4Fun2, HUMAnN |
| ML & longitudinal | `ml classify`, `ml longitudinal` | random forest, XGBoost, native time-series |
| Visualization & reports | `viz`, `report` | native figures and HTML provenance reports |

👉 **Full method catalog** — every backend, version, container, and tradeoff —
is in [docs/methods.md](docs/methods.md).

## Install

Choose the path that matches how you want to use `microsuite`:

| You want to... | Use this path |
| --- | --- |
| Try the CLI or run examples from this repo | [Developer/local install](docs/installation.md#developerlocal-install) |
| Use functions in notebooks or scripts | [Python SDK install](docs/installation.md#python-sdk-users) |
| Run one command at a time on local files | [CLI users](docs/installation.md#cli-users) |
| Avoid installing FastQC, VSEARCH, R, or QIIME 2 locally | [Docker users](docs/installation.md#docker-users) |
| Run full reproducible workflows | [Nextflow/HPC users](docs/installation.md#nextflow-and-hpc-users) |
| Run BioProject-scale jobs on cheap cloud VMs | [Cloud/Spot VM users](docs/installation.md#cloudspot-vm-users) |

A core Python install gives you the `microsuite` CLI/SDK and native table
analysis only. External tools stay separate: QIIME 2 commands need QIIME 2,
R-backed methods need `Rscript` plus R packages, and taxonomy profilers need
their own databases. Full instructions, optional extras, Docker images, and
verification commands are in [docs/installation.md](docs/installation.md).

## CLI in practice

`microsuite` commands fall into three layers:

- **Method commands** (`tax_classify`, `denoise`, `trim`, ...) — run one named
  task with an explicit `--backend`.
- **Workflow commands** (`microsuite workflow ...`) — run complete pipelines.
- **Building blocks** (`import`, `diversity`, `ordination`, `viz`, `qiime`) —
  lower-level steps for custom runs.

Commands overwrite outputs only when `--force` is supplied.

A few representative commands:

```bash
# Trim paired reads with fastp
microsuite trim \
  --backend fastp \
  --read1 sample_R1.fastq.gz --read2 sample_R2.fastq.gz \
  --output1 trimmed_R1.fastq.gz --output2 trimmed_R2.fastq.gz \
  --html qc/fastp.html --json-report qc/fastp.json

# Denoise to ASVs with QIIME 2 DADA2
microsuite denoise \
  --backend qiime2-dada2 \
  --demux demux.qza \
  --output-table table.qza --output-rep-seqs rep-seqs.qza --output-stats stats.qza \
  --trunc-len 150

# Check a Cutadapt primer against a deterministic FASTQ sample
microsuite primer-check \
  --input sample_R1.fastq.gz \
  --front '^TCAGNNNNNNNNNNGGATTAGATACCCTGGTAGT' \
  --mode error -o qc/primer-check.json

# Native table analysis (no external tools needed)
microsuite import tsv table.tsv --metadata metadata.tsv --taxonomy taxonomy.tsv -o table.h5ad
microsuite diversity alpha table.h5ad --metric shannon -o alpha.tsv
microsuite diversity beta table.h5ad --metric bray-curtis -o beta.tsv
microsuite ordination pcoa beta.tsv -o pcoa.tsv
microsuite viz barplot table.h5ad --level genus -o barplot.png

# Classify taxonomy and test differential abundance
microsuite tax_classify --backend qiime2 --rep-seqs rep-seqs.qza --classifier classifier.qza -o taxonomy.qza
microsuite diff_abundance --backend ancombc --table table.h5ad --group treatment -o diff.tsv

# Generate an HTML provenance report
microsuite report --backend native --run-dir runs/table-summary -o report.html
```

Run `microsuite methods` to list every available method/backend. For the full
command set see [docs/api-cli.md](docs/api-cli.md).

The full `moving-pictures-qiime2` workflow reproduces the QIIME 2 Moving
Pictures tutorial (raw read download, import, DADA2, classifier training,
diversity/taxonomy/composition) and requires an activated QIIME 2 amplicon
environment. See
[docs/moving-pictures-qiime2-parity.md](docs/moving-pictures-qiime2-parity.md).

## Reproducibility: runtime logs

External-tool commands can write a minimal provenance bundle with `--run-dir`:

```bash
microsuite qc \
  --backend fastqc \
  --input sample_R1.fastq.gz \
  --output-dir qc/fastqc \
  --threads auto \
  --run-dir runs/fastqc/sample_R1
```

The run directory contains:

- `command.txt`: shell-quoted command line
- `stdout.log` and `stderr.log`: captured process streams
- `events.jsonl`: command start/end/timeout events
- `run.json`: structured task, backend, command, timing, and exit metadata
- `microsuite-results.json`: a stable result-bundle manifest for downstream
  consumers such as OmicScribe, listing executions plus semantic artifact
  records with `kind`, `label`, `path`, `format`, `task`, and `backend`

Commands that expose `--threads` accept a positive integer; method-oriented
external wrappers also accept `--threads auto`, which uses the detected CPU
count minus one reserved core.

## Documentation

- [docs/methods.md](docs/methods.md) — full method/backend catalog and validation status
- [docs/installation.md](docs/installation.md) — install by user type
- [docs/three-api-roadmap.md](docs/three-api-roadmap.md) — architecture
- [docs/toolbox.md](docs/toolbox.md) — how the toolbox is organized
- [docs/data-attribution.md](docs/data-attribution.md) — demo data and citations
- [docs/dada2.md](docs/dada2.md) — DADA2 backend guidance
- [docs/external-integration-tests.md](docs/external-integration-tests.md) — optional real-tool tests

## Development

```bash
uv run ruff check .
uv run ruff format --check .
uv run ty check
uv run pytest
uv build
```

GitHub Actions runs the same quality gate on Python 3.11 and 3.12.
