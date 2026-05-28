# Nextflow API

The Nextflow API is for full reproducible workflows.

Use it for large or publication-oriented pipelines. The current Nextflow entry
point runs the `amplicon_qiime2` workflow:

```bash
nextflow run workflows/nextflow/main.nf \
  -profile docker \
  --workflow amplicon_qiime2 \
  --manifest manifest.tsv \
  --metadata metadata.tsv \
  --classifier classifier.qza \
  --outdir results
```

Profile examples:

```bash
nextflow run workflows/nextflow/main.nf -profile local \
  --workflow amplicon_qiime2 \
  --manifest manifest.tsv \
  --metadata metadata.tsv \
  --classifier classifier.qza \
  --outdir results

nextflow run workflows/nextflow/main.nf -profile docker \
  --workflow amplicon_qiime2 \
  --manifest manifest.tsv \
  --metadata metadata.tsv \
  --classifier classifier.qza \
  --outdir results

nextflow run workflows/nextflow/main.nf -profile singularity \
  --workflow amplicon_qiime2 \
  --manifest manifest.tsv \
  --metadata metadata.tsv \
  --classifier classifier.qza \
  --outdir results
```

The Nextflow layer owns:

- multi-step orchestration
- resume/caching
- local, Docker, Singularity, HPC, and cloud profiles
- sample batching
- external tool environments

The first workflow is `amplicon_qiime2`:

```text
manifest -> FastQC -> MultiQC
manifest + reads -> QIIME 2 import -> DADA2 -> taxonomy -> phylogeny -> diversity -> report
```

## Manifest contract

Raw-read workflows should use a tab-separated manifest with one row per sample.
The required columns are:

| Column | Required | Description |
| --- | --- | --- |
| `sample_id` | yes | Stable sample identifier. Must match metadata sample IDs exactly. |
| `read1` | yes | Path to the single-end read file or forward read file. |
| `read2` | paired-end only | Path to the reverse read file. Leave absent or empty for single-end data. |
| `platform` | no | Sequencing platform label, such as `illumina`. |
| `layout` | no | `single-end` or `paired-end`; inferred from `read2` when absent. |

Rules:

- paths may be absolute or relative to the manifest file
- `sample_id` values must be unique
- single-end rows must not set `read2`
- paired-end rows must set both `read1` and `read2`
- metadata files must use the same sample IDs as the manifest
- workflow parameters should state whether reads are raw, trimmed, or already
  demultiplexed

Example:

```tsv
sample_id	read1	read2	layout
S1	reads/S1_R1.fastq.gz	reads/S1_R2.fastq.gz	paired-end
S2	reads/S2_R1.fastq.gz	reads/S2_R2.fastq.gz	paired-end
```

Runtime parameters:

| Parameter | Default | Description |
| --- | --- | --- |
| `--threads` | `2` | CPUs used by QIIME 2 DADA2 and taxonomy classification. |
| `--trim_left`, `--trunc_len` | `0`, `0` | Single-end DADA2 trimming and truncation. |
| `--trim_left_f`, `--trunc_len_f` | `0`, `0` | Paired-end forward-read DADA2 trimming and truncation. |
| `--trim_left_r`, `--trunc_len_r` | `0`, `0` | Paired-end reverse-read DADA2 trimming and truncation. |
| `--sampling_depth` | `1000` | Sampling depth for `qiime diversity core-metrics-phylogenetic`. |

The local profile expects `fastqc`, `multiqc`, and `qiime` on `PATH`. The Docker
profile assigns process-specific images for FastQC, MultiQC, QIIME 2, and the
microsuite report step. The Singularity profile expects matching `.sif` files
under `containers/singularity/`.

Continuous integration runs the workflow with Nextflow `-stub-run` against tiny
FASTQ fixtures. This validates the executable process graph without downloading
QIIME 2 databases or running heavy external tools.

Native statistics and AnnData operations should stay in the Python SDK and CLI
backends, not in Nextflow process scripts.

Current status:

- `workflows/nextflow/main.nf` orchestrates the `amplicon_qiime2` process graph.
- local, Docker, and Singularity profiles exist.
- module files contain runnable commands plus Nextflow stubs for CI smoke tests.
- default Python tests validate the workflow files statically; GitHub Actions
  runs a Nextflow stub smoke test.
