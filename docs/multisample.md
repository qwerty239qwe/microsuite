# Multisample runs and concurrency

microsuite has three ways to process many samples. They share one **sample
manifest** format and one **concurrency model**; this guide explains both and
helps you pick a path.

## The sample manifest

Batched runs are driven by a tab-separated manifest with these columns:

| Column | Meaning |
| --- | --- |
| `sample_id` | Unique sample name (no whitespace). Names the outputs. |
| `layout` | `PE` (paired-end) or `SE` (single-end). |
| `read1` | Forward (or single) FASTQ path. |
| `read2` | Reverse FASTQ path; empty for single-end. |

Example:

```tsv
sample_id	layout	read1	read2
sampleA	PE	/data/acc/sampleA_R1.fastq.gz	/data/acc/sampleA_R2.fastq.gz
sampleB	SE	/data/acc/sampleB.fastq.gz	
```

`scripts/run_fastp_multiqc.sh` writes this file — run it with `--manifest-only`
to generate the manifest and stop. The Nextflow `amplicon_qiime2` workflow reads
the same `sample_id`/`read1`/`read2` columns (it ignores the extra `layout`
column and treats an empty `read2` as single-end), so a manifest built for the
script also works with the workflow.

## Concurrency: two independent axes

**Within a sample — tool threads.** `microsuite trim`/`qc` pass `--threads` to
the underlying tool (`fastp --thread`, `cutadapt -j`, ...). `--threads auto`
uses the detected CPU count minus one reserved core. fastp scales sublinearly
past roughly 4–8 threads, so piling every core onto a single sample wastes most
of them.

**Across samples — parallel jobs.** A bare `microsuite trim`/`qc` call processes
**one** sample; the CLI itself has no cross-sample parallelism. Two paths add
it:

- `scripts/run_fastp_multiqc.sh --jobs N` runs up to `N` samples at once via
  `xargs -P`.
- The Nextflow workflow runs one process per sample and the executor schedules
  them concurrently (local = host CPUs; HPC/cloud = the cluster).

**Rule of thumb:** total cores used ≈ `jobs × threads`. For many small samples,
prefer more `jobs` with fewer `threads` each; for a few large samples, fewer
jobs with more threads.

## Which path to use

| Situation | Use | What it gives you |
| --- | --- | --- |
| One sample, or a few local files | `microsuite trim` / `qc` directly | Simplest; one task at a time |
| Batched trimming on a single machine | `scripts/run_fastp_multiqc.sh --jobs N` | Cross-sample parallelism on one box |
| Reproducible multisample pipeline (HPC/cloud) | Nextflow `amplicon_qiime2` (`--trim`) | Scheduler parallelism, pinned containers, provenance |

## Examples

Single sample, one CLI call:

```bash
microsuite trim --backend fastp \
  --read1 reads_R1.fastq.gz --read2 reads_R2.fastq.gz \
  --output1 trimmed_R1.fastq.gz --output2 trimmed_R2.fastq.gz \
  --html sample.fastp.html --json-report sample.fastp.json --threads auto
```

Batched trimming on one machine (4 samples at a time, 4 threads each):

```bash
bash scripts/run_fastp_multiqc.sh ACCESSION --input-dir /data/ACCESSION --jobs 4 --threads 4
```

Reproducible pipeline with fastp trimming via Nextflow:

```bash
nextflow run workflows/nextflow/main.nf -profile docker \
  --workflow amplicon_qiime2 \
  --manifest manifest.tsv --metadata metadata.tsv --classifier classifier.qza \
  --trim true
```
