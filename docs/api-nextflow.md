# Nextflow API

The Nextflow API is for full reproducible workflows.

Use it for large or publication-oriented pipelines. The current Nextflow entry
point is a skeleton for the planned `amplicon_qiime2` workflow:

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

The first planned workflow is `amplicon_qiime2`:

```text
manifest -> qiime2 dada2 -> taxonomy -> phylogeny -> diversity -> report
```

Native statistics and AnnData operations should stay in the Python SDK and CLI
backends, not in Nextflow process scripts.

Current status:

- `workflows/nextflow/main.nf` exists.
- local, Docker, and Singularity profiles exist.
- module files are placeholders and are not yet production implementations.
- default tests validate the skeleton statically rather than running Nextflow.
