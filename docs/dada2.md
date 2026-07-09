# DADA2 denoising

microsuite targets DADA2 `1.40.0` on Bioconductor `3.23` for direct R workflows
and q2-dada2 `2026.4.0` for QIIME 2 artifact workflows. The wrapper preserves
the legacy default: single-end DADA2 unless `paired=True`, `--paired`, or
`mode="paired"` is supplied.

## Backend choice

Use `qiime2-dada2` when the input is already a QIIME 2 demultiplexed reads
artifact and the outputs should remain QIIME artifacts. This backend exposes
q2-dada2 `denoise-single`, `denoise-paired`, `denoise-ccs`, and
`denoise-pyro`; PacBio CCS data should use `mode="ccs"`.

Use `dada2-r` when the input is a FASTQ directory and the desired outputs are
microsuite-native files: an ASV table TSV, representative sequence FASTA, and
denoising stats TSV.

## Running DADA2 in a container

The `dada2-r` backend needs R and the `dada2` package. You do **not** have to
install R locally — run the R step in a container instead. There are three ways:

**1. `--runtime docker` (recommended).** microsuite runs locally (pip-installed)
and delegates only the R step to the `r-dada2` container:

```bash
microsuite denoise --backend dada2-r --runtime docker \
  --demux reads_dir --mode paired \
  --output-table asv_table.tsv --output-rep-seqs rep_seqs.fasta \
  --output-stats stats.tsv --output-plot-dir plots
```

Override the image with `--image IMAGE` or the `MICROSUITE_R_DADA2_IMAGE`
environment variable (default: `ghcr.io/qwerty239qwe/microsuite/r-dada2:latest`).

**2. The whole CLI in one container.** The `microsuite-dada2` image bundles the
CLI and the packaged `dada2_denoise.R`, so the entire command runs inside one
container — mount your working directory:

```bash
docker run --rm -v "$PWD:/work" -w /work \
  ghcr.io/qwerty239qwe/microsuite/microsuite-dada2:latest \
  denoise --backend dada2-r --demux reads_dir --mode paired \
  --output-table asv_table.tsv --output-rep-seqs rep_seqs.fasta \
  --output-stats stats.tsv
```

**3. Nextflow.** The `amplicon_qiime2` workflow runs its DADA2 in the QIIME 2
container automatically; see [multisample runs and concurrency](multisample.md).

### Pull the prebuilt image; don't build it

```bash
docker pull ghcr.io/qwerty239qwe/microsuite/r-dada2:latest
```

Building the DADA2 image locally installs `dada2` through BiocManager, which is
slow and network-fragile — low CPU/memory during the Bioconductor downloads can
look like a frozen build. Pull the published GHCR image instead. For a
reproducible run, pin a `sha-<commit>` tag rather than `latest`.

### Container boundaries

Two images back the DADA2 path, and they are not the same thing:

- **`r-dada2`** — R + `dada2` only (no microsuite). This is what
  `--runtime docker` invokes for the R step.
- **`microsuite-dada2`** — the microsuite CLI + R backend in one image. This is
  the "run the whole command in Docker" path (option 2 above).

See [containers.md](containers.md) for how these images are built, published to
GHCR, and pinned.

## Primers, ITS, and multi-run studies

For ITS workflows, remove primers first with an existing trimming backend such
as Cutadapt or Trim Galore, then run DADA2 on the trimmed reads. Do not rely on
fixed-length truncation alone for variable-length ITS amplicons.

Large multi-run studies are not fully automated yet. The recommended pattern is
to denoise each sequencing run separately with identical parameters where
possible, merge the sequence tables, then remove chimeras and assign taxonomy
globally.

## Diagnostics

For q2-dada2 runs, `--output-base-transition-stats` writes the DADA2 base
transition stats artifact. Supplying both `--output-base-transition-stats` and
`--output-base-transition-plot` also runs `qiime dada2 plot-base-transitions`
and writes a `.qzv` visualization.
