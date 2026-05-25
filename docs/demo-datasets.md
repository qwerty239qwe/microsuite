# Demo datasets

## Moving Pictures

The preferred real-world demo dataset for 0.1.0 is the QIIME 2 Moving Pictures
tutorial dataset. It is small enough for quick demos, but rich enough to show
metadata-aware microbiome workflows across body sites, subjects, and timepoints.
The data are not owned by this project.

Citation:

Caporaso, J. G., Lauber, C. L., Costello, E. K. et al. Moving pictures of the
human microbiome. Genome Biology 12, R50 (2011).
https://doi.org/10.1186/gb-2011-12-5-r50

The QIIME 2 Moving Pictures tutorial describes its tutorial files as a small
subset of the Caporaso et al. study selected for quick local runs. See
[docs/data-attribution.md](data-attribution.md) for fixture and download
attribution details.

Fetch the tiny bundled fixture:

```bash
microsuite data fetch moving-pictures -o data/moving-pictures
```

Fetch the real public tutorial inputs and artifacts:

```bash
microsuite data fetch moving-pictures -o data/moving-pictures-real --full
```

The full download writes:

- `sample-metadata.tsv`
- `emp-single-end-sequences.zip`
- `85_otus.qza`
- `ref-taxonomy.qza`
- `table.qza`
- `taxonomy.qza`

`emp-single-end-sequences.zip`, `85_otus.qza`, and `ref-taxonomy.qza` are used
by the raw-read QIIME 2 parity workflow:

```bash
microsuite workflow moving-pictures-qiime2 \
  --out runs/moving-pictures-qiime2 \
  --force
```

This workflow expects `qiime` to be available from an activated QIIME 2
amplicon environment. It keeps downloaded tutorial files local; do not commit
them to the repository.

Import the real feature table with taxonomy:

```bash
uv sync --extra qza --extra dev
uv run microsuite import qza data/moving-pictures-real/table.qza \
  --metadata data/moving-pictures-real/sample-metadata.tsv \
  --taxonomy-artifact data/moving-pictures-real/taxonomy.qza \
  -o runs/moving-pictures-real/table.h5ad \
  --force
```

Inspect or extract an artifact without QIIME 2 installed:

```bash
uv run microsuite qiime inspect data/moving-pictures-real/table.qza
uv run microsuite qiime extract data/moving-pictures-real/taxonomy.qza \
  -o data/moving-pictures-real/taxonomy-payload \
  --force
```

These downloaded files are intentionally treated as local demo data rather than
source files.

## Gut To Soil

The QIIME 2 gut-to-soil tutorial is available as an optional real-data
integration dataset. The fetcher downloads small derived artifacts rather than
raw reads, so it can exercise real QIIME 2 archives without rerunning denoising
or taxonomy classification.

Fetch the tutorial artifacts:

```bash
microsuite data fetch gut-to-soil -o data/gut-to-soil
```

The download writes:

- `sample-metadata.tsv`
- `asv-table-ms2.qza`
- `taxonomy.qza`

Import the feature table with taxonomy:

```bash
uv sync --extra qza --extra dev
uv run microsuite import qza data/gut-to-soil/asv-table-ms2.qza \
  --metadata data/gut-to-soil/sample-metadata.tsv \
  --taxonomy-artifact data/gut-to-soil/taxonomy.qza \
  -o runs/gut-to-soil/table.h5ad \
  --force
```

The optional integration test for this dataset is gated by
`MICROSUITE_RUN_REAL_DATA_TESTS=1` because it downloads files from the public
QIIME 2 tutorial site.
