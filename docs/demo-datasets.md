# Demo datasets

## Moving Pictures

The preferred real-world demo dataset for 0.1.0 is the QIIME 2 Moving Pictures
tutorial dataset. It is small enough for quick demos, but rich enough to show
metadata-aware microbiome workflows across body sites, subjects, and timepoints.

Fetch the tiny bundled fixture:

```bash
microsuite data fetch moving-pictures -o data/moving-pictures
```

Fetch the real public tutorial artifacts:

```bash
microsuite data fetch moving-pictures -o data/moving-pictures-real --full
```

The full download writes:

- `sample-metadata.tsv`
- `table.qza`
- `taxonomy.qza`

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
