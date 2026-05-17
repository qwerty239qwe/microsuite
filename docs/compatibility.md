# Compatibility

Core installation supports TSV import and downstream analysis without heavy
ecosystem dependencies.

Optional extras:

- `biom`: enables BIOM import through `biom-format`
- `qza`: enables QIIME2 artifacts that contain BIOM feature tables
- `r`: documents ANCOM-BC support through external `Rscript`

If an optional dependency is unavailable, the related CLI command exits with an
actionable error instead of breaking core workflows.

The real Moving Pictures downloader uses official QIIME 2 tutorial artifacts:
`sample-metadata.tsv`, `table.qza`, and `taxonomy.qza`. Importing the downloaded
feature table requires the `qza` extra because the artifact contains a BIOM table.

QIIME 2 support intentionally avoids requiring a full QIIME 2 installation.
`microsuite qiime inspect` and `microsuite qiime extract` operate directly on the
artifact ZIP structure and work for `.qza` and `.qzv` files. Feature table import
supports QIIME 2 artifacts whose payload contains BIOM or TSV table data.
