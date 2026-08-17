# LEfSe

microsuite 0.4.0 hardens the existing LEfSe backend around Bioconductor
`lefser` 1.22.0. It remains a two-class biomarker method; use MaAsLin 3 when the
analysis needs multivariable fixed effects or random effects.

## Basic use

```bash
microsuite diff_abundance \
  --backend lefse \
  --table table.h5ad \
  --group treatment \
  --reference control \
  --output lefse.tsv \
  --runtime docker
```

The Python API exposes the same controls through `microsuite.api.diff_abundance`.
The group must contain exactly two non-empty levels. `--reference` fixes the
first factor level; if omitted, microsuite uses the lexicographically first
level and records that decision. Positive scores indicate enrichment in the
recorded `comparison` class and negative scores indicate enrichment in the
recorded `reference` class.

## Controls and reproducibility

The backend exposes:

- `--seed` (default `1234`)
- `--kruskal-threshold` and `--wilcoxon-threshold` (default `0.05`)
- `--lda-threshold` (default `2.0`)
- `--p-adjust-method` (`none`, `holm`, `hochberg`, `hommel`, `bonferroni`,
  `BH`, `BY`, or `fdr`; default `none`)
- `--trim-names` for terminal feature labels

`lefser` contains a random-number step for sparse/tied values. microsuite calls
`set.seed()` immediately before the fit and sorts the result deterministically.
The output schema is always `features` and `scores`, including an empty result.

## Subclass/block designs

`--subclass column` enables LEfSe's Wilcoxon subclass consistency stage. Every
subclass level must occur in both classes: this is a crossed blocking/replicate
factor. It is not a nested design and it is not equivalent to an lme4 random
effect such as `(1 | subject)`. For longitudinal or nested mixed models, use
MaAsLin 3 instead.

## Input and output contract

LEfSe accepts non-negative counts or relative abundances. Counts are converted
with `lefser::relativeAb`; declared CLR tables, negative/non-finite values, and
zero-library samples are rejected before R runs. Supplying terminal taxonomy
nodes avoids the collinearity created by including both ancestors and their
descendants.

`--output lefse.tsv` writes the stable result table and
`lefse.tsv.params.json`, which records the resolved reference/comparison,
subclass, seed, thresholds, and adjustment method. Docker runs additionally
write `lefse_container.json`. `--force` replaces existing results only after a
new result has passed schema validation.
