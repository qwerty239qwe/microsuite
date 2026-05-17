# Data Attribution

## Moving Pictures

`microsuite` uses Moving Pictures data only for examples, tests, and demos. The
dataset is not owned by this project.

The scientific source for the Moving Pictures study is:

Caporaso, J. G., Lauber, C. L., Costello, E. K., Berg-Lyons, D., Gonzalez, A.,
Stombaugh, J., Knights, D., Gajer, P., Ravel, J., Fierer, N., Gordon, J. I.,
and Knight, R. Moving pictures of the human microbiome. Genome Biology 12, R50
(2011). https://doi.org/10.1186/gb-2011-12-5-r50

The full demo artifacts fetched by:

```bash
microsuite data fetch moving-pictures --full
```

come from the QIIME 2 Moving Pictures tutorial data URLs. The QIIME 2 tutorial
describes those data as a small subset of the Caporaso et al. study selected for
quick tutorial runs.

The tiny fixture bundled under
`src/microsuite/data/fixtures/moving_pictures_small/` is a miniature test
fixture for CI and examples. It uses Moving Pictures-style sample identifiers and
metadata fields, but is reduced for deterministic tests and is not a replacement
for the original dataset.

When using this demo in documentation, notebooks, or derived examples, cite the
Caporaso et al. paper and mention that the tutorial artifacts are distributed by
QIIME 2.
