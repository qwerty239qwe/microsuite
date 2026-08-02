# Pinned SpiecEasi SparCC evidence

> **Generated evidence; never hand-edit.** Regenerate every TSV with the commands
> below, then review and update the recorded hashes and metrics together.

These fixtures are clean-room benchmark evidence for microsuite's native SparCC
implementation. No SpiecEasi source is copied into this repository.

## Oracle provenance

- Repository: `https://github.com/zdk123/SpiecEasi`
- Git commit: `faed6a4476fe0a8dc701ea15cbdfe98d56ce6704`
- `DESCRIPTION` version: `1.99.0`
- License of the external oracle: GPL (>= 3)
- Capture date: 2026-08-02
- R: `R version 4.6.0 (2026-04-24)`
- VGAM: `1.1-14` (reported by `packageVersion()` as `1.1.14`)
- MASS: `7.3-65` (reported by `packageVersion()` as `7.3.65`)
- R RNG: `Mersenne-Twister`, `Inversion`, `Rejection`

`capture_reference.R` refuses to run unless the checkout is exactly the pinned
commit and its `DESCRIPTION` reports version `1.99.0`. It sources only
`R/normalization.R`, `R/mvdistributions.R` (for `cor2cov()`), and `R/spaRcc.R`
from that checkout. `sparcc()` runs with its defaults (`iter=20`,
`inner_iter=10`, `th=0.1`) after `set.seed()` for each capture.

## Regeneration

Run from the microsuite repository root:

```bash
git clone https://github.com/zdk123/SpiecEasi.git /tmp/microsuite-spieceasi-faed6a4
git -C /tmp/microsuite-spieceasi-faed6a4 checkout --detach faed6a4476fe0a8dc701ea15cbdfe98d56ce6704
uv run python tests/fixtures/sparcc/generate_inputs.py --output-dir tests/fixtures/sparcc
Rscript tests/fixtures/sparcc/capture_reference.R /tmp/microsuite-spieceasi-faed6a4 tests/fixtures/sparcc
uv run pytest tests/test_sparcc_reference_fixtures.py -q
```

The Python generator uses `numpy.random.default_rng(10010)`. Both count tables
have 400 samples and ten features. Their latent log-abundance correlation is the
positive-definite matrix in `latent_correlation.tsv`, built from three sparse
factors with strong and moderate positive and negative pairs. Dense sample
depths are log-uniform from 500 to 20,000. The zero fixture starts at depths 250
to 1,000 and deterministically censors the lowest observed counts (seeded random
tie-breaking) until exactly 19% of cells are zero. The committed zero table has
post-censoring depths 204 to 997 (median 644). The inner fixture contains 24
strictly positive six-feature Dirichlet compositions.

Tables are samples by features: the first column is the sample ID and the header
is feature order. Correlation and latent matrices are features by features, with
the same order on rows and columns. Reference TSVs retain R's full double
precision and must not be rounded.

## Frozen measurements

All MAEs use the 45 strict upper-triangle coefficients. The reference median is
the elementwise median of captures at seeds 10010, 10011, and 10012.

| Dataset | Seed pair | Reference-to-reference MAE |
| --- | --- | ---: |
| dense | 10010 / 10011 | 0.0029436360626401118 |
| dense | 10010 / 10012 | 0.0027414467647566330 |
| dense | 10011 / 10012 | 0.0024112772744486320 |
| zero | 10010 / 10011 | 0.0096626842339306830 |
| zero | 10010 / 10012 | 0.0106018845520708500 |
| zero | 10011 / 10012 | 0.0092446187415693710 |

| Dataset | CLR-Pearson to reference median MAE | Reference median to latent truth MAE | CLR-Pearson to latent truth MAE |
| --- | ---: | ---: | ---: |
| dense | 0.09383232047077433 | 0.06363335180611032 | 0.07589352699121227 |
| zero | 0.09221407177885772 | 0.07541929022563038 | 0.08866756901499570 |

Here CLR-Pearson means Pearson correlation of sample-wise
`clr(counts + 1)`, the pre-fix microsuite behavior.

The maximum reference seed MAE is `0.01060188455207085`. The outer parity
criterion was frozen before implementation using the plan's rule:

```text
allowed_mae = max(0.02, 5 * maximum_reference_seed_mae)
            = 0.05300942276035425
```

This numeric threshold must not be loosened without recapturing the oracle
evidence and documenting the reason.

Reference-median edges at `abs(correlation) >= 0.3` are:

- dense: `feature_01--feature_02` (0.773271143559),
  `feature_01--feature_03` (-0.593417967208),
  `feature_02--feature_03` (-0.561864623797),
  `feature_03--feature_04` (0.439317225271),
  `feature_05--feature_06` (0.389012626284),
  `feature_05--feature_07` (-0.541920626524),
  `feature_06--feature_07` (-0.452285932850),
  `feature_08--feature_09` (-0.558311018872),
  `feature_08--feature_10` (0.310870215957), and
  `feature_09--feature_10` (-0.353740528014).
- zero: `feature_01--feature_02` (0.625525581116),
  `feature_01--feature_03` (-0.407921373684),
  `feature_02--feature_03` (-0.405004127229),
  `feature_03--feature_04` (0.301574172807),
  `feature_05--feature_07` (-0.481956882276),
  `feature_06--feature_07` (-0.342810273694), and
  `feature_08--feature_09` (-0.331258247379).

## SHA-256

The validation test pins these hashes as well. `README.md` is excluded to avoid
a self-referential hash.

<!-- HASHES_START -->
| File | SHA-256 |
| --- | --- |
| `generate_inputs.py` | `1f807d24b6d5279097fa73de9c6ed19dc67dbb044b5649843aa347940bb4cd63` |
| `capture_reference.R` | `aea4b40cde4d089bfd4b23a6a8704a95811dfdc7ab9acd28c59405d734dbf9a3` |
| `dense_counts.tsv` | `6f3f5d9591ab6a6ad94c7a66948be2dd3cf1ee0e52862a3b2db85b57707d76c8` |
| `zero_counts.tsv` | `741273ed55fd75e5f07bcce6ff7b3d3efee4342faecd828e1a281b5974c940da` |
| `inner_compositions.tsv` | `dd0e611b85863fe88a7f2e66a9c30b014032d9e9c530bb1b97ffe5f4ecbd79d9` |
| `latent_correlation.tsv` | `b1aa6e2b318cb94efa622a7328492626241390b15553651ef73cd5ac9101be38` |
| `inner_initial_reference_cor.tsv` | `ec12047a7c9954f3be4e421714158120143faf86c872805ca09dbdb56d0703b9` |
| `inner_reference_cor.tsv` | `bac0fb0f1e291a7fd0418598101717ee8683ebda1a30a4f2b1316af06ba2750f` |
| `dense_reference_cor_seed_10010.tsv` | `9f4cf0eaf3715bb3d9620ba1973754c880aefbcaf3ccdfb7449534e750dda007` |
| `dense_reference_cor_seed_10011.tsv` | `2e84673c572f48de0bf02feeb1d1984abed3a138481acb51c335104d9a2f1e42` |
| `dense_reference_cor_seed_10012.tsv` | `7edf4d3ed2c0ab10f613e4a59410509a4852b503b1b855073a9aaeece4e8a0fb` |
| `zero_reference_cor_seed_10010.tsv` | `0661bc47815604fb91edd7945801bde23974dedcb12c6e559804254510f9c335` |
| `zero_reference_cor_seed_10011.tsv` | `ec46876418136c32487ac7d4dab2ff4e7286e8645563c83c5987c4213bc4dbd6` |
| `zero_reference_cor_seed_10012.tsv` | `2fe541ab65bfab9aa81875c344b3db4a57df729428a2dbbd66f979a5a254af7b` |
<!-- HASHES_END -->
