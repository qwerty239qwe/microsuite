# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.4.1] - 2026-08-20

### Fixed

- MaAsLin 3 now accepts character categorical fixed effects with more than two
  levels after metadata crosses the TSV boundary. Modeled categorical columns
  are restored as factors; the first sorted level is the default baseline.
- `diff_abundance(..., backend="maaslin3", reference="column,level")` and the
  matching CLI `--reference` option select explicit baselines, including
  semicolon-delimited references for multiple columns. Syntax, duplicate
  columns, metadata columns, formula membership, categorical types, and levels
  are validated instead of silently falling back to another baseline.
- MaAsLin 3 factor conversion and reference logging are limited to modeled
  fixed effects, so unrelated sample metadata is no longer retyped or reported
  as if it affected a contrast.

## [0.4.0] - 2026-08-18

### Added

- MaAsLin 3 1.4.0 as a new, independently containerized differential-abundance
  backend for 0.4.0. It accepts a complete lme4 formula or separate fixed and
  random formula terms, preserves the corrected `TSS` + `LOG` defaults, and
  exposes `normalization`, `transform`, `min_prevalence`, and `min_abundance`.
  The output directory retains MaAsLin 3's upstream `all_results.tsv` and adds
  distinct `abundance_results.tsv` and `prevalence_results.tsv` contracts.
- A hardened LEfSe 1.22.0 backend contract: deterministic seeds, explicit
  reference-class ordering, optional crossed subclass/block designs, configurable
  Kruskal-Wallis, Wilcoxon, LDA, and p-adjustment thresholds, strict abundance
  scale/matrix validation, stable `features`/`scores` output, and a parameter
  manifest recording score orientation. The dedicated Bioconductor 3.23 image
  is version-pinned and its build runs a real planted-signal smoke test.

### Fixed

- LEfSe no longer accepts CLR or negative values, silently chooses factor order
  from R's TSV parsing, or relies on `lefser`'s random tied-value handling without
  a seed. Invalid or nested-only subclass layouts now fail before analysis rather
  than yielding an empty or misleading blocked comparison.
- Balanced class-by-subclass designs no longer crash in `lefser` 1.22.0 when its
  internal `apply()` simplifies equal-length sample-index vectors into a matrix;
  the wrapper preserves those vectors as a list and runs the same Wilcoxon
  consistency algorithm.
- The original direct R/container entrypoint
  `lefse.R counts.tsv metadata.tsv group_col output.tsv` remains accepted and
  maps to the hardened defaults, alongside the new parameter-JSON entrypoint.
- Manual heavy-image validation no longer cancels the automatic Docker CI run
  for the same `main` commit; concurrency remains enabled independently for
  push, pull-request, and manually dispatched runs.
- GitHub Actions use their Node 24-compatible releases, removing deprecated
  Node 20 runtimes and the retired `setup-java` v4 integration from CI.

## [0.3.0] - 2026-08-17

### Added

- `microsuite diversity adonis` / `microsuite.api.adonis2`: formula-based
  multi-term PERMANOVA. Wilkinson formulas with `+`, `:`, `*`, and `/`, sequential
  (type I) or marginal (type II/III) sums of squares via `--by terms|margin`,
  and a `(1 | group)` term that restricts which permutations are drawn rather
  than fitting a variance component. Unlike `permute::how()`, group exchange
  also works when group sizes are unbalanced (groups swap only with groups of
  the same size). `tests/test_adonis.py` carries an opt-in parity test
  (`Rscript` on `PATH`) comparing native `Df`, `SumOfSqs`, `R2` and `F`
  against `vegan::adonis2` 2.7.5 for an additive multi-term model and a
  model with an interaction.
- Vegan-backed formula beta-diversity significance via `adonis2` and the
  microsuite `anosim2` compatibility entry point (`vegan::anosim`), with
  blocked permutations through `--strata`.
- A dedicated `r-ecology` Docker image with build-time adonis2/anosim2 smoke
  tests and heavy-image CI registration.
- A QIIME 2 `adonis` backend through `diversity_test`, with formula,
  permutation-count, and parallel-job forwarding; it writes QIIME `.qzv`
  visualizations and does not expose restricted `strata` permutations.
- `microsuite batch correct` — batch effect correction with five backends:
  `mmuphin` (default), `combat-seq`, `conqur`, `plsda-batch`, and `metadict`,
  each in its own container image. Corrected tables record their scale in
  `uns["microsuite"]["value_type"]` as `counts`, `relative`, or `clr`.
- Count-requiring commands (`diff_abundance --backend ancombc/aldex2`,
  `rarefy`, `normalize`) now refuse tables whose recorded scale they cannot
  consume. Tables without a recorded scale are unaffected, so no existing
  pipeline changes behaviour.

### Fixed

- R/DADA2 denoising now accepts `--error-estimation-function noqual` for
  archived FASTQs whose uniform quality scores cannot support the default
  quality-dependent LOESS error model; the resolved estimator is recorded in
  the DADA2 parameter manifest.
- `batch correct`'s default `--runtime docker` image name is now a
  declared field per backend (`BatchBackend.image`) rather than derived by
  interpolating `--backend` verbatim, fixing a broken image pull for
  `combat-seq` and `plsda-batch` (their images build as `r-batch-combatseq`
  and `r-batch-plsdabatch`, without the hyphen).
- `mmuphin` and `metadict` now TSS-normalize their corrected output so the
  `value_type="relative"` label they record is actually true of the data
  (previously both likely returned count-scale values under a `relative`
  stamp; see docs/batch_correction.md Section 4).
- `batch correct` now rejects NA values in `--batch-col`, `--covariates`,
  and `--target-col`, and rejects a backend response containing NA or
  (for `counts`-declared backends) non-integer values, instead of silently
  propagating them downstream.
- `metadict.R` now refuses to run if `--batch-col` names a column other
  than a pre-existing `batch` column in the metadata, instead of silently
  overwriting that column.
- `adonis2`'s `(1 | group)` restricted permutation now raises instead of
  silently returning `p = 1.0` when every plot has a distinct size, so no
  between-plot exchange is possible; it warns when the achievable exchange
  count is coarser than the requested `--permutations`, and reports the
  achievable count as an output column.
- `adonis2` now raises when the distance matrix contains sample IDs missing
  from metadata (naming the missing IDs), instead of silently analysing
  whichever samples happened to match and reporting a smaller `n_samples`.
- `adonis2`'s `--within` is now validated up front and rejected when there
  is no `(1 | group)` term or `--blocks` to shuffle within, instead of
  being silently ignored.
- Vegan `adonis2`/`anosim2` now validate restricted-permutation feasibility:
  all-singleton `--strata` designs are rejected, partial singleton blocks and
  coarse permutation spaces warn, and output records both
  `requested_permutations` and `effective_permutations`. Formula columns that
  are constant within every stratum also warn because their permutation
  p-values are uninformative. The vegan path now rejects distance-matrix
  samples missing from metadata instead of silently dropping them.

### Note

- Three of the five `batch correct` backends — `conqur`, `plsda-batch`, and
  `metadict` — have R scripts whose signatures are derived from published
  package documentation, not from running the packages against real data;
  no container engine was available while they were written. See "How
  proven is each backend" in docs/batch_correction.md before trusting their
  early outputs.

## [0.2.2] - 2026-08-03

### Fixed

- **Released `network infer --backend sparcc` outputs were mislabeled.** In all
  released versions, the backend calculated CLR-Pearson correlations while
  labeling them as SparCC. It now uses the native SparCC estimator with
  Dirichlet normalization and iterative pair exclusion. **Rerun every network
  previously produced with `--backend sparcc`;** corrected edge weights can
  differ. The backend does not calculate bootstrap significance, so `p_value`
  remains `NaN` and must not be interpreted as evidence of significance.

## [0.2.1] - 2026-08-02

### Fixed

- **`diff_abundance --backend maaslin2` produced results confounded by
  sequencing depth.** The backend passed a raw count table to MaAsLin 2 with
  `normalization = "NONE"`, so every feature in a deeply sequenced sample
  carried a systematically larger value and the model attributed that depth
  difference to whatever covariate correlated with it. Output remained
  well formed with plausible p-values, so nothing surfaced the error. Now uses
  `normalization = "TSS"`, MaAsLin 2's own default.
- **`diff_abundance --backend lefse` had the same defect.** LEfSe expects
  relative abundances; handed raw counts, `lefser` emits a warning and
  continues, yielding LDA scores driven by library size. Counts are now
  converted with `lefser::relativeAb` before testing.

`aldex2` and `ancombc` were checked and are unaffected — both handle
compositionality and library size internally.

**If you ran either backend on 0.2.0 or earlier, re-run those analyses.** The
severity depends on how much sequencing depth varies across your samples and
whether that variation correlates with the group being tested.

## [0.2.0] - 2026-07-30

### Added

- **mothur backends.** Three new entry points covering mothur's MiSeq SOP:
  - `microsuite cluster --backend mothur` — alignment-based OTU clustering via
    twelve chained mothur commands (`unique.seqs` → `align.seqs` → `screen.seqs`
    → `filter.seqs` → `pre.cluster` → `chimera.vsearch` → `remove.seqs` →
    `dist.seqs` → `cluster` → `make.shared` → `get.oturep` → `degap.seqs`),
    producing a feature table, a `.shared` sidecar, and unaligned
    representative sequences.
  - `microsuite tax_classify --backend mothur` — naive Bayes classification via
    `classify.seqs`, with per-OTU consensus taxonomy via `classify.otu` when an
    OTU list is supplied. Takes a reference FASTA and taxonomy file pair
    (`--taxonomy-reference` / `--taxonomy-map`) rather than `--classifier`.
  - `microsuite workflow mothur` — the full SOP from a paired FASTQ directory
    to an OTU table plus taxonomy.
- `containers/mothur` — a mothur 1.48.5 image, added to the container build
  matrix.
- A `mothur-smoke` CI job that runs the real twelve-command pipeline on a
  generated two-sample dataset and asserts on the resulting community
  structure, not merely on exit status.
- Reference-data guidance in `docs/mothur.md`. mothur's reference alignment and
  trainset are user-supplied by design, matching the kraken2 and metaphlan
  convention.

### Fixed

- `runtime.run_command` decoded subprocess output as strict, platform-locale
  text. Tools emitting bytes that are not valid UTF-8 — mothur's `align.seqs`
  progress bar, among others — aborted the run with `UnicodeDecodeError`
  instead of executing, and identical tool output decoded differently on
  Windows and Linux. Output is now decoded as UTF-8 with replacement. This
  affects every backend, since all of them route through this function.

### Notes

Sample identity in the mothur backend comes from the count table produced by
`make.contigs`, not from sequence labels. This differs from the `vsearch` and
`usearch` backends, which parse sample IDs out of FASTA headers via
`--sample-delimiter` / `--sample-field`. Pass `--count-table` when invoking
`cluster --backend mothur` directly; `microsuite workflow mothur` threads it
for you.

## [0.1.0]

Initial release. Never tagged; recorded here for continuity.
