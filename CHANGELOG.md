# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0] - 2026-08-05

### Added

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
