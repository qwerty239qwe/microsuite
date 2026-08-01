# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
