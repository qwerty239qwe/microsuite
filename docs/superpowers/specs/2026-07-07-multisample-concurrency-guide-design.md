# Multisample & Concurrency Guide (`docs/multisample.md`) — Design

- **Date:** 2026-07-07
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** While bringing fastp trimming into both the Nextflow workflow
  (sub-project B) and a portable script (sub-project A), the analysis surfaced
  that microsuite has three distinct multisample paths that share one manifest
  format and one concurrency model, but none of it is documented in one place.
  Users (and agents) reinvent batching because the design intent is implicit.
  This sub-project (C of three) captures that knowledge as a doc.

## Scope

Sub-project **C** of three (A = portable script, B = Nextflow fastp module —
both merged). This spec covers **C only**, and C is **documentation only**: one
new guide plus cross-links plus a presence test. No source-code or workflow
behavior changes.

### Out of scope for C
- Any change to the CLI, Python package, script, or Nextflow workflow.
- A new `microsuite manifest` command or CLI `--help` edits (considered and
  deliberately deferred; docs-only was chosen).

## Verified context

- **Shared manifest:** `scripts/run_fastp_multiqc.sh` writes a
  `sample_id / layout / read1 / read2` TSV (layout ∈ {PE, SE}, SE rows have an
  empty `read2`); the Nextflow `amplicon_qiime2` workflow reads a
  `sample_id / read1 / read2` TSV via `splitCsv(header: true, sep: '\t')`. Both
  key by `sample_id` and treat an absent/empty `read2` as single-end. The two
  formats are compatible (the workflow simply ignores the extra `layout`
  column).
- **Concurrency model (from the session analysis):**
  - Within a sample: tool threads — `fastp --thread`, `cutadapt -j`, etc.;
    `microsuite ... --threads auto` = detected CPUs − 1. fastp scales sublinearly
    past ~4–8 threads.
  - Across samples: only two places — `scripts/run_fastp_multiqc.sh --jobs N`
    (`xargs -P`) and Nextflow's scheduler (one process per sample). Plain
    `microsuite trim`/`qc` are single-sample; there is no Python-level
    cross-sample parallelism.
- Docs live under `docs/`; there is no existing multisample/concurrency/manifest
  reference. `docs/api-cli.md` mentions `--threads auto` only in passing.

## Design

### New file: `docs/multisample.md`

Sections, in order:

1. **The sample manifest.** Document the shared TSV: columns
   `sample_id`, `layout`, `read1`, `read2`; `layout` is `PE` or `SE`; SE rows
   leave `read2` empty. State that `scripts/run_fastp_multiqc.sh` produces this
   file (`--manifest-only`) and the Nextflow `amplicon_qiime2` workflow consumes
   the same `sample_id/read1/read2` shape (it ignores the extra `layout`
   column), so a manifest made for one path works with the other. Include a
   small example manifest with one PE and one SE row.

2. **Concurrency model.** Explain the two independent axes:
   - *Within a sample* — tool threads; `--threads auto` = CPUs − 1; fastp's
     sublinear scaling past ~4–8 threads.
   - *Across samples* — only `scripts/run_fastp_multiqc.sh --jobs N` and
     Nextflow parallelize across samples; the bare CLI is one sample per call.
   Give the practical rule: **total cores ≈ jobs × threads**, and for many small
   samples prefer more `jobs` with fewer `threads` each.

3. **Which path to use** — a decision table:

   | Situation | Use |
   |---|---|
   | One sample, or a few local files | `microsuite trim` / `qc` CLI directly |
   | Batched trimming on a single machine | `scripts/run_fastp_multiqc.sh --jobs N` |
   | Reproducible multisample pipeline (HPC/cloud, containers) | Nextflow `amplicon_qiime2` (`--trim`) |

   One line each on what the path gives (parallelism source, reproducibility,
   containers).

4. **Worked examples.** One concrete invocation per path:
   - single CLI: `microsuite trim --backend fastp --read1 ... --output1 ...`
   - script: `bash scripts/run_fastp_multiqc.sh ACC --input-dir DIR --jobs 4 --threads 4`
   - Nextflow: `nextflow run workflows/nextflow/main.nf -profile docker --workflow amplicon_qiime2 --manifest ... --metadata ... --classifier ... --trim true`

### Cross-links (one line each, pointing to `docs/multisample.md`)

- `README.md` — near the "Three ways to use it" table.
- `docs/methods.md` — near the trim/QC entries.
- `docs/installation.md` — near the Nextflow/HPC and Docker rows.
- `docs/api-nextflow.md` — where the manifest/`amplicon_qiime2` workflow is described.

Each cross-link is additive (a sentence or table cell); no existing content is
removed or restructured.

## Testing

New `tests/test_multisample_docs.py` (string/file-presence, matching the repo's
doc/skeleton test convention):

1. `docs/multisample.md` exists and contains the section anchors/headers for the
   manifest schema, the concurrency model, and the decision guide (assert a few
   distinctive substrings, e.g. `sample_id`, `jobs × threads` (or `jobs x
   threads`), `amplicon_qiime2`).
2. `README.md` links to `docs/multisample.md` (assert the relative path
   `multisample.md` appears in the README), so the primary cross-link cannot
   silently rot.

No behavior to test beyond presence; docs are prose.

## Success criteria

1. `docs/multisample.md` exists with the four sections (manifest, concurrency,
   decision guide, examples), the example manifest, and the decision table.
2. The manifest schema and the `jobs × threads` guidance match the as-built
   behavior of `scripts/run_fastp_multiqc.sh` (A) and the Nextflow workflow (B).
3. README and at least `methods.md` / `installation.md` / `api-nextflow.md`
   link to the new guide.
4. `tests/test_multisample_docs.py` asserts the doc exists (with its section
   substrings) and that README links to it; the full offline suite stays green.

## Open questions / follow-ups (not blocking C)

- If a `microsuite manifest build` command is added later (deferred ergonomic
  option), this guide's manifest section becomes its reference and should
  cross-link to it.
