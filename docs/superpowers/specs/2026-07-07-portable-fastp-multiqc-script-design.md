# Portable `run_fastp_multiqc.sh` Script — Design

- **Date:** 2026-07-07
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** A codex-authored `run_fastp_multiqc.sh` (living in the user's
  Dropbox) batched fastp trimming by looping the single-sample `microsuite trim`
  CLI, then summarized with MultiQC. It was accurate about the microsuite CLI but
  environment-bound: it scraped fastp `.deb`s via `apt`/`dpkg-deb` into `/tmp`,
  hardcoded a WSL `/mnt/c/...` path, bootstrapped a `/tmp` uv env, defaulted to
  overwriting outputs, and ran samples sequentially. This sub-project (A of
  three) brings a portable, parallel, reviewed version into the repo.

## Scope

Sub-project **A** of three (B = Nextflow fastp module, already merged; C =
multisample API/CLI/workflow docs + ergonomics). This spec covers **A only**: a
single self-contained bash script plus its tests. No Python package changes.

### Out of scope for A
- Any change to the `microsuite` CLI or Python package.
- The Nextflow workflow (sub-project B, done).
- Documentation of the shared manifest schema / concurrency guidance (C).
- Actually running fastp/multiqc in CI (those stay opt-in/manual; the real
  trimming is exercised by the user on real data).

## Repo conventions this must follow (verified)

- `docs/installation.md` documents that `trim --backend fastp` requires `fastp`
  on `PATH`, with Docker images as the alternative for avoiding local installs.
  The script follows this: check tools on `PATH`, error clearly otherwise — no
  deb-scraping.
- The CLI convention (README): "Commands overwrite outputs only when `--force`
  is supplied." The script must match: no overwrite by default.
- `microsuite trim --backend fastp` is single-sample (`--read1/--read2`,
  `--output1/--output2`, `--html`, `--json-report`, `--threads`, `--run-dir`,
  `--force`); `microsuite qc --backend multiqc` takes `--input-dir`,
  `--output-dir`, `--run-dir`, `--force`. These flag names are already correct in
  codex's script and are kept.

## Design

### Location and shape

New file `scripts/run_fastp_multiqc.sh` (new top-level `scripts/` dir),
`#!/usr/bin/env bash`, `set -euo pipefail`.

**Kept from the original:** the argument parser; the inline `python3` manifest
builder (regex PE/SE grouping with duplicate/orphan detection, writes
`fastq_manifest.tsv` with columns `sample_id / layout / read1 / read2`); the
per-sample `microsuite trim --backend fastp` → single `microsuite qc --backend
multiqc` structure.

**Removed / fixed:**
- Delete the `apt download fastp libisal2` + `dpkg-deb -x` bootstrap and the
  `LD_LIBRARY_PATH=.../x86_64-linux-gnu` export.
- Delete the `/tmp` uv-env bootstrap (`UV_PROJECT_ENVIRONMENT`, `uv sync`,
  `uv pip install multiqc`) — microsuite is required on `PATH`, so the script
  just invokes it.
- Delete the hardcoded `/mnt/c/Users/<user>/...` input fallback and all WSL
  assumptions.
- Flip the overwrite default: today `FORCE=1` (always `--force`). Change to
  no-overwrite by default; add an opt-in `--force` flag and drop `--no-force`.
- Delete the MultiQC temp-dir + `rsync --no-perms` / `cp --no-preserve` dance
  (it existed only to dodge container permission quirks; `cp --no-preserve` is
  GNU-only and breaks on macOS). MultiQC writes straight to the final dir.

### Dependency preflight

Near the top (after arg parsing, before work), check `command -v` for
`microsuite`, `fastp`, and `multiqc`. If any is missing, print an actionable
message naming the missing tool and pointing to `docs/installation.md` (fastp:
"install fastp or use the microsuite fastp container"), then exit non-zero.

### CLI surface

```
run_fastp_multiqc.sh ACCESSION [options]
  --input-root DIR     Root containing accession dirs (default: ./data)
  --input-dir DIR      FASTQ dir to use directly (overrides --input-root)
  --output-root DIR    Output root (default: ./results)
  --jobs N             Samples to trim concurrently (default: 1)
  --threads T          fastp threads per sample (default: 4)
  --force              Overwrite existing outputs (default: off)
  --manifest-only      Detect layout, write the manifest, then exit
  --help
```

Help text notes: total cores used ≈ `jobs × threads`; `--jobs 1` (default)
reproduces one-sample-at-a-time behavior.

### Cross-sample parallelism

Each sample's fastp run is a self-contained unit. After the manifest is written,
the script feeds the manifest's data rows to `xargs -P "$JOBS" -L 1` invoking a
worker (`run_one_sample`, a bash function in the same file, called via
`bash -c '... run_one_sample "$@"' _ <fields>` with the needed vars exported), so
up to `JOBS` samples trim concurrently, each with `--threads "$THREADS"`. Each
sample still writes its own `logs/${sample}.fastp.{stdout,stderr}.log`. `xargs`
returns non-zero if any worker exits non-zero; with `set -o pipefail` the script
propagates the failure. MultiQC runs **once, after** all fastp jobs complete
(it is an aggregation over `fastp_reports/`, not per-sample).

### Portability

- `#!/usr/bin/env bash`; no GNU-only flags. MultiQC writes directly to
  `${OUT_DIR}/multiqc` (no temp-dir/rsync/`cp --no-preserve`).
- The manifest builder stays inline `python3` (portable; python3 is present
  wherever microsuite is installed).
- Works on macOS and Linux.

## Testing

New `tests/test_run_fastp_multiqc_script.py`:

1. **Manifest correctness** — run
   `bash scripts/run_fastp_multiqc.sh acc --input-dir <fixture> --output-root <tmp> --manifest-only`
   against a tiny fixture FASTQ dir holding one PE pair (`x_R1.fastq.gz`,
   `x_R2.fastq.gz`) and one SE file (`y.fastq.gz`); assert the emitted
   `fastq_manifest.tsv` has exactly one PE row (correct read1/read2) and one SE
   row (empty read2), with the right `sample_id`/`layout`. Empty gzipped files
   suffice (manifest detection is filename-based).
2. **Missing-tool guard** — run the script with `PATH` scrubbed of `fastp` (a
   minimal `PATH` containing only python3/bash/microsuite/multiqc stand-ins, or
   monkeypatched via a shim dir) and assert non-zero exit + a stderr message that
   names `fastp` and references installation docs. (If constructing a clean PATH
   is impractical in the environment, assert the guard via a shim directory that
   shadows `command -v` for `fastp` returning empty — implementer picks the
   robust approach; the contract is: missing fastp ⇒ clear non-zero failure.)
3. **Syntax + lint** — `bash -n scripts/run_fastp_multiqc.sh` must pass; and if
   `shellcheck` is on PATH, it must pass with no errors (skip cleanly if
   `shellcheck` is absent).

No test runs fastp or multiqc. CI verifies manifest logic, the missing-tool
guard, and script syntax; the actual trimming is exercised manually on real data.

## Success criteria

1. `scripts/run_fastp_multiqc.sh` exists, is `bash -n`-clean and shellcheck-clean.
2. No `apt`/`dpkg`/`/tmp`-uv/`/mnt/c`/`LD_LIBRARY_PATH`/`rsync --no-perms`/
   `cp --no-preserve` remain.
3. Missing `microsuite`/`fastp`/`multiqc` ⇒ actionable non-zero error citing docs.
4. `--force` is opt-in; absent it, the script does not pass `--force` to microsuite.
5. `--jobs N` trims up to N samples concurrently via `xargs -P`; `--jobs 1` is
   the default and matches prior one-at-a-time behavior; MultiQC runs once after.
6. `--manifest-only` produces a correct PE/SE manifest, asserted by the test.
7. Runs on macOS and Linux (no GNU-only constructs).

## Open questions / follow-ups (not blocking A)

- Whether `--threads` default (4) and `--jobs` guidance belong in sub-project C's
  concurrency documentation (they do; C will cross-reference).
- Whether to also register the script under a `microsuite`-adjacent help pointer
  (C's ergonomics scope).
