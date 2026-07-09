# DADA2 naming-contract test (P4) — Design

- **Date:** 2026-07-09
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** Codex's ERP120510 friction list, complaint #10 (wrapper/backend
  naming contract not aligned/tested end-to-end). **P4** (final) of the DADA2
  roadmap (see [[dada2-improvement-roadmap]]); P1–P3 merged. Also folds in the
  deferred lowercase-`_r1/_r2` `IGNORECASE` cleanup from P2's review.

## Scope

P4 locks the sample-naming contract between the microsuite `dada2-r` wrapper and
the R backend, via: a shared naming-contract corpus, a CI-runnable Python unit
test that pins `_expected_sample_ids` to it, an opt-in integration test that runs
the **real** `dada2-r` backend and asserts the ASV columns equal the intended
sample IDs, and the `IGNORECASE` fix so lowercase read suffixes pass.

### Out of scope for P4
- Other deferred minors (directory-as-file validation nuance, broad
  FASTQ-artifact substring) — leave as logged follow-ups; P4 stays focused on the
  naming contract.
- Any change to the R backend's own detection logic (it is the source of truth;
  P4 makes the Python side and the tests match it).

## Verified context

- P2 added `_expected_sample_ids(input_dir, *, paired)` in
  `src/microsuite/methods/denoise.py`. `_READ_PATTERNS[0]` (`R[12]`) and `[2]`
  (bare `[12]`) currently lack `re.IGNORECASE`; `[1]` (`read[12]`) and `[3]`
  (`forward|reverse`) have it. The R backend applies `ignore.case=TRUE` globally
  (`dada2_denoise.R:24-33`), so it accepts lowercase `_r1/_r2`; Python does not →
  the deferred drift.
- Single mode: the R backend keeps the FULL FASTQ stem
  (`file_path_sans_ext` twice); paired mode strips the read suffix. P2's
  `_expected_sample_ids` already mirrors this (mode-aware).
- The repo has an opt-in integration convention:
  `pytest.mark.skipif(os.environ.get("MICROSUITE_RUN_EXTERNAL_INTEGRATION") != "1", ...)`
  plus a per-test tool `shutil.which` skip (`tests/integration/test_external_tools.py`).
- A DADA2-**learnable** read generator already exists in
  `.github/workflows/docker.yml` (the `microsuite-dada2` smoke, ~lines 296-318):
  5000 reads, position-dependent jittered quality, quality-correlated
  substitutions, so `learnErrors` fits (uniform/too-clean reads make it return a
  NULL error matrix). P4 reuses this approach.

## Design

### Component 1 — the IGNORECASE fix (`denoise.py`)

Add `re.IGNORECASE` to `_READ_PATTERNS[0]` (the `R[12]` pattern) so paired
lowercase `_r1/_r2` strips like the R backend's global `ignore.case=TRUE`. Add it
to `_READ_PATTERNS[2]` too for uniformity (a no-op there — the bare-digit pattern
has no letters — but keeps the set consistent). No other logic changes.

### Component 2 — shared naming-contract corpus (`tests/naming_contract_cases.py`)

A module of test data, importable by both tests:

```python
@dataclass(frozen=True)
class NamingCase:
    label: str
    filenames: tuple[str, ...]
    paired: bool
    expected: frozenset[str]

CASES: tuple[NamingCase, ...] = ( ... )
```
Cases cover, at minimum: `_R1/_R2`, `_R1_001/_R2_001`, lowercase `_r1/_r2`,
`_read1/_read2`, `_1/_2`, `_forward/_reverse` (all `paired=True`, expected =
stripped sample); and single-end `sample.fastq.gz`, `sample_R1.fastq.gz`
(`paired=False`, expected = full stem, mirroring the R single-mode behavior).
Each `expected` is the set the R backend would produce for that mode.

### Component 3 — CI-runnable Python unit test (`tests/test_naming_contract.py`)

Parametrized over `CASES`: create the empty FASTQ files in `tmp_path`, call
`_expected_sample_ids(tmp_path, paired=case.paired)`, assert
`== set(case.expected)`. Runs in normal CI (no R needed). This is the guard that
would have caught the lowercase-`r` drift; after Component 1 the lowercase cases
pass.

### Component 4 — opt-in real-dada2 end-to-end test (`tests/integration/test_dada2_naming_contract_live.py`)

Gated by `MICROSUITE_RUN_EXTERNAL_INTEGRATION=1`; skips unless `Rscript` is on
PATH (local runtime). A small learnable-read generator (ported from
`docker.yml`'s approach) writes gzipped FASTQs. Two end-to-end cases — the ones
that most exercise the contract:

1. **Single-end** — one sample `sampleS.fastq.gz` (suffix-free) → run
   `denoise(backend="dada2-r", mode="single", demux=dir, ...)`, assert the ASV
   table header columns `== {"sampleS"}`.
2. **Paired** — `sampleP_R1.fastq.gz` / `sampleP_R2.fastq.gz` with **overlapping**
   forward/reverse learnable reads (F = template[:150], R = revcomp of an
   overlapping window) so `mergePairs` succeeds → run
   `denoise(..., mode="paired", ...)`, assert the ASV columns `== {"sampleP"}`
   (the read suffix stripped — the exact case codex's #2 originally broke).

Both read the intended sample id from the shared corpus concept (the expected id
is asserted directly). The test uses `validate=True` (default) so P2's
`_validate_dada2_asv_samples` also exercises the real ASV table. Skips cleanly in
CI (no R) and is never run in the normal suite.

**Primary implementation risk:** DADA2 must produce ASVs on the fixture.
Single-end reuses the proven `docker.yml` generator. For paired, the forward and
reverse reads must overlap enough to merge; the implementer must tune the
generator (overlap window + learnable quality) until a real paired run yields a
non-empty ASV table. If, after genuine effort, a reliably-learnable *paired*
fixture proves infeasible in a tiny test, the paired end-to-end case may be
reduced to single-end only, with paired naming still fully covered by Component
3 — but this fallback must be explicit in the code/report, not silent.

## Success criteria

1. Lowercase `_r1/_r2` paired FASTQs derive the same sample id as `_R1/_R2` (the
   `IGNORECASE` fix), verified by a corpus case in the CI unit test.
2. `tests/test_naming_contract.py` pins `_expected_sample_ids` to the shared
   corpus across all conventions and runs in normal CI.
3. `tests/integration/test_dada2_naming_contract_live.py` runs the real
   `dada2-r` backend and asserts the ASV table's sample columns equal the
   intended IDs; it skips cleanly without `Rscript`/opt-in and never runs in the
   default suite.
4. The full offline suite stays green.

## Open questions / follow-ups (not blocking P4)

- A `--runtime docker` variant of the e2e (run the real dada2 in the container)
  could be added later; P4's e2e uses local `Rscript` for simplicity.
- The remaining deferred minors (directory-as-file validation; broad artifact
  substring) stay as logged follow-ups.
