# 0.3.0 — Differential-Abundance Unification and Batch Correction — Design

- **Date:** 2026-08-02
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** An audit of the consumer project, a private research project that
  consumes microsuite across ~21 study accessions. Every workaround script,
  container, and written complaint in that repo was traced back to the microsuite
  code it works around. Findings: `.superpowers/sdd/consumer-audit.md`.

## Background

The audit produced 19 findings. One — `maaslin2` running with
`normalization = "NONE"` on raw counts — was a shipped correctness defect and
was released separately as 0.2.1. This spec covers the next tranche.

Two facts from the audit shape everything below.

**The differential-abundance surface is split.** `microsuite diffab ancombc`
exposes `--fix-formula` and `--rand-formula`; `diff_abundance()` accepts a bare
`group: str` and nothing else. The consumer project uses **both surfaces in one
script** — microsuite's ancombc for one method, a hand-written 76-line R script
inside microsuite's own `r-diffab-maaslin2` container for the other, reached by
overriding the entrypoint. They wanted our container, not our wrapper.

So the most-cited gap is not missing modelling power. The power exists; it is
unreachable from the surface people start at.

**There is no batch correction of any kind.** The project merges 21 run tables,
names primer region, DNA extraction, sequencing platform, cohort, site, and
taxonomy version as confounders in `docs/HOMD_ABUNDANCE_QC.md`, and corrects for
none of them. `run_id` survives only as a provenance column.

## Scope

| # | Deliverable |
|---|---|
| 1 | Formulas as a first-class part of `diff_abundance`, across every backend that supports them |
| 2 | `microsuite batch_correct` — a new verb with two backends |
| 3 | Simpson metric disambiguation |
| 4 | Zero-depth sample guard in alpha diversity |

Items 4 and 5 are unrelated to the rest and ride along because both are
silent-wrong-result defects costing roughly an afternoon each, and item 4 is a
*reversal*: it inverts which group appears more diverse.

### Out of scope, with reasons

- **MMUPHin `lm_meta`.** It is meta-analytic differential abundance, not batch
  correction, and it runs MaAsLin **2** internally. Shipping it beside a future
  MaAsLin 3 backend would invite the assumption that it inherits v3's mixed
  effects and prevalence testing. It does not. Deferred until that can be
  surfaced honestly.
- **MaAsLin 3.** Moved to the 0.4.0 plan so the 0.3.0 release can stabilize the
  formula surface, batch-correction contract, and existing backends first. Its
  separate abundance/prevalence outputs and new container need a dedicated
  output contract and integration smoke test.
- **ALDEx2 `aldex.glm`.** ALDEx2 does support model matrices, so it *could* join
  the formula surface. microsuite currently calls `aldex.ttest`, and extending it
  is a separate piece of work with its own validation burden.
- **Diversity ergonomics and provenance** — audit items 3, 5, 6, 7, 11, 12.
  A coherent 0.4.0 of its own: tidy multi-metric alpha output, `--run-dir` on the
  native commands, Sørensen/Aitchison/UniFrac, PCA/t-SNE, alpha plots.
- **Publishing to PyPI** (audit item 15). Real, large, and orthogonal.

## Component 1 — The unified formula surface

`diff_abundance()` gains two parameters:

```python
fix_formula: str | None = None,
rand_formula: str | None = None,
```

`group` remains and keeps working. When only `group` is given, it is sugar for
a single-term fixed formula — existing callers are unaffected. Supplying both
`group` and `fix_formula` is an error rather than a silent precedence rule,
because guessing which one the user meant is exactly the class of quiet wrongness
this codebase keeps producing.

### Backends split by what they actually are

| Backend | Formulas | Mechanism |
|---|---|---|
| `ancombc` | yes | already has them; `params.json` |
| `maaslin2` | yes | migrate to `params.json` |
| `aldex2` | no | `aldex.ttest` is a two-group test, not a regression |
| `lefse` | no | LEfSe is a two-class biomarker method |

This is not an arbitrary tiering. ALDEx2 and LEfSe are genuinely two-group
methods as microsuite invokes them; pretending otherwise by accepting a formula
and ignoring it would be worse than refusing.

Passing `fix_formula` or `rand_formula` to `aldex2` or `lefse` therefore
**raises**, via the `reject_options` helper in `methods/_dispatch.py` — the same
mechanism `tax_classify` uses for mothur-only options.

### The R script convention

`ancombc.R` already takes `counts.tsv metadata.tsv params.json output.tsv` and
reads `fix_formula`, `rand_formula`, `group`, and `reference` from the JSON.
`maaslin2.R` takes a positional `group_col`. Migrate `maaslin2.R` to the
`params.json` convention.

`aldex2.R` and `lefse.R` keep their positional signature. They gain no
parameters, so a params file would be ceremony.

### `diffab ancombc` stays

It is adopted and wired with full option pass-through in the consumer's
`run_differential_abundance.sh:328-340`. It becomes a thin alias over the
unified path, keeping its current flags. No deprecation warning in 0.3.0 —
breaking a working integration to tidy a namespace is not worth it.

## Component 2 — `microsuite batch_correct`

A new verb: table in, corrected table out. Any downstream method consumes the
result. Keeping it a separate verb rather than a flag on other commands means it
is composable and testable on its own.

| Backend | Method | When |
|---|---|---|
| `mmuphin` | MMUPHin `adjust_batch` | Default. ComBat extended to zero-inflated microbiome profiles, covariate-controlled so biological signal is preserved while study effects are removed. |
| `combat-seq` | sva `ComBat_seq` | Negative-binomial model on raw counts, returns integer counts. Genuinely different assumptions, and the standard where downstream methods require counts. |

Both are R backends, following the established per-backend image pattern.

### The output-type hazard, and how it is handled

This is the part of the design most likely to cause a silent wrong result.

`adjust_batch` returns a corrected **abundance** matrix. `ComBat_seq` returns
**integer counts**. Downstream methods disagree about what they want: ANCOM-BC
and ALDEx2 expect counts; MaAsLin 2 and LEfSe normalize internally; the native
`normalize` command assumes counts.

Feeding a non-count matrix to a count-expecting method produces a complete,
plausible, wrong result — the same shape as the 0.2.1 defect.

So the corrected table records what it is. `batch_correct` writes
`uns["microsuite"]["value_type"]` on the output AnnData as `counts` or
`abundance`, alongside the existing provenance.

Exactly these call sites read the key and raise when it says `abundance`:

- `diff_abundance` backends `ancombc` and `aldex2` — both require integer counts
- `rarefy` — subsampling a non-count matrix is meaningless
- `normalize --method clr` and `--method relative` — already-relative input would
  be normalized twice

`maaslin2` and `lefse` do **not** check, because both normalize internally and
accept either.

Absent the key — any table written before 0.3.0, or by any other command — the
check is skipped and behaviour is unchanged. The key is an assertion when
present, never a requirement.

Batch correction is also **not** a substitute for modelling batch as a covariate.
`docs/batch_correction.md` will say so plainly, because the audit shows this
project has 21 runs and no correction of any kind, and the tempting move is to
correct once and stop thinking about it.

## Component 3 — Simpson disambiguation

`diversity/alpha.py` computes `simpson = 1 - dominance(...)` (Gini-Simpson) and
`simpson_d = dominance(...)` (Simpson's D). No microsuite document says so;
`grep -rn simpson docs/` returns nothing. vegan and mothur use the opposite
convention, and the two quantities are monotonically inverted — a misread
reverses which group looks more diverse.

The consumer maintains a private remap table
(`scripts/python/finalize_microsuite_diversity.py:13-14`) precisely because of
this.

- Add `gini_simpson` as an explicit alias for the current `simpson`.
- Document both in `docs/methods.md`, naming the scikit-bio convention microsuite
  follows and the vegan/mothur convention it differs from.
- Keep `simpson` working unchanged. Renaming it would silently change results for
  every existing caller, which is the very failure being fixed.

## Component 4 — Zero-depth sample guard

A sample with zero total counts yields `shannon = NaN` and `sobs = 0.0` in the
same row, with no warning. The consumer strips such samples before microsuite
sees them and records `reason = zero_total_count`
(`scripts/python/prepare_microsuite_metadata.py:23-46`).

The check goes in `alpha_diversity`'s per-sample entry point, once, rather than
in each metric function — every metric degrades on an empty vector and they
should not each carry the same guard. It raises naming the offending samples,
with a message pointing at filtering or rarefaction depth.

Raising rather than warning is consistent with the rest of the codebase, and a
NaN that flows into a downstream mean is not recoverable by the time anyone
notices it.

`beta_diversity` has the same exposure — a zero-depth sample gives an undefined
Bray-Curtis distance — but is left alone here. Its failure is not silent in the
same way (the distance matrix carries visible NaNs rather than a plausible
number), and changing both at once widens the blast radius of a trivial fix.
Recorded as a follow-up.

## Error handling

| Condition | Behaviour |
|---|---|
| `group` and `fix_formula` both supplied | Raise. Do not guess precedence. |
| Formula given to `aldex2` or `lefse` | Raise via `reject_options`, naming the backend and the unsupported option. |
| Formula references a column absent from `obs` | Raise before invoking R, naming the column and listing available ones. |
| Count-expecting method given an `abundance` table | Raise, naming the backend that produced it and what it emitted. |
| Zero-depth samples in alpha diversity | Raise, naming the samples. |
| `rand_formula` on a backend that supports fixed but not random effects | Raise. No backend in this spec is in that position, but the check keeps a future one honest. |

## Testing

Command construction and option rejection are unit-tested with mocked
subprocess, matching the existing R-backend tests.

One thing needs real execution, because the mothur work established that mocked
subprocess tests verify only that we construct the commands we *intended*:

- **A `batch_correct` smoke test** on a synthetic two-batch dataset with a known
  batch effect and a known biological signal, asserting the batch effect shrinks
  and the biological signal survives. A correction that flattens everything would
  otherwise look like success.

Both follow `tests/integration/test_mothur_smoke.py`: deterministically generated
data, no committed fixtures, assertions on the biology.

## Documentation

- `docs/methods.md` — both `batch_correct` backends; Simpson convention; the
  formula-capability table.
- `docs/batch_correction.md` — new. When correction is appropriate, why it is not
  a substitute for modelling batch as a covariate, and the count-versus-abundance
  distinction.
- `CHANGELOG.md` — note that `maaslin2` results from before 0.2.1 are not
  comparable.
