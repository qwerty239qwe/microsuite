# 0.3.0 — Batch Effect Correction and Diagnosis — Design

- **Date:** 2026-08-05
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Supersedes:** Component 2 of
  `docs/superpowers/specs/2026-08-02-diffab-unification-and-batch-correction-design.md`.
  Components 1, 3, and 4 of that spec are unaffected.
- **Origin:** mBatchNet — "an interactive web server for diagnosis, correction,
  and benchmarking of batch effects in microbiome data", *Bioinformatics* 2026,
  <https://academic.oup.com/bioinformatics/article/42/7/btag538/8739521>,
  source at <https://github.com/gilmore307/mBatchNet>.

## Background

The 2026-08-02 consumer audit found that the consumer project merges 21 run
tables, names primer region, DNA extraction, sequencing platform, cohort, site,
and taxonomy version as confounders in its own QC document, and corrects for
none of them. `run_id` survives only as a provenance column. microsuite has no
batch correction of any kind.

The 2026-08-02 spec answered that with two backends. mBatchNet raises the target:
it benchmarks twelve correction methods, five of them microbiome-specific, and
pairs every correction with a diagnostic panel. This spec adopts the
microbiome-specific five plus ComBat-seq, and the diagnostic half, and drops
benchmarking.

The change in ambition also exposes a defect in the earlier design, described
under "The output-scale hazard" below: a two-value output contract cannot
describe what these six methods actually return.

## Scope

| # | Deliverable |
|---|---|
| 1 | `microsuite batch correct` — six backends, six containers |
| 2 | The `value_type` output contract and the downstream guards that read it |
| 3 | `microsuite batch diagnose` — adonis2 variance partition, pre versus post |
| 4 | Generalizing the R-backend runner out of `diffab/` |

Target release: **0.3.0**, which is written in `CHANGELOG.md` but not tagged and
not published. This work folds into that entry rather than opening 0.4.0.

### Out of scope, with reasons

- **Benchmarking.** mBatchNet's third pillar ranks methods against each other by
  the diagnostic metrics. Once correction and diagnosis both exist, that is a
  loop over the two; building it before them would fix the metric surface before
  it has been used once.
- **The six general-purpose baselines** — ComBat, limma `removeBatchEffect`,
  BMC, FSQN, FAbatch, RUV-III-NB. They are in mBatchNet as baselines to beat,
  not as recommendations. ComBat-seq is included because it is the standard
  where downstream methods require integer counts, which is exactly
  microsuite's situation with ANCOM-BC and ALDEx2.
- **MMUPHin `lm_meta`.** Meta-analytic differential abundance, not batch
  correction, and it runs MaAsLin 2 internally. Deferred for the reason given in
  the 2026-08-02 spec.

## Component 1 — `microsuite batch correct`

A new `batch` sub-app with two subcommands, `correct` and `diagnose`, matching
the `diffab`/`network`/`table` pattern in `cli/app.py`. Table in, corrected
table out; any downstream method consumes the result. Keeping correction a
separate verb rather than a flag on other commands makes it composable and
testable on its own.

### Backends

| Backend | Implementation | Covariates | Requires target label | Emits |
|---|---|---|---|---|
| `mmuphin` (default) | Bioconductor MMUPHin `adjust_batch` | yes | no | `relative` |
| `combat-seq` | Bioconductor sva `ComBat_seq` | yes, via `covar_mod` | no | `counts` |
| `conqur` | GitHub `wdl2459/ConQuR` | yes | no | `counts` |
| `plsda-batch` | GitHub `EvaYiwenWang/PLSDAbatch` | no | yes | `clr` |
| `metadict` | GitHub MetaDICT | yes | no | `relative` |
| `debias-m` | pip `debiasm` (PyTorch) | no | yes | `relative` |

`mmuphin` is the default: it is covariate-controlled, needs no outcome label, and
is the method the 2026-08-02 spec already chose.

Passing `--covariates` to `plsda-batch` or `debias-m`, or `--target-col` to a
backend that does not use it, raises through the existing `reject_options`
helper in `methods/_dispatch.py`. Accepting an option and ignoring it produces a
complete, plausible, wrong result — the failure mode this codebase keeps
reproducing.

### Containers

Six per-backend images, following the existing `r-diffab-*` convention:
`containers/r-batch-mmuphin`, `-combatseq`, `-conqur`, `-plsdabatch`,
`-metadict`, and `containers/py-debiasm`. Each carries a build-time smoke run of
a tiny real dataset through the backend script — not a package import — and
fails the build if it cannot produce a non-empty result. The pattern is
`containers/r-diffab-ancombc/Dockerfile`. All six register in the heavy-image CI
job.

Per-backend images rather than one shared image: three of the five R methods
install from GitHub with no release tags, and a broken install in one must not
be able to break the others.

**Those three GitHub installs pin to commit SHAs, not to a branch.** They can
otherwise change under us with no change on our side. The build-time smoke is
what makes such a change loud instead of silent.

### The R script convention

Every backend script takes `counts.tsv metadata.tsv params.json corrected.tsv`
and reads its options from the JSON — the convention `ancombc.R` already uses.
No positional option passing, so adding a backend option never reorders an
argument list.

`py/debiasm_run.py` takes the same four positional arguments inside the
`py-debiasm` image, so the two runtimes present one interface to the caller.

## Component 2 — The output-scale hazard

This is the part of the design most likely to cause a silent wrong result, and
the reason the 2026-08-02 contract needed widening.

The six backends do not agree on what they return. `ComBat_seq` and ConQuR
return integer counts. `adjust_batch`, MetaDICT, and DEBIAS-M return
non-integer abundances. PLSDA-batch returns data in **CLR log-ratio space**.

Downstream methods disagree in turn: ANCOM-BC and ALDEx2 require integer counts;
MaAsLin 2 and LEfSe normalize internally and accept either; `rarefy` and
`normalize` assume counts.

The 2026-08-02 spec proposed a two-value `value_type` of `counts` or
`abundance`. That cannot express CLR. Under it, `plsda-batch` output would be
labelled `abundance`, and `normalize --method clr` on it would CLR-transform
already-CLR data and return a full, plausible, wrong table.

So the contract is three-valued. `batch correct` writes, on the output AnnData:

```python
uns["microsuite"]["value_type"]      # "counts" | "relative" | "clr"
uns["microsuite"]["batch_correct"]   # {backend, batch, covariates, target, image, digest}
```

Exactly these call sites read `value_type` and raise:

| Call site | Rejects |
|---|---|
| `diff_abundance` backend `ancombc` | `relative`, `clr` |
| `diff_abundance` backend `aldex2` | `relative`, `clr` |
| `rarefy` | `relative`, `clr` |
| `normalize --method relative` | `relative`, `clr` |
| `normalize --method clr` | `clr` |

`maaslin2` and `lefse` do **not** check, because both normalize internally and
accept either.

Absent the key — any table written before 0.3.0, or by any command other than
`batch correct` — the check is skipped and behaviour is unchanged. The key is an
assertion when present, never a requirement. The error names the backend that
produced the table and what it emitted, so the message points at the cause
rather than at the symptom.

### Supervised backends and label leakage

`plsda-batch` and `debias-m` fit using the phenotype labels. Correcting a table
with the outcome label and then testing that same outcome inflates significance:
the correction has already moved the data toward the labels it will be scored
against.

Both backends therefore require `--target-col`, and both emit a warning naming
the column. `docs/batch_correction.md` states the hazard in prose, because a
runtime warning scrolls past and a document does not. Neither backend is
removed — this is the method as published and as benchmarked in mBatchNet, and
refusing to ship it would not stop anyone from doing the same thing by hand.

## Component 3 — `microsuite batch diagnose`

Correction without measurement is faith. `diagnose` answers one question in one
table: **did the batch effect shrink, and did the biological signal survive?**

A method that flattens every difference in the data scores perfectly on any
batch metric read alone. So no metric is reported alone.

### adonis2 variance partition — the headline

The primary diagnostic is PERMANOVA variance partitioning via `vegan::adonis2`,
run on a multi-term formula and reported per term:

```
~ batch + group [+ covariates]
```

This already exists in microsuite. `beta_significance(..., backend="vegan",
method="adonis2", formula=...)` returns a per-term `r_squared` column
(`diversity/r/beta_significance.R:107`) through the `r-ecology` image. `diagnose`
composes it; it adds no R code and no image.

`--compare-to CORRECTED.h5ad` is the main mode, not an extra: R² is interpreted
as a pair. The output puts pre, post, and delta side by side, one row per term.

| term | r2_pre | r2_post | delta |
|---|---|---|---|
| batch | 0.184 | 0.021 | −0.163 |
| group | 0.052 | 0.049 | −0.003 |

A correction succeeds when the `batch` row falls and the `group` row holds. Both
rows are always printed, so the pair cannot be read apart.

### Supporting metrics

Native, in numpy/scipy, adding no dependency and no container:

- kNN batch-mixing entropy
- batch silhouette, computed on a `beta_diversity` distance matrix
- per-principal-component ANOVA variance fraction, for batch and for group

Through the existing `r-ecology` image, which gains `r-lme4`:

- pRDA variance partition, via `vegan::rda` and `varpart`
- PVCA

Reused rather than reimplemented: PERMANOVA, ANOSIM, and PERMDISP already exist
in `beta_significance`; ordination already exists under `microsuite ordination`.
`diagnose` calls them.

Output is a tidy TSV, one row per metric per term, with a `stage` column of
`pre`, `post`, or `delta`.

## Component 4 — Generalizing the R runner

`diffab/_runner.py:invoke_r_backend` carries the bind-mount layout, the
caller-UID execution, and the image-digest sidecar that every containerized R
backend needs. It is hardcoded to scripts in `microsuite.diffab.r` and to
`resolve_diffab_image`.

Move it to `runtime/r_backend.py`, parametrized on the script package and the
resolved image. `diffab/_runner.py` becomes a shim, so every existing diffab
caller is unchanged and no diffab test changes.

Without this, `batch/` copies about 110 lines of mount and provenance logic, and
the two copies drift.

`runtime/container.py` gains `resolve_batch_image(backend, override)`, mirroring
`resolve_diffab_image`, with the `MICROSUITE_R_BATCH_<BACKEND>_IMAGE` override
convention the other resolvers already use.

## Code layout

```
src/microsuite/batch/
  __init__.py
  backends.py          # per-backend argument building and the capability table
  diagnostics.py       # native kNN entropy, silhouette, per-PC ANOVA
  r/mmuphin.R
  r/combat_seq.R
  r/conqur.R
  r/plsda_batch.R
  r/metadict.R
  r/prda_pvca.R        # runs in r-ecology
  py/debiasm_run.py    # runs in py-debiasm
src/microsuite/methods/batch_correct.py     # public API
src/microsuite/methods/batch_diagnose.py
src/microsuite/cli/batch_cmd.py
src/microsuite/runtime/r_backend.py         # moved from diffab/_runner.py
```

`batch/backends.py` holds one capability record per backend — covariate support,
target requirement, emitted `value_type`, R package name for the local-runtime
error message. The dispatch reads that record; it does not branch per backend in
prose. Adding a seventh backend is then a table entry, a script, and a
container.

## Error handling

| Condition | Behaviour |
|---|---|
| `--batch-col` missing from `obs` | Raise, naming the column and listing available ones |
| `--covariates` on `plsda-batch` or `debias-m` | Raise via `reject_options`, naming the backend |
| `--target-col` absent for `plsda-batch` or `debias-m` | Raise, explaining that the method is supervised |
| `--target-col` supplied to an unsupervised backend | Raise via `reject_options` |
| A covariate perfectly confounded with batch | Raise before invoking the backend, naming both columns |
| Count-expecting method given a `relative` or `clr` table | Raise, naming the backend that produced it and what it emitted |
| `normalize --method clr` on a `clr` table | Raise |
| `diagnose --compare-to` on tables with different samples | Raise, naming the difference |
| Backend script exits non-zero | Existing `run_command` failure path, unchanged |

## Testing

The 2026-08-02 handover records the lesson that governs this section: mocked
subprocess tests verify only that we construct the commands we *intended*. The
mothur work shipped fourteen defects past a green suite, every one of them a
complete, well-formed, wrong result.

**Unit, mocked subprocess:** command construction for each of the six backends;
`reject_options` behaviour for every unsupported option; every `value_type`
guard, including the absent-key skip.

**Diagnostics, against analytically known answers:** perfectly mixed batches
give maximum kNN entropy and a batch silhouette near zero; perfectly separated
batches give minimum entropy and silhouette near one. These need no container
and no fixture.

**Integration smoke, per backend.** Deterministically generated two-batch data
with a known batch effect and a known biological signal, no committed fixtures,
modelled on `tests/integration/test_mothur_smoke.py`. It asserts on the biology:
adonis2 batch R² falls below a threshold and group R² is retained above one. An
exit-status-only smoke test would pass on a correction that flattened the entire
table.

**Fixtures.** Real captured stdout per backend, never hand-written, under the
"these are evidence, never edit them" rule from `tests/fixtures/mothur/README.md`.
Captured in the configuration the code actually uses, because the mothur work
found four commands that behave differently on toy input than on real input.

**One deliberate-break check.** After the `value_type` guards are wired, feed a
`clr` table to `ancombc` and confirm the raise fires, then revert. A guard that
passes against input it should reject is worse than no guard.

## Documentation

- `docs/batch_correction.md` — new. When correction is appropriate; why it is
  not a substitute for modelling batch as a covariate; the count-versus-relative
  -versus-CLR distinction and which downstream methods care; the supervised
  -backend leakage hazard; how to read the adonis2 pre/post table.
- `docs/methods.md` — the six backends, their capability table, and the
  `diagnose` metric list.
- `CHANGELOG.md` — fold into the unreleased 0.3.0 entry.

The batch-correction document says plainly that correction is not a substitute
for modelling batch as a covariate. The audit shows a project with 21 runs and
no correction of any kind, and the tempting move on acquiring a correction verb
is to correct once and stop thinking about it.

## Implementation order

Each phase ships independently and leaves the tree releasable.

1. **Foundation.** `runtime/r_backend.py` extraction, the `value_type` contract
   and its guards, `resolve_batch_image`, the `batch` sub-app, and the
   `mmuphin` and `combat-seq` backends with their two containers.
2. **The GitHub-sourced backends.** `conqur`, `plsda-batch`, `metadict`. Highest
   build risk, so they follow a foundation that is already proven.
3. **`debias-m`.** A PyTorch image and a second runtime family.
4. **`batch diagnose`.** adonis2 composition first, then the supporting metrics.

Phases 1–2 and phases 3–4 are likely two implementation plans rather than one.
