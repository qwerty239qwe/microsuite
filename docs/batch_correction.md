# Batch effect correction

`microsuite batch correct` adjusts a feature table so that samples processed
under different technical conditions become more comparable, before they are
merged into one analysis. It wraps five R backends behind a single CLI
surface and stamps every corrected table with the scale it now holds, so that
downstream commands can refuse to misread it.

This document explains when correction is the right tool, when it is the
wrong one, how to read the guardrails the code enforces, and how to tell
whether a correction actually worked.

## 1. When correction is appropriate

Batch correction is for the case where **the thing you want to compare and
the thing that changed the numbers are different columns in your metadata**.
Sequencing runs, flow cells, extraction kits, library-prep batches,
processing centers, and sometimes whole cohorts each imprint a technical
signature on a feature table that has nothing to do with biology. When you
merge tables produced under several such conditions into one analysis, that
signature is now mixed in with whatever biological effect you actually care
about.

The concrete case this project was built for: an audit merged **21 run
accessions** into a single feature table, with `run_id` retained in sample
metadata purely as provenance — not as something anyone intended to model or
report on. Twenty-one runs is enough that "just eyeball it" stops being a
plan. `run_id` is exactly the kind of column that should predict almost
nothing about a sample's biology and yet, unmerged and uncorrected, often
predicts a great deal about its read counts.

Correction is appropriate when:

- Batch is **confounded with nothing you care about** — every batch level
  contains a representative mix of the biological groups under study, not
  one batch per group. (`run_batch_correction` checks this and raises if a
  covariate is perfectly confounded with batch; it cannot check whether your
  *outcome* is confounded with batch, because the outcome may not be passed
  as a covariate at all.)
- You need **one merged table** to run diversity, ordination, or
  differential-abundance analyses across all batches at once, rather than a
  per-batch analysis that never gets pooled.
- The technical variable is not the effect under study. If `run_id` were
  itself the variable of scientific interest, correcting it away would
  remove the answer, not the noise.

## 2. Why correction is not a substitute for modelling batch as a covariate

Correction and covariate modelling solve related but different problems, and
conflating them is the most common way batch handling goes wrong.

**Correction adjusts the table.** A backend fits a model of batch effects on
the abundances and subtracts (or otherwise removes) its estimate of that
effect, producing a new table. Every downstream command that reads the
corrected table — diversity, ordination, differential abundance — inherits
whatever the correction did or didn't remove, with no further opportunity to
account for uncertainty in that removal.

**Modelling batch as a covariate adjusts the inference.** Where your design
permits it — for example, `microsuite diff_abundance --backend maaslin2` with
a `--fix-formula`-style multivariable model, or `microsuite diversity
beta-significance --backend vegan --method adonis2 --formula "run_id +
group"` — batch enters the statistical model directly, alongside the term
you're testing. The test itself accounts for batch, rather than trusting that
an upstream correction step already removed it perfectly.

These are not mutually exclusive. Where your design allows it, model batch as
a covariate **instead of, or in addition to,** correcting the table. Modelling
costs you nothing extra once you have the metadata column, and it does not
depend on trusting a correction backend's assumptions about what "batch
effect" means for your data.

**The failure this section exists to prevent:** running `batch correct`
once, writing the corrected table to disk, and then never thinking about
batch again for the rest of the analysis. A correction is an estimate, not a
guarantee. Every result downstream of a corrected table should still be read
with the question "did the correction actually work for *this* comparison?"
in mind — see Section 6.

## 3. Backend table

| Backend | Emits | Covariates | Target | Container image |
| --- | --- | --- | --- | --- |
| `mmuphin` (default) | `relative` | yes | no | `r-batch-mmuphin` |
| `combat-seq` | `counts` | yes | no | `r-batch-combatseq` |
| `conqur` | `counts` | **required** | no | `r-batch-conqur` |
| `plsda-batch` | `clr` | **no** | **required** | `r-batch-plsdabatch` |
| `metadict` | `relative` | yes | no | `r-batch-metadict` |

Notes on specific rows:

- **`plsda-batch`** is the only supervised backend in this release. It
  cannot accept `--covariates` (the CLI rejects the option outright) and it
  requires `--target-col`, because the PLSDA model it fits needs an outcome
  label to condition on. See Section 5 before using it.
- **`metadict`** accepts `--covariates`, but internally, an empty covariate
  list is routed through a constant placeholder column rather than passed to
  the R package as-is. See the callout at the end of Section 3 for why.
- All five backends run in their own container image
  (`ghcr.io/qwerty239qwe/microsuite/<container image>` from the table
  above — note that `combat-seq` and `plsda-batch` build as
  `r-batch-combatseq` and `r-batch-plsdabatch`, without the hyphen, matching
  the `containers/` directory names and CI), selected automatically from
  `--backend`; `--runtime docker` uses the image, `--runtime local` expects a
  local `Rscript` with the corresponding R package installed (`microsuite
  batch correct --help` prints the install hint per backend).

Example invocation:

```bash
microsuite batch correct table.h5ad \
  --output table.corrected.h5ad \
  --batch-col run_id \
  --backend mmuphin \
  --covariates body_site \
  --runtime docker
```

### How proven is each backend

Two of the five backends — `mmuphin` and `combat-seq` — wrap R packages
(`MMUPHin`, `sva`) with stable CRAN/Bioconductor releases and documented,
stable function signatures.

`conqur` is the one backend that **requires** `--covariates`. ConQuR is
conditional by construction: it removes batch effects while holding the named
variables fixed, and with nothing to condition on its design matrix is
degenerate. Naming the biological variable of interest is also the right
scientific default — a correction that is not told what to preserve is how the
signal gets removed along with the batch effect. This constraint was found by
running the package, not from its documentation: the container smoke failed
inside `model.matrix` with "contrasts can be applied only to factors with 2 or
more levels".

**The other three — `conqur`, `plsda-batch`, and `metadict` — are sourced
directly from GitHub, with no release tags to pin against.** Their R scripts
in this repository were written from the packages' published documentation
(man pages and, where the man pages left gaps, the package source itself),
**not from running the packages**. No container engine was available during
the implementation of this feature, so none of these three scripts has ever
been executed against the real package. The container build-time smoke
tests (`containers/r-batch-conqur`, `containers/r-batch-plsdabatch`,
`containers/r-batch-metadict`) are, as of this writing, these scripts' first
real execution — but that build only happens on a **manually dispatched**
heavy-image build (`workflow_dispatch` with `build-heavy-containers=true`);
it does not run on ordinary pull-request or push-to-main CI. Whatever
manually dispatched run first builds those images is the first evidence
that they work at all.

If you are deciding how much to trust a result, weight it accordingly:
`mmuphin` and `combat-seq` rest on tested, versioned packages; `conqur`,
`plsda-batch`, and `metadict` rest on a careful reading of documentation that
has not yet been checked against running code by this project. Treat their
early outputs with more scrutiny, and prefer `mmuphin` or `combat-seq` when
either would satisfy your design.

### MetaDICT: two implementation details that affect anyone calling it directly

These two points are specific to the `metadict` backend and matter beyond
microsuite, because they are properties of the upstream `MetaDICT` package
(`BoYuan07/MetaDICT`) itself, not of this wrapper:

1. **MetaDICT mislabels output columns when samples are not already grouped
   by batch.** Reading `R/MetaDICT.R:201-203` at the pinned commit
   `5b052877328c05e7337e4ce2789a9f48fdecbd9b` (see
   `containers/r-batch-metadict/Dockerfile`) in the package source shows
   that `MetaDICT()` builds its corrected count matrix by binding together
   one column-block per batch, in the order batches are first encountered —
   and then reassigns the *original*, unsorted input's column names onto
   that reordered matrix. If your samples are not already contiguous by
   batch, every sample's corrected values come back labeled with a
   different sample's name — the table is complete and well-formed, and
   silently wrong. `src/microsuite/batch/r/metadict.R` works around this by
   sorting samples into batch-contiguous blocks before calling `MetaDICT()`
   and restoring the caller's original sample order by name afterward. If
   you call `MetaDICT()` yourself outside microsuite, you are exposed to
   this unless you replicate that sort/restore step.

2. **`metadict` with no `--covariates` does not call MetaDICT with an empty
   covariate list.** Passing a zero-length covariate vector reaches
   `rbind()` on zero-column metadata frames inside the package, which
   collapses to zero rows and then crashes on the following assignment.
   Rather than surface that crash, `metadict.R` substitutes a single
   constant placeholder column as the sole "covariate" whenever the caller
   passes none, so the package always sees a one-column frame. This keeps
   MetaDICT running with an effectively intercept-only design when you ask
   for none — but this is **microsuite's own inference about the package's
   internals**, not documented or endorsed MetaDICT behaviour, and it has
   not been checked against a real run of the package.

## 4. The output scale contract

Every corrected table records, under `adata.uns["microsuite"]["value_type"]`,
which of three scales it now holds:

- **`counts`** — integer (or count-like) abundances. `combat-seq` and
  `conqur` emit this scale.
- **`relative`** — proportions that sum to roughly 1 per sample.
  `mmuphin` and `metadict` emit this scale.
- **`clr`** — centered log-ratio values. `plsda-batch` emits this scale.

Downstream commands that need one of these scales assert on it before
running, rather than silently computing something meaningless:

| Command | Requires |
| --- | --- |
| `diff_abundance --backend ancombc` | `counts` |
| `diff_abundance --backend aldex2` | `counts` |
| `rarefy` | `counts` |
| `normalize --method relative` | `counts` |
| `normalize --method total-sum` | `counts` |
| `normalize --method clr` | `counts` or `relative` |
| `normalize --method prevalence-filter` | no restriction |

An **unmarked table is never refused.** Every table written before this
feature shipped, and every table produced by a command that isn't `batch
correct`, has no recorded `value_type` — and the guard treats "unknown" the
same as "anything goes." This is deliberate: the contract must not change the
behaviour of any pipeline that predates it.

If you run a guarded command against a table of the wrong scale, you see an
error like this:

```
MicrobiomeSuiteError: normalize --method relative requires a table of type
counts, but this table is 'clr'. It was produced by 'batch correct --backend
plsda-batch', which emits 'clr'. Use a backend that accepts 'clr', or
correct a different way.
```

What to do about it: either pick a different downstream command that accepts
the scale you have (here, skip renormalizing a CLR table — it is already
transformed), or pick a different `batch correct` backend whose output scale
matches what you need downstream. There is no "force past this" flag; the
guard exists because rerunning `normalize --method relative` on a CLR table
produces a plausible-looking float matrix that means nothing.

## 5. Supervised backends and label leakage

`plsda-batch` is a **supervised** backend: it fits its correction model using
the outcome labels you pass through `--target-col`, not just the batch
labels. (A second supervised backend, DEBIAS-M, is planned for a later
release and is not available yet.)

This creates a specific hazard. If you correct a table using `group` as the
target, and then test for a difference in `group` on the corrected table,
the correction step has already seen and used the exact labels the
downstream test is trying to detect. The correction can partially fit to the
outcome itself, and the resulting test statistics run optimistic — the
correction step **inflates significance** for any test of the same variable
it was given as `--target-col`. This is a form of the classic double-dipping
problem: fit once using the labels, test again using the same labels, and
the second test is no longer honest.

What to do instead:

- **Hold the correction out of the test.** If you must use `plsda-batch`,
  target a variable you are *not* going to test downstream — for example,
  correct using a technical covariate as the nominal target if your design
  allows it — or split your data so the correction is fit on a portion that
  is not used to evaluate the outcome.
- **Use an unsupervised backend instead.** `mmuphin`, `combat-seq`,
  `conqur`, and `metadict` do not take target labels at all, and so cannot
  leak your outcome into the correction step. For most exploratory or
  confirmatory work where the biological group is exactly what you plan to
  test, prefer one of these four.

The CLI enforces the mechanical half of this (you cannot omit `--target-col`
for `plsda-batch`, and you cannot pass it to a backend that doesn't accept
it) but it cannot know which variable you plan to test downstream. That part
is a design decision only you can make.

## 6. How to check that correction worked

Until a dedicated `batch diagnose` command lands, check a correction's
effect directly with beta-diversity significance testing, run twice — once
before correction, once after:

```bash
# Before
microsuite diversity beta table.h5ad --metric bray-curtis -o beta.before.tsv
microsuite diversity beta-significance beta.before.tsv \
  --metadata metadata.tsv --backend vegan --method adonis2 \
  --formula "run_id + group" --runtime docker -o adonis.before.tsv

# After
microsuite diversity beta table.corrected.h5ad --metric bray-curtis -o beta.after.tsv
microsuite diversity beta-significance beta.after.tsv \
  --metadata metadata.tsv --backend vegan --method adonis2 \
  --formula "run_id + group" --runtime docker -o adonis.after.tsv
```

Compare the two `adonis.tsv` outputs term by term. The rule is a pair, not a
single number:

- **The batch term's R² must fall.** `run_id` should explain noticeably less
  variance after correction than before. If it doesn't, the correction did
  not do its job.
- **The biological term's R² must hold.** `group` should explain roughly as
  much variance after correction as before — not more, not dramatically
  less.

**A correction that improves the first number by destroying the second is a
failure, not a success.** It is easy to build a transformation that flattens
`run_id`'s R² to nearly zero by also flattening everything else — including
the biology you were trying to preserve. Always read both terms from the
same formula, before and after, and treat "batch shrank *and* group held" as
the only passing outcome. If batch shrinks but group also collapses, try a
different backend, add covariates that better capture the confound, or
reconsider whether these batches should be merged at all.
