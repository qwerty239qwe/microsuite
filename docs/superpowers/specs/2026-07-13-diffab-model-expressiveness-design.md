# ANCOM-BC model expressiveness (Round-4 H) — Design

- **Date:** 2026-07-13
- **Status:** Approved (design), pending implementation plan
- **Origin:** Round-4 complaints **#1** (ANCOM-BC cannot represent repeated-measures
  designs — one `--group`, no `rand_formula`) and **#2** (ANCOM-BC2 controls and
  CPU parallelism hard-coded away). First diffab sub-project (**H**) of round-4;
  I (containerization), J (result contract) follow. See [[microsuite-round4-roadmap]].

## Scope

Make `microsuite diffab ancombc` express real ANCOM-BC2 models: fixed-effect
formulas with interactions, a subject random intercept, explicit factor reference
levels, and the previously hard-coded controls (prevalence/library cutoffs,
structural zeros, adjustment, global/pairwise/trend/dunnet, pseudo-sensitivity,
worker count). Validate the model (columns exist, full-rank design) and record the
resolved configuration in provenance. This makes ANCOM-BC a genuine
repeated-measures option (the payoff: the oral pipeline could retire its bespoke
MaAsLin2 Docker detour — see [[oral-pipeline-microsuite-refactor]]).

### Out of scope for H
- Containerization / `--runtime`/`--image` (#3, #4) — **I**.
- Standardized cross-backend result contract (#5) — **J**.
- Formula parity for the aldex2/maaslin2/lefse backends — separate follow-up.

## Verified current state

`diffab ancombc <table> --group <col>` → `run_ancombc(adata, group, output)` writes
counts.tsv + metadata.tsv and calls `ancombc.R` with **4 positional args**;
`ancombc.R` calls `ancombc2(fix_formula=group, group=group, p_adj="BH", prv_cut=0,
lib_cut=0, struc_zero=FALSE, neg_lb=FALSE, global=FALSE)` — no rand_formula, all
controls hard-coded, prevalence/library filtering forced off.

## Design

### Component 1 — CLI (`cli/diffab_cmd.py`, `ancombc`)

New options (all optional; keep `--group`):
- `--fix-formula TEXT` — R fixed-effects formula RHS, e.g. `"visit_code*hygiene_status"`.
  Default: the `--group` value. `--group` still sets ANCOM-BC2's `group=`
  (structural-zero / global test target).
- `--rand-formula TEXT` — random-effects formula, e.g. `"(1|subject_code)"`.
- `--reference "col=level"` (repeatable) — set a factor's reference (baseline) level.
- Controls, defaulting to **ANCOM-BC2 native defaults** (documented change from the
  old force-all-off): `--prv-cut FLOAT=0.10`, `--lib-cut INT=0`,
  `--struc-zero/--no-struc-zero` (default off), `--neg-lb/--no-neg-lb` (off),
  `--p-adj-method TEXT=BH`, `--global/--pairwise/--trend/--dunnet` (bool, off),
  `--pseudo-sens/--no-pseudo-sens` (default on, ANCOM-BC2 default), `--n-cl INT=1`.
- Existing: `--output/-o`, `--force`, `--run-dir`, `--timeout`.

### Component 2 — Python wrapper (`diffab/ancombc.py`)

`run_ancombc(adata, *, group=None, fix_formula=None, rand_formula=None, reference=None,
prv_cut=0.10, lib_cut=0, struc_zero=False, neg_lb=False, p_adj_method="BH",
global_test=False, pairwise=False, trend=False, dunnet=False, pseudo_sens=True,
n_cl=1, output, run_dir=None, timeout=None)`:
- Resolve `fix_formula = fix_formula or group`; require at least one of them.
- Light pre-check: the columns referenced by `--group`/`--reference` exist in
  `adata.obs` (raise `MicrobiomeSuiteError`); full formula-column + rank validation
  is authoritative in R.
- Write counts.tsv + metadata.tsv + a **`params.json`** (all resolved options).
  Replaces positional args (too many now). Call `ancombc.R <counts> <metadata>
  <params.json> <output>`.

`params.json` keys: `fix_formula, rand_formula, group, reference` (dict col→level),
`prv_cut, lib_cut, struc_zero, neg_lb, p_adj_method, global, pairwise, trend,
dunnet, pseudo_sens, n_cl`.

### Component 3 — R backend (`diffab/r/ancombc.R`)

- Read the 4 args (counts, metadata, params.json, output); parse params via a small
  dependency-free JSON reader or `jsonlite` if available (guard with
  `requireNamespace`); prefer a minimal hand-parse to avoid a new dependency, else
  jsonlite.
- **Relevel** each `reference` factor (`metadata[[col]] <- relevel(factor(...), ref=level)`);
  error if the level is absent.
- **Validate the design**: for `fix_formula` (and rand), `all.vars` must exist in
  metadata (clean stop naming the missing column); build
  `model.matrix(as.formula(paste("~", fix_formula)), metadata)` and stop with a
  clear "rank-deficient / confounded design" message when
  `qr(mm)$rank < ncol(mm)` (catches the phase*time confound #1 hit).
- Call `ancombc2(...)` with all resolved params (rand_formula passed only when set).
- Write the result table (as today) **and** a provenance JSON beside the output.

### Component 4 — provenance (`<output_dir>/ancombc_provenance.json`)

R writes the resolved config: `fix_formula`, `rand_formula`, `group`, `reference`
(as applied), every effective control value, the factor reference levels actually
used (first level per model factor), and `ancombc_version`/`r_version`. This is the
diffab analogue of the dada2 provenance manifest (A).

## Testing

- **Offline (unit):** monkeypatch `subprocess.run`/`shutil.which`; assert
  `run_ancombc` writes `params.json` with the resolved formulas + controls, that
  `--fix-formula` overrides `--group` and `--group` is the default, that a
  `--reference` for a missing obs column raises, and that the R command line is
  `ancombc.R counts metadata params.json output`. CLI smoke via `CliRunner` with a
  stubbed backend.
- **R parse/validate unit (skips without Rscript):** feed a tiny params.json +
  metadata and assert the rank-deficiency check stops on a confounded design and
  passes on a full-rank one (run the validation section via `Rscript -e`), guarded
  like the dada2 live tests.
- **Opt-in live (ANCOMBC importable + `MICROSUITE_RUN_EXTERNAL_INTEGRATION=1`):**
  a small dataset with a subject random intercept runs end-to-end and writes the
  result + provenance; a rank-deficient design errors cleanly.

## Success criteria

1. `diffab ancombc` accepts `--fix-formula`/`--rand-formula` (repeated-measures
   designs) and the ANCOM-BC2 controls, with ANCOM-BC2 native defaults (not
   force-all-off); `--group` still works as before.
2. Configured formula columns are validated to exist, and a rank-deficient
   fixed-effects design is rejected with a clear message before fitting.
3. `--reference` sets factor baselines; the resolved formulas, controls, and
   reference levels are recorded in `ancombc_provenance.json`.
4. Offline plumbing/validation unit-tested; the live ANCOM-BC path covered by an
   opt-in test. Full suite + `ty check` + ruff + format green.

## Open questions / follow-ups (not blocking H)

- MaAsLin2/ALDEx2/LEfSe formula parity (repeated measures for other backends) —
  a later item; the oral MaAsLin2 retirement needs either ANCOM-BC (this) or a
  maaslin2 rand-effect wrapper.
- Trend-test contrast ordering and pairwise output shape feed into J's
  standardized result contract.
