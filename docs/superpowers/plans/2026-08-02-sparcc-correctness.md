# SparCC Correctness Implementation Plan

> **For Luna:** Execute this plan task by task. Keep each task independently
> reviewable and do not combine the mathematical core, public API wiring, and
> live-reference validation into one commit.

**Goal:** Replace the backend currently labeled `sparcc`—which is exactly
Pearson correlation of pseudocount-CLR values—with a reproducible native
implementation of the SparCC correlation estimator, benchmarked against
SpiecEasi's implementation at pinned commit
`faed6a4476fe0a8dc701ea15cbdfe98d56ce6704`.

**Architecture:** Put all SparCC mathematics in a new dependency-free internal
module, `microsuite.methods._sparcc`. Keep `methods/network.py` responsible only
for AnnData filtering, backend dispatch, and edge-list formatting. The pure
core returns covariance and correlation matrices and knows nothing about
AnnData, pandas, CLI arguments, thresholds for emitting edges, or file I/O.

**Tech stack:** Python 3.11/3.12, NumPy, pytest, existing AnnData/pandas network
surface, R 4.6 + VGAM only in the live SpiecEasi reference job.

**Primary references:** Friedman & Alm, 2012,
<https://doi.org/10.1371/journal.pcbi.1002687>; SpiecEasi
`R/spaRcc.R` at commit `faed6a4476fe0a8dc701ea15cbdfe98d56ce6704`.

## Why this is a correctness fix

The current `sparcc_network()` calls `_filtered_feature_matrix(...,
transform="clr")`, then `_edge_list_from_matrix(..., method="pearson")`. A
benchmark against SpiecEasi showed:

| Dataset | Coefficient MAE | Maximum error | Current / reference truth MAE |
| --- | ---: | ---: | ---: |
| Variable-depth lognormal counts | 0.1040 | 0.2696 | 0.1136 / 0.0330 |
| 18.9% zero-inflated low-depth counts | 0.1067 | 0.1674 | 0.1400 / 0.0752 |

The current result matched CLR-Pearson to machine precision (`< 8e-16`).
SpiecEasi seed-to-seed MAE was only `0.0013–0.0055`, so the disagreement is
algorithmic, not Monte Carlo noise.

## Review of the handover

`docs/superpowers/HANDOVER-2026-08-02.md` contains no SparCC-specific design or
partial implementation. Its immediate-next-step section concerns the separate
differential-abundance/batch-correction 0.3.0 spec. The applicable handover
constraint is its central lesson: tests of intended calls are insufficient;
reference-backed and result-level tests are mandatory. This plan therefore
starts by capturing real SpiecEasi results and ends with a live pinned-reference
CI job.

## Scope and non-goals

In scope:

- faithful SparCC correlation and covariance estimation;
- count validation, Dirichlet zero replacement, basis-variance reconstruction,
  iterative exclusion, and median aggregation;
- deterministic results through an explicit seed;
- the existing edge-list interface, prevalence filtering, edge threshold, and
  `top_n` behavior;
- captured reference evidence and a live SpiecEasi parity job;
- public Python/CLI tuning options for the three SparCC controls.

Out of scope:

- `sparccboot()` permutation/bootstrap p-values; `p_value` remains `NaN`;
- changes to the separate `spieceasi` or `flashweave` backends;
- a new runtime dependency or copying GPL-licensed SpiecEasi source into the
  Python package;
- version bump, tagging, or release publication;
- optimization beyond straightforward vectorized Dirichlet draws and small
  dense matrix algebra.

## Fixed design decisions

1. **Keep the backend name `sparcc`.** This is an in-place correctness repair,
   not a new backend. Do not retain the old CLR-Pearson behavior under that name.
   Users who want it already have `native-correlation --transform clr`.
2. **Clean-room implementation boundary.** Use the paper's equations and the
   behavior contract below. Do not copy SpiecEasi source, comments, variable
   names, or control-flow structure into microsuite. SpiecEasi is GPL >= 3; it
   is a benchmark oracle, not vendored implementation material.
3. **Reproducible by default.** Add `seed=0` to the Python API and
   `--sparcc-seed 0` to the CLI. Construct one local
   `numpy.random.Generator`; never use or mutate NumPy's global RNG.
4. **Reference-compatible defaults.** `iterations=20`, `inner_iterations=10`,
   `exclusion_threshold=0.1`, and `pseudocount=1.0` match SpiecEasi's estimator.
   In SparCC, `pseudocount` is the Dirichlet concentration offset
   `alpha = count + pseudocount`; it is not a pre-CLR deterministic addition.
5. **Raw counts only.** Accept integer-valued arrays even when stored as floats.
   Reject negative, non-finite, genuinely fractional, all-zero-sample, and
   all-zero-feature inputs with `MicrobiomeSuiteError`. Reject fewer than two
   samples or three retained features; the sparse-correlation approximation is
   underdetermined below that boundary.
6. **Preserve the output table schema.** Columns remain `source`, `target`,
   `weight`, `abs_weight`, `p_value`, `method`, `backend`. For this backend both
   `method` and `backend` must be `sparcc`; `p_value` remains `NaN`.
7. **No hidden fallback to CLR-Pearson.** Singular basis systems use
   `numpy.linalg.pinv`, matching SpiecEasi's generalized-inverse fallback.
   Any other invalid numeric state raises clearly.
8. **Fixtures are evidence.** Reference output files are generated by the pinned
   implementation, never edited by hand, and accompanied by provenance and
   regeneration commands.

## Mathematical contract

Given a samples-by-features count matrix:

1. For each outer iteration, draw one Dirichlet composition per sample using
   concentration `counts + pseudocount`.
2. CLR-transform each positive composition and compute the sample covariance
   with `ddof=1`.
3. Build the Aitchison variation matrix
   `T[i,j] = var(clr_i) + var(clr_j) - 2*cov(clr_i, clr_j)`.
4. Solve the SparCC basis-variance linear system. Clamp basis variances below
   `1e-4`, matching the reference.
5. Reconstruct covariance by
   `Cov[i,j] = 0.5 * (V[i] + V[j] - T[i,j])`, symmetrize it, convert to
   correlation, clamp reconstructed coefficients to `[-1, 1]` as the pinned
   reference does, then rebuild covariance from the bounded correlation and
   the basis standard deviations.
6. In each inner iteration, find the largest absolute unexcluded off-diagonal
   correlation. Stop when it is at or below `exclusion_threshold`; otherwise
   exclude that symmetric pair, update the basis system for both features, and
   recompute basis variances and correlations. Use deterministic lexicographic
   tie-breaking.
7. Take elementwise medians of outer correlation and covariance estimates.
   Rebuild the final covariance from the median correlation and the square roots
   of the median covariance diagonal so `corr(final_covariance)` equals the
   returned correlation.

## File map

| Path | Responsibility |
| --- | --- |
| `src/microsuite/methods/_sparcc.py` | New pure NumPy estimator and validation. No AnnData/pandas/CLI. |
| `src/microsuite/methods/network.py` | Filter counts, call estimator, convert correlation matrix to edge list. |
| `src/microsuite/cli/network_cmd.py` | Expose SparCC iterations, inner iterations, exclusion threshold, and seed. |
| `tests/test_sparcc_core.py` | Pure mathematical unit and invariant tests. |
| `tests/test_network_method.py` | Backend/API/CLI/schema/filtering tests. |
| `tests/fixtures/sparcc/` | Pinned count inputs, normalized inner input, reference matrices, capture script, provenance README. |
| `tests/integration/test_sparcc_spieceasi_reference.py` | Gated live comparison against external pinned source. |
| `.github/workflows/ci.yml` | Dedicated real-reference job. |
| `docs/methods.md` | Replace “CLR approximation” claim with the actual algorithm and limitations. |
| `CHANGELOG.md` | Unreleased correctness warning and rerun guidance. |

---

## Task 0: Freeze the SpiecEasi reference evidence

Do this before writing production code. It prevents acceptance thresholds from
being moved to fit the implementation.

**Files:**

- Create: `tests/fixtures/sparcc/README.md`
- Create: `tests/fixtures/sparcc/generate_inputs.py`
- Create: `tests/fixtures/sparcc/capture_reference.R`
- Create: `tests/fixtures/sparcc/dense_counts.tsv`
- Create: `tests/fixtures/sparcc/zero_counts.tsv`
- Create: `tests/fixtures/sparcc/inner_compositions.tsv`
- Create: `tests/fixtures/sparcc/inner_initial_reference_cor.tsv`
- Create: `tests/fixtures/sparcc/inner_reference_cor.tsv`
- Create: three reference correlation files per count dataset, for seeds
  `10010`, `10011`, and `10012`
- Create: `tests/test_sparcc_reference_fixtures.py`

### Steps

- [ ] Write `generate_inputs.py` using `np.random.default_rng(10010)`.
  Generate two 400-sample, 10-feature tables from a positive-definite latent
  correlation matrix containing strong and moderate positive and negative
  pairs. One table has variable multinomial depths; the other has low depths
  and deterministic 18–20% dropout. Also write a small strictly positive
  composition matrix for deterministic `sparccinner` capture.
- [ ] Write `capture_reference.R`. It must fail unless the checked-out
  SpiecEasi source commit equals
  `faed6a4476fe0a8dc701ea15cbdfe98d56ce6704`; source only its normalization,
  covariance utility, and SparCC files; run `sparcc()` with defaults for the
  three seeds; run `sparccinner()` once on the committed composition matrix;
  capture both its initial pre-exclusion reconstruction and its completed
  inner-loop result; and write TSV matrices without rounding.
- [ ] In `README.md`, record SpiecEasi commit, its `DESCRIPTION` version, R
  version, VGAM version, exact commands, SHA-256 hashes, matrix orientation,
  and the rule: “generated evidence; never hand-edit.”
- [ ] Capture the files by executing the real reference code.
- [ ] Write fixture tests that verify shapes, symmetry, diagonal ones, finite
  values, hashes/provenance presence, and exact feature order. These tests do
  not test microsuite yet.
- [ ] Calculate and record before implementation:
  - reference-to-reference off-diagonal MAE across the three seeds;
  - current CLR-Pearson-to-reference-median MAE;
  - reference edge sets at `abs(correlation) >= 0.3`;
  - truth MAE for the generated latent model.
- [ ] Define the outer parity acceptance from reference variability, not from
  candidate output:
  `allowed_mae = max(0.02, 5 * maximum_reference_seed_mae)`.
  Record the resulting numeric value in the README. Do not loosen it later
  without recapturing evidence and explaining why.
- [ ] Run `tests/test_sparcc_reference_fixtures.py`; expected PASS.
- [ ] Commit only evidence and evidence-validation tests.

### Achieved when

The repository contains reproducible, provenance-stamped SpiecEasi ground truth
and frozen numerical acceptance criteria. No production SparCC code has changed.

---

## Task 1: Add input validation and Dirichlet normalization

**Files:**

- Create: `src/microsuite/methods/_sparcc.py`
- Create: `tests/test_sparcc_core.py`

**Interfaces introduced:**

```python
@dataclass(frozen=True)
class SparCCResult:
    covariance: np.ndarray
    correlation: np.ndarray

def estimate_sparcc(
    counts: np.ndarray,
    *,
    iterations: int = 20,
    inner_iterations: int = 10,
    exclusion_threshold: float = 0.1,
    pseudocount: float = 1.0,
    seed: int = 0,
) -> SparCCResult: ...
```

Private helper names may vary, but keep each responsibility separate:
validation, Dirichlet normalization, variation matrix, basis solve, covariance
reconstruction, exclusion loop, and outer aggregation.

### Steps

- [ ] Write failing validation tests for wrong dimensionality, fewer than two
  samples, fewer than three features, NaN/Inf, negative counts, fractional
  values, all-zero samples/features, zero/negative pseudocount, iterations `<1`,
  inner iterations `<1`, and threshold outside `[0,1]`.
- [ ] Test that float arrays containing exact integer values are accepted and
  the caller's array is not mutated.
- [ ] Test Dirichlet draws are positive, rows sum to one, repeat exactly for the
  same local seed, differ for a different seed, and do not touch global RNG
  state.
- [ ] Implement only validation, the result dataclass, RNG construction, and
  Dirichlet normalization. Use vectorized gamma draws followed by row
  normalization; do not loop through `Generator.dirichlet` per row.
- [ ] Run `tests/test_sparcc_core.py`; the validation/Dirichlet subset passes,
  later estimator tests may still be absent.
- [ ] Commit.

### Achieved when

The stochastic input layer has a precise counts contract and deterministic,
isolated RNG behavior. No network/backend code calls it yet.

---

## Task 2: Implement and verify the deterministic algebra

**Files:**

- Modify: `src/microsuite/methods/_sparcc.py`
- Modify: `tests/test_sparcc_core.py`

### Steps

- [ ] Add failing tests for CLR centering and Aitchison variation: row CLR means
  are zero; variation is symmetric, nonnegative within floating tolerance, and
  has a zero diagonal.
- [ ] Add hand-calculated small-matrix tests for the initial basis system and
  covariance reconstruction. Assert `ddof=1`, symmetry, diagonal variances,
  correlation diagonal ones, and coefficients bounded by one.
- [ ] Add a singular-system test proving the generalized-inverse path returns
  finite symmetric output. Patch or instrument `np.linalg.solve` so this test
  proves the fallback was exercised rather than merely accepting any output.
- [ ] Add the deterministic reference test: run the new inner algebra with
  zero exclusions on `inner_compositions.tsv` and compare to the appropriate
  `inner_initial_reference_cor.tsv` capture. Require `rtol=1e-10`,
  `atol=1e-12`.
- [ ] Implement the variation matrix, basis system construction/solve,
  `1e-4` variance floor, covariance reconstruction, symmetrization, and
  correlation clipping.
- [ ] Run the core tests and fixture tests; expected PASS.
- [ ] Commit.

### Achieved when

For the same deterministic composition/variation input, microsuite and
SpiecEasi agree to floating-point precision before iterative exclusions.

---

## Task 3: Implement the iterative exclusion loop

**Files:**

- Modify: `src/microsuite/methods/_sparcc.py`
- Modify: `tests/test_sparcc_core.py`

### Steps

- [ ] Write failing tests that pin selection of the largest absolute
  off-diagonal pair, symmetric exclusion, decrement/update of both involved
  basis equations, lexicographic tie-breaking, threshold stop, and hard stop at
  `inner_iterations`.
- [ ] Test that an already excluded pair cannot be selected again.
- [ ] Test positive and negative strong pairs—the selection criterion is
  absolute correlation.
- [ ] Compare the complete inner-loop output on `inner_compositions.tsv` with
  `inner_reference_cor.tsv` using `rtol=1e-9`, `atol=1e-11`.
- [ ] Implement the loop without porting SpiecEasi control-flow text. Represent
  exclusions explicitly as a symmetric boolean matrix or set of `(i,j)` pairs;
  never depend on flattened R/Fortran indices.
- [ ] Deliberately disable exclusion once and verify the complete-inner parity
  test fails; restore it before committing.
- [ ] Run core tests; expected PASS.
- [ ] Commit.

### Achieved when

The entire deterministic inner estimator, including strong-pair removal, agrees
with the reference. A mutation that skips exclusion is caught.

---

## Task 4: Implement outer Monte Carlo aggregation

**Files:**

- Modify: `src/microsuite/methods/_sparcc.py`
- Modify: `tests/test_sparcc_core.py`
- Create: `tests/test_sparcc_parity.py`

### Steps

- [ ] Add failing tests for same-seed reproducibility, different-seed variation,
  result shapes, finite/symmetric matrices, exact unit correlation diagonal,
  coefficient bounds, and `corr(covariance) == correlation`.
- [ ] Implement `estimate_sparcc`: run the Dirichlet + inner estimator exactly
  `iterations` times, take elementwise medians, and rebuild final covariance
  from median correlation plus standard deviations from median covariance.
- [ ] In parity tests, compare microsuite's fixed-seed result to the median of
  the three SpiecEasi matrices. Require off-diagonal MAE no greater than the
  frozen `allowed_mae` from Task 0 for both dense and zero-containing datasets.
- [ ] Assert the new estimator is not CLR-Pearson: on the dense fixture, its
  reference MAE must be at most one third of the recorded current
  CLR-Pearson-to-reference MAE.
- [ ] Assert result-level biology on the generated latent model: new SparCC
  truth MAE is at least 10% lower than CLR-Pearson truth MAE, and all strong
  planted pairs have the correct sign. This guard band was frozen after Task 0
  showed the pinned reference median improves truth MAE by about 15–16%; use
  frozen count input and never regenerate random data inside this test.
- [ ] Deliberately replace the estimator result with CLR-Pearson and verify the
  parity and truth tests fail; restore it.
- [ ] Run core + parity tests; expected PASS.
- [ ] Commit.

### Achieved when

The public pure estimator is reproducible, internally consistent, materially
closer to SpiecEasi than CLR-Pearson, and recovers the planted associations.

---

## Task 5: Wire the estimator into the network backend

**Files:**

- Modify: `src/microsuite/methods/network.py`
- Modify: `tests/test_network_method.py`

### Steps

- [ ] Replace the existing label-only SparCC test with a failing numerical test
  over a small committed count table. Prove the backend calls
  `estimate_sparcc`, not `_transform_counts(..., "clr")` or `pearsonr`.
- [ ] Add a dedicated `_filtered_feature_counts` path. Preserve
  `_filtered_feature_matrix` unchanged for native correlation. Apply prevalence
  filtering on raw counts before SparCC validation and preserve retained
  AnnData feature order.
- [ ] Add a correlation-matrix-to-edge-list helper rather than routing SparCC
  through `_edge_list_from_matrix`, which recomputes Pearson correlations and
  p-values. It must emit each upper-triangle pair once, apply
  `min_abs_weight`, sort identically, honor `top_n`, use `method="sparcc"`,
  `backend="sparcc"`, and `p_value=NaN`.
- [ ] Extend `sparcc_network()` and `network()` with the four tuning arguments
  and forward them unchanged.
- [ ] Test prevalence filtering, feature names/order, threshold boundary,
  empty result schema, `top_n`, output ordering, parameter forwarding, and
  deterministic repeated calls.
- [ ] Test fractional/normalized AnnData fails with an actionable “raw counts”
  message instead of silently estimating a network.
- [ ] Run `tests/test_network_method.py`, core, and parity tests; expected PASS.
- [ ] Commit.

### Achieved when

`network(backend="sparcc")` emits edges from the actual SparCC correlation
matrix, while `native-correlation` behavior remains byte-for-byte unchanged.

---

## Task 6: Expose the reproducibility and tuning contract in the CLI

**Files:**

- Modify: `src/microsuite/cli/network_cmd.py`
- Modify: `tests/test_network_method.py`
- Modify if API-signature assertions require it: `tests/test_api_facade.py`

**New options:**

- `--sparcc-iterations` (default `20`, minimum `1`)
- `--sparcc-inner-iterations` (default `10`, minimum `1`)
- `--sparcc-exclusion-threshold` (default `0.1`, range `[0,1]`)
- `--sparcc-seed` (default `0`)

### Steps

- [ ] Add CLI tests proving defaults and explicit values reach `network()`.
  Prefer monkeypatching the module-level function and inspecting typed values;
  this task tests Typer wiring, not mathematics.
- [ ] Add an end-to-end CLI test that writes two files with the same seed and
  confirms identical bytes, then uses a different seed and confirms at least
  one unrounded weight differs when `--min-abs-weight 0`.
- [ ] Update `--pseudocount` help to state its backend-specific meanings and
  that SparCC requires a value greater than zero.
- [ ] Add the four options and forward them. Do not expose bootstrap/p-value
  options in this task.
- [ ] Run CLI/network tests; expected PASS.
- [ ] Commit.

### Achieved when

Users can reproduce and tune the actual estimator through both Python and CLI,
and unchanged calls retain stable defaults.

---

## Task 7: Add a live pinned SpiecEasi result-level check

Captured fixtures protect default tests; this task proves the capture and the
Python implementation still agree with the real external source.

**Files:**

- Create: `tests/integration/test_sparcc_spieceasi_reference.py`
- Create: `tests/integration/run_spieceasi_reference.R`
- Modify: `.github/workflows/ci.yml`

### Steps

- [ ] Write the integration test gated by
  `MICROSUITE_RUN_SPARCC_REFERENCE=1`, `Rscript`, VGAM, and
  `SPIECEASI_SOURCE_DIR`. It runs real SpiecEasi on both committed count tables,
  parses the matrices, and applies the same frozen MAE plus edge-sign checks as
  the fixture parity test.
- [ ] Make the R runner verify the exact source git commit before sourcing any
  file. It must print R, VGAM, and SpiecEasi versions/commit into captured test
  output for diagnosis.
- [ ] Add a dedicated `sparcc-reference` CI job:
  - checkout microsuite;
  - checkout `zdk123/SpiecEasi` at the exact commit into a separate directory;
  - set up Python and R;
  - install VGAM (do not install or vendor all of SpiecEasi);
  - run only the gated integration test with the source directory environment
    variable.
- [ ] Assert coefficients and planted association signs, not merely exit code or
  matrix existence.
- [ ] Run the job locally if the environment supports it, then verify the real
  GitHub Actions job is green.
- [ ] Commit.

### Achieved when

Every PR runs microsuite and the pinned real SpiecEasi implementation on the
same counts and compares scientific results. A green mocked test cannot mask a
wrong algorithm.

---

## Task 8: Correct documentation, changelog, and full validation

**Files:**

- Modify: `docs/methods.md`
- Modify: `CHANGELOG.md`
- Modify: `README.md` only if its network summary needs qualification
- Modify: `tests/test_container_skeletons.py` only if its documentation
  contract assertion needs updated wording

### Steps

- [ ] Change the methods table from “CLR correlation approximation” and
  “SparCC-style” to “native SparCC estimator,” naming Dirichlet normalization,
  iterative exclusion, fixed-seed reproducibility, and the absence of
  bootstrap p-values.
- [ ] Document that input must be raw nonnegative counts and explain the four
  tuning options. State that `p_value` is `NaN`; do not imply statistical
  significance.
- [ ] Add an `[Unreleased] / Fixed` changelog entry: all `sparcc` networks from
  released versions were CLR-Pearson under a wrong backend label and should be
  rerun. Do not choose a release version in this plan.
- [ ] Run focused tests with JUnit XML and report the exact count:

  ```bash
  uv run pytest tests/test_sparcc_core.py tests/test_sparcc_parity.py \
    tests/test_sparcc_reference_fixtures.py tests/test_network_method.py \
    --junitxml=/tmp/sparcc-focused.xml
  ```

- [ ] Run the complete local quality gate:

  ```bash
  uv run ruff check .
  uv run ruff format --check .
  uv run ty check
  uv run pytest --junitxml=/tmp/sparcc-full.xml
  git diff --check
  ```

- [ ] Compare any failures against the current documented platform baseline,
  but do not dismiss new failures as pre-existing without checking `main`.
- [ ] Mutation audit before final review:
  1. substitute CLR-Pearson for the estimator—parity test must fail;
  2. skip iterative exclusion—inner reference test must fail;
  3. ignore `--sparcc-seed`—CLI reproducibility test must fail;
  4. transpose the count matrix—shape/order or parity test must fail.
  Revert each mutation immediately after proving the failure.
- [ ] Run the live reference job and inspect CI status. All required jobs must
  be green before declaring completion.
- [ ] Commit documentation and final test adjustments.

### Achieved when

The implementation, public contract, documentation, changelog, default tests,
mutation checks, and live external reference all agree that `sparcc` means the
SparCC estimator. The correction is ready for review but not tagged or released.

---

## Task/commit boundaries

Each row is a review boundary. Luna should stop and resolve review findings
before advancing.

| Task | Suggested commit | Reviewer focus |
| --- | --- | --- |
| 0 | `test(sparcc): capture pinned SpiecEasi reference results` | Evidence provenance; thresholds frozen before implementation |
| 1 | `feat(sparcc): add count validation and Dirichlet normalization` | Input/RNG contract; no global randomness |
| 2 | `feat(sparcc): implement basis variance reconstruction` | Equations, `ddof=1`, symmetry, pseudoinverse |
| 3 | `feat(sparcc): implement iterative pair exclusion` | Absolute max, symmetric state, stop/tie behavior |
| 4 | `feat(sparcc): aggregate reproducible correlation estimates` | Median aggregation; covariance/correlation consistency; parity |
| 5 | `fix(network): use real SparCC estimates for sparcc backend` | No accidental Pearson recomputation; edge schema |
| 6 | `feat(cli): expose SparCC tuning and seed options` | Backward-compatible defaults; typed forwarding |
| 7 | `ci: compare SparCC against pinned SpiecEasi` | Real result assertions; exact external revision |
| 8 | `docs: document corrected SparCC backend` | Rerun warning; no significance overclaim |

## Final acceptance checklist

- [ ] Current CLR-Pearson implementation fails the new parity test.
- [ ] New deterministic inner implementation matches pinned SpiecEasi within
  `1e-9` relative tolerance.
- [ ] New outer implementation satisfies the predeclared variability-derived
  MAE bound on dense and zero-containing counts.
- [ ] New SparCC improves latent-truth MAE by at least 10% over CLR-Pearson on
  the frozen simulation and gets planted strong-edge signs right.
- [ ] Same seed gives byte-identical CLI output; different seed changes at least
  one coefficient.
- [ ] Negative/fractional/non-finite/all-zero invalid inputs fail clearly.
- [ ] Prevalence filtering, feature IDs, thresholding, sorting, `top_n`, and
  empty output schema remain correct.
- [ ] `native-correlation`, `spieceasi`, and `flashweave` focused tests are
  unchanged and green.
- [ ] Bootstrap p-values are not claimed; `p_value` is `NaN`.
- [ ] No SpiecEasi GPL source is copied or vendored.
- [ ] Ruff, formatting, type checking, full pytest, `git diff --check`, and the
  live reference CI job are green.
- [ ] No version bump, tag, push, or release occurs without user direction.
