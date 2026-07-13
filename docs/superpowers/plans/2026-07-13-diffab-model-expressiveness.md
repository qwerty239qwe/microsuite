# ANCOM-BC Model Expressiveness (Round-4 H) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `diffab ancombc` express real ANCOM-BC2 models — fixed/random formulas, interactions, reference levels, and the previously hard-coded controls — with model validation (columns exist, full-rank design) and resolved-config provenance.

**Architecture:** The Python wrapper writes a `params.json` (replacing 4 positional args) that the R backend consumes; the R backend relevels reference factors, validates the design (base R, before the ANCOMBC requirement), calls `ancombc2` with all resolved params, and writes the result + an `ancombc_provenance.json`.

**Tech Stack:** Python 3.12, anndata/pandas, Typer, pytest; R (ANCOMBC + jsonlite).

## Global Constraints

- Defaults are **ANCOM-BC2 native** (documented change from the old force-all-off): `prv_cut=0.10`, `lib_cut=0`, `struc_zero=False`, `neg_lb=False`, `p_adj_method="BH"`, `global/pairwise/trend/dunnet=False`, `pseudo_sens=True`, `n_cl=1`.
- `fix_formula` defaults to the `--group` value; require at least one. `--group` still sets ANCOM-BC2's `group=` (structural-zero / global-test target).
- The `params.json` contract (Python writes, R reads): keys `fix_formula` (str), `rand_formula` (str|null), `group` (str|null), `reference` (obj col→level), `prv_cut, lib_cut` (num), `struc_zero, neg_lb, global, pairwise, trend, dunnet, pseudo_sens` (bool), `p_adj_method` (str), `n_cl` (int).
- `global` is a Python keyword → the wrapper/CLI parameter is `global_test`; the `params.json` key is `"global"`.
- R command line: `ancombc.R <counts.tsv> <metadata.tsv> <params.json> <output.tsv>`.
- Fatal → `MicrobiomeSuiteError` (`microsuite._errors`). Full suite + `ty check` + `ruff check .` + `ruff format --check .` all green.
- `from __future__ import annotations` at the top of touched Python modules.

---

### Task 1: Python wrapper — `run_ancombc` params.json + validation

**Files:**
- Modify: `src/microsuite/diffab/ancombc.py`
- Test: `tests/test_diffab_ancombc.py` (new)

**Interfaces:**
- Produces: `run_ancombc(adata, *, output, group=None, fix_formula=None, rand_formula=None, reference=None, prv_cut=0.10, lib_cut=0, struc_zero=False, neg_lb=False, p_adj_method="BH", global_test=False, pairwise=False, trend=False, dunnet=False, pseudo_sens=True, n_cl=1, run_dir=None, timeout=None) -> None`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_diffab_ancombc.py
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diffab.ancombc import run_ancombc


def _adata() -> ad.AnnData:
    rng = np.random.default_rng(0)
    X = rng.integers(0, 40, size=(6, 4)).astype(float)
    obs = pd.DataFrame(
        {"phase": ["pre", "pre", "post", "post", "pre", "post"],
         "subject": ["s1", "s2", "s1", "s2", "s3", "s3"]},
        index=[f"S{i}" for i in range(6)],
    )
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"F{i}" for i in range(4)]))


def _capture(monkeypatch) -> dict:
    captured: dict = {}

    def fake_run_command(command, **kw):
        captured["command"] = command
        captured["params"] = json.loads(Path(command[4]).read_text())
        # counts/metadata written to command[2]/command[3]
        captured["counts_exists"] = Path(command[2]).exists()

    monkeypatch.setattr("microsuite.diffab.ancombc.run_command", fake_run_command)
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    return captured


def test_run_ancombc_params_defaults(tmp_path, monkeypatch) -> None:
    captured = _capture(monkeypatch)
    run_ancombc(_adata(), output=tmp_path / "out.tsv", group="phase")
    p = captured["params"]
    assert p["fix_formula"] == "phase"  # defaulted from --group
    assert p["group"] == "phase"
    assert p["rand_formula"] is None
    # ANCOM-BC2 native defaults, not force-all-off
    assert p["prv_cut"] == 0.10 and p["lib_cut"] == 0
    assert p["pseudo_sens"] is True and p["p_adj_method"] == "BH"
    assert p["global"] is False and p["n_cl"] == 1
    assert captured["command"][1].endswith("ancombc.R")
    assert captured["counts_exists"]


def test_run_ancombc_fix_formula_overrides_group_and_rand(tmp_path, monkeypatch) -> None:
    captured = _capture(monkeypatch)
    run_ancombc(
        _adata(), output=tmp_path / "out.tsv", group="phase",
        fix_formula="phase*subject", rand_formula="(1|subject)",
        prv_cut=0.2, struc_zero=True, global_test=True, n_cl=4,
        reference={"phase": "pre"},
    )
    p = captured["params"]
    assert p["fix_formula"] == "phase*subject"  # overrides group
    assert p["rand_formula"] == "(1|subject)"
    assert p["reference"] == {"phase": "pre"}
    assert p["prv_cut"] == 0.2 and p["struc_zero"] is True
    assert p["global"] is True and p["n_cl"] == 4


def test_run_ancombc_requires_a_formula(tmp_path, monkeypatch) -> None:
    _capture(monkeypatch)
    with pytest.raises(MicrobiomeSuiteError, match="fix-formula|group"):
        run_ancombc(_adata(), output=tmp_path / "out.tsv")


def test_run_ancombc_missing_reference_column_raises(tmp_path, monkeypatch) -> None:
    _capture(monkeypatch)
    with pytest.raises(MicrobiomeSuiteError, match="nope"):
        run_ancombc(_adata(), output=tmp_path / "out.tsv", group="phase",
                    reference={"nope": "x"})
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_diffab_ancombc.py -v`
Expected: FAIL (`run_ancombc` has no `fix_formula`/params.json behavior).

- [ ] **Step 3: Rewrite `run_ancombc`**

Add `import json`. Replace the function:

```python
def run_ancombc(
    adata: ad.AnnData,
    *,
    output: Path,
    group: str | None = None,
    fix_formula: str | None = None,
    rand_formula: str | None = None,
    reference: dict[str, str] | None = None,
    prv_cut: float = 0.10,
    lib_cut: int = 0,
    struc_zero: bool = False,
    neg_lb: bool = False,
    p_adj_method: str = "BH",
    global_test: bool = False,
    pairwise: bool = False,
    trend: bool = False,
    dunnet: bool = False,
    pseudo_sens: bool = True,
    n_cl: int = 1,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    resolved_fix = fix_formula or group
    if not resolved_fix:
        raise MicrobiomeSuiteError("Provide --fix-formula or --group for ANCOM-BC.")
    reference = reference or {}
    obs_cols = set(adata.obs.columns)
    referenced = ([group] if group else []) + list(reference)
    missing = [c for c in referenced if c not in obs_cols]
    if missing:
        raise MicrobiomeSuiteError(f"Metadata columns not found in obs: {missing}")

    rscript = shutil.which("Rscript")
    if rscript is None:
        raise MicrobiomeSuiteError(
            "ANCOM-BC requires external Rscript and the R packages 'ANCOMBC' and "
            "'jsonlite'. Install R and those packages, then rerun this command."
        )

    with TemporaryDirectory() as temp_dir:
        temp = Path(temp_dir)
        counts_path = temp / "counts.tsv"
        metadata_path = temp / "metadata.tsv"
        params_path = temp / "params.json"

        pd.DataFrame(
            dense_counts(adata).T, index=adata.var_names, columns=adata.obs_names
        ).to_csv(counts_path, sep="\t")
        pd.DataFrame(adata.obs).to_csv(metadata_path, sep="\t")

        params = {
            "fix_formula": resolved_fix,
            "rand_formula": rand_formula,
            "group": group,
            "reference": reference,
            "prv_cut": prv_cut,
            "lib_cut": lib_cut,
            "struc_zero": struc_zero,
            "neg_lb": neg_lb,
            "p_adj_method": p_adj_method,
            "global": global_test,
            "pairwise": pairwise,
            "trend": trend,
            "dunnet": dunnet,
            "pseudo_sens": pseudo_sens,
            "n_cl": n_cl,
        }
        params_path.write_text(json.dumps(params, indent=2), encoding="utf-8")

        run_command(
            [
                rscript,
                str(ANCOMBC_SCRIPT),
                str(counts_path),
                str(metadata_path),
                str(params_path),
                str(output),
            ],
            failure_message="ANCOM-BC failed.",
            run_dir=run_dir,
            log=CommandLog(
                task="diff_abundance",
                backend="ancombc",
                inputs={"fix_formula": resolved_fix, "rand_formula": rand_formula or ""},
                outputs={"output": str(output)},
            ),
            timeout=timeout,
        )
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_diffab_ancombc.py -v`
Expected: PASS (4).

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/diffab/ancombc.py tests/test_diffab_ancombc.py
git commit -m "feat(diffab): ancombc wrapper writes params.json (formulas + controls + reference)"
```

---

### Task 2: CLI — expose the model + controls

**Files:**
- Modify: `src/microsuite/cli/diffab_cmd.py`
- Test: `tests/test_diffab_ancombc.py` (append a CLI smoke test)

**Interfaces:**
- Consumes: `run_ancombc` (Task 1).
- Produces: `diffab ancombc` CLI with the new options.

- [ ] **Step 1: Write the failing test** (append)

```python
def test_cli_ancombc_threads_options(tmp_path, monkeypatch) -> None:
    from typer.testing import CliRunner

    from microsuite.cli.app import app
    from microsuite.io.h5ad import write_h5ad

    src = tmp_path / "t.h5ad"
    write_h5ad(_adata(), src)
    captured: dict = {}

    def fake_run_ancombc(adata, **kwargs):
        captured.update(kwargs)

    monkeypatch.setattr("microsuite.cli.diffab_cmd.run_ancombc", fake_run_ancombc)
    r = CliRunner().invoke(
        app,
        ["diffab", "ancombc", str(src), "--group", "phase",
         "--rand-formula", "(1|subject)", "--reference", "phase=pre",
         "--prv-cut", "0.2", "--global", "--n-cl", "4", "-o", str(tmp_path / "o.tsv")],
    )
    assert r.exit_code == 0, r.stdout
    assert captured["group"] == "phase"
    assert captured["rand_formula"] == "(1|subject)"
    assert captured["reference"] == {"phase": "pre"}
    assert captured["prv_cut"] == 0.2
    assert captured["global_test"] is True and captured["n_cl"] == 4
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_diffab_ancombc.py -k cli -v`
Expected: FAIL (unknown option `--rand-formula`).

- [ ] **Step 3: Rewrite the `ancombc` command**

Add a `--reference` parser helper and the options:

```python
def _parse_reference(values: list[str] | None) -> dict[str, str]:
    reference: dict[str, str] = {}
    for item in values or []:
        if "=" not in item:
            raise typer.BadParameter(f"--reference must be col=level, got: {item}")
        col, level = item.split("=", 1)
        reference[col] = level
    return reference


@app.command("ancombc")
def ancombc(
    table: Annotated[Path, typer.Argument(help="Input .h5ad table.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
    group: Annotated[str | None, typer.Option("--group", help="obs column (ANCOM-BC2 group=; also the default fix-formula).")] = None,
    fix_formula: Annotated[str | None, typer.Option("--fix-formula", help="Fixed-effects R formula RHS, e.g. 'visit*hygiene'. Default: --group.")] = None,
    rand_formula: Annotated[str | None, typer.Option("--rand-formula", help="Random-effects formula, e.g. '(1|subject_code)'.")] = None,
    reference: Annotated[list[str] | None, typer.Option("--reference", help="Factor reference level as col=level (repeatable).")] = None,
    prv_cut: Annotated[float, typer.Option("--prv-cut", help="Prevalence cutoff (ANCOM-BC2 default 0.10).")] = 0.10,
    lib_cut: Annotated[int, typer.Option("--lib-cut", help="Library-size cutoff.")] = 0,
    struc_zero: Annotated[bool, typer.Option("--struc-zero/--no-struc-zero", help="Detect structural zeros.")] = False,
    neg_lb: Annotated[bool, typer.Option("--neg-lb/--no-neg-lb", help="Classify structural zeros by lower bound.")] = False,
    p_adj_method: Annotated[str, typer.Option("--p-adj-method", help="Multiple-testing adjustment.")] = "BH",
    global_test: Annotated[bool, typer.Option("--global/--no-global", help="Global test across group levels.")] = False,
    pairwise: Annotated[bool, typer.Option("--pairwise/--no-pairwise", help="Pairwise group comparisons.")] = False,
    trend: Annotated[bool, typer.Option("--trend/--no-trend", help="Trend test across ordered levels.")] = False,
    dunnet: Annotated[bool, typer.Option("--dunnet/--no-dunnet", help="Dunnett-type comparisons to reference.")] = False,
    pseudo_sens: Annotated[bool, typer.Option("--pseudo-sens/--no-pseudo-sens", help="Pseudo-count sensitivity analysis.")] = True,
    n_cl: Annotated[int, typer.Option("--n-cl", help="Worker processes.")] = 1,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
    run_dir: Annotated[Path | None, typer.Option("--run-dir", help="Write runtime logs here.")] = None,
    timeout: Annotated[float | None, typer.Option("--timeout", help="Command timeout in seconds.")] = None,
) -> None:
    adata = read_h5ad(ensure_input(table))
    run_ancombc(
        adata,
        output=prepare_output(output, force=force),
        group=group,
        fix_formula=fix_formula,
        rand_formula=rand_formula,
        reference=_parse_reference(reference),
        prv_cut=prv_cut,
        lib_cut=lib_cut,
        struc_zero=struc_zero,
        neg_lb=neg_lb,
        p_adj_method=p_adj_method,
        global_test=global_test,
        pairwise=pairwise,
        trend=trend,
        dunnet=dunnet,
        pseudo_sens=pseudo_sens,
        n_cl=n_cl,
        run_dir=run_dir,
        timeout=timeout,
    )
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_diffab_ancombc.py -v` (all pass). Sanity: `uv run python -c "from typer.testing import CliRunner; from microsuite.cli.app import app; print(CliRunner().invoke(app,['diffab','ancombc','--help']).stdout)"` shows `--fix-formula`, `--rand-formula`, `--reference`.

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/cli/diffab_cmd.py tests/test_diffab_ancombc.py
git commit -m "feat(diffab): expose ANCOM-BC2 model + controls on the ancombc CLI"
```

---

### Task 3: R backend — params.json, relevel, validate, provenance

**Files:**
- Modify: `src/microsuite/diffab/r/ancombc.R`
- Test: `tests/test_diffab_ancombc_rscript.py` (new; skips without Rscript/jsonlite)

**Interfaces:**
- Consumes: the `params.json` contract (Task 1).
- Produces: validated `ancombc2` run + `ancombc_provenance.json` beside the output.

- [ ] **Step 1: Write the failing test** (validation runs in base R before the ANCOMBC requirement, so it is testable with just Rscript+jsonlite)

```python
# tests/test_diffab_ancombc_rscript.py
from __future__ import annotations

import json
import shutil
import subprocess
from importlib.resources import files
from pathlib import Path

import pytest

SCRIPT = str(files("microsuite.diffab.r").joinpath("ancombc.R"))


def _have_r_jsonlite() -> bool:
    if shutil.which("Rscript") is None:
        return False
    r = subprocess.run(
        ["Rscript", "-e", 'quit(status = !requireNamespace("jsonlite", quietly = TRUE))'],
        capture_output=True,
    )
    return r.returncode == 0


pytestmark = pytest.mark.skipif(
    not _have_r_jsonlite(), reason="Rscript with jsonlite not available"
)


def _write_inputs(tmp_path: Path, params: dict) -> tuple[Path, Path, Path, Path]:
    counts = tmp_path / "counts.tsv"
    counts.write_text("\ts1\ts2\ts3\ts4\nF1\t5\t1\t8\t2\nF2\t0\t3\t1\t4\n", encoding="utf-8")
    meta = tmp_path / "meta.tsv"
    # phase and time are perfectly confounded -> rank-deficient when combined
    meta.write_text(
        "\tphase\ttime\tsubject\ns1\tpre\t0\ta\ns2\tpre\t0\tb\ns3\tpost\t7\ta\ns4\tpost\t7\tb\n",
        encoding="utf-8",
    )
    pj = tmp_path / "params.json"
    pj.write_text(json.dumps(params), encoding="utf-8")
    out = tmp_path / "out.tsv"
    return counts, meta, pj, out


def test_ancombc_r_rejects_rank_deficient_design(tmp_path: Path) -> None:
    counts, meta, pj, out = _write_inputs(
        tmp_path, {"fix_formula": "phase + time", "reference": {}}
    )
    r = subprocess.run(["Rscript", SCRIPT, str(counts), str(meta), str(pj), str(out)],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "rank" in (r.stdout + r.stderr).lower()


def test_ancombc_r_reports_missing_formula_column(tmp_path: Path) -> None:
    counts, meta, pj, out = _write_inputs(tmp_path, {"fix_formula": "nope", "reference": {}})
    r = subprocess.run(["Rscript", SCRIPT, str(counts), str(meta), str(pj), str(out)],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "nope" in (r.stdout + r.stderr)
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_diffab_ancombc_rscript.py -v`
Expected: FAIL or SKIP. If Rscript+jsonlite present: FAIL (current script takes 4 positional args differently and requires ANCOMBC first, so it errors on ANCOMBC before the rank check). If absent: SKIP — that is acceptable; Task 3 is still verified by the parse check in Step 4 and the opt-in live path.

- [ ] **Step 3: Rewrite `ancombc.R`**

```r
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: ancombc.R counts.tsv metadata.tsv params.json output.tsv")
}
counts_path <- args[[1]]
metadata_path <- args[[2]]
params_path <- args[[3]]
output_path <- args[[4]]

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("R package 'jsonlite' is required to read the ANCOM-BC parameters.")
}
params <- jsonlite::fromJSON(params_path)

counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_path, row.names = 1, check.names = FALSE)

fix_formula <- params$fix_formula
rand_formula <- if (is.null(params$rand_formula)) NULL else params$rand_formula
group_col <- if (is.null(params$group)) NULL else params$group

# Apply reference levels (base R; before the ANCOMBC requirement so validation is
# reachable without the heavy package installed).
reference <- params$reference
if (length(reference)) {
  for (col in names(reference)) {
    if (!(col %in% colnames(metadata))) stop(sprintf("Reference column not found: %s", col))
    lvl <- reference[[col]]
    metadata[[col]] <- factor(metadata[[col]])
    if (!(lvl %in% levels(metadata[[col]]))) {
      stop(sprintf("Reference level '%s' not found in column '%s'.", lvl, col))
    }
    metadata[[col]] <- stats::relevel(metadata[[col]], ref = lvl)
  }
}

check_formula_columns <- function(fml) {
  if (is.null(fml) || !nzchar(fml)) return(invisible(NULL))
  vars <- all.vars(stats::as.formula(paste("~", fml)))
  missing <- setdiff(vars, colnames(metadata))
  if (length(missing)) {
    stop(sprintf("Formula references unknown metadata columns: %s", paste(missing, collapse = ", ")))
  }
}
check_formula_columns(fix_formula)
check_formula_columns(rand_formula)

# Full-rank fixed-effects design (rejects confounded models, e.g. phase + time).
mm <- stats::model.matrix(stats::as.formula(paste("~", fix_formula)), data = metadata)
rank <- qr(mm)$rank
if (rank < ncol(mm)) {
  stop(sprintf(
    "Fixed-effects design is rank deficient (rank %d < %d columns): the model is confounded. Simplify fix_formula or drop a collinear term.",
    rank, ncol(mm)
  ))
}

# Group-dependent tests need a group column.
group_tests <- c(params$global, params$pairwise, params$trend, params$dunnet)
if (any(unlist(group_tests)) && is.null(group_col)) {
  stop("--global/--pairwise/--trend/--dunnet require --group.")
}

if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
  stop("R package 'ANCOMBC' is required. Install it with BiocManager::install('ANCOMBC').")
}
if (!("ancombc2" %in% getNamespaceExports("ANCOMBC"))) {
  stop("This command requires ANCOM-BC2; update the ANCOMBC package.")
}

fit <- ANCOMBC::ancombc2(
  data = counts,
  meta_data = metadata,
  fix_formula = fix_formula,
  rand_formula = rand_formula,
  group = group_col,
  p_adj_method = params$p_adj_method,
  prv_cut = params$prv_cut,
  lib_cut = params$lib_cut,
  struc_zero = isTRUE(params$struc_zero),
  neg_lb = isTRUE(params$neg_lb),
  global = isTRUE(params$global),
  pairwise = isTRUE(params$pairwise),
  trend = isTRUE(params$trend),
  dunnet = isTRUE(params$dunnet),
  pseudo_sens = isTRUE(params$pseudo_sens),
  n_cl = params$n_cl
)
write.table(fit$res, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)

# Resolved-config provenance beside the output.
model_factors <- Filter(is.factor, metadata)
provenance <- list(
  fix_formula = fix_formula,
  rand_formula = if (is.null(rand_formula)) NA else rand_formula,
  group = if (is.null(group_col)) NA else group_col,
  reference = reference,
  controls = list(
    p_adj_method = params$p_adj_method, prv_cut = params$prv_cut, lib_cut = params$lib_cut,
    struc_zero = isTRUE(params$struc_zero), neg_lb = isTRUE(params$neg_lb),
    global = isTRUE(params$global), pairwise = isTRUE(params$pairwise),
    trend = isTRUE(params$trend), dunnet = isTRUE(params$dunnet),
    pseudo_sens = isTRUE(params$pseudo_sens), n_cl = params$n_cl
  ),
  factor_reference_levels = lapply(model_factors, function(x) levels(x)[1]),
  ancombc_version = as.character(utils::packageVersion("ANCOMBC")),
  r_version = R.version.string
)
prov_path <- file.path(dirname(output_path), "ancombc_provenance.json")
writeLines(jsonlite::toJSON(provenance, auto_unbox = TRUE, null = "null", pretty = TRUE), prov_path)
```

- [ ] **Step 4: Verify parse + run the validation tests**

Run: `command -v Rscript && Rscript -e 'invisible(parse("src/microsuite/diffab/r/ancombc.R")); cat("parse OK\n")' || echo "no Rscript — skipped"`.
Then `uv run pytest tests/test_diffab_ancombc_rscript.py -v` (passes where Rscript+jsonlite exist; skips otherwise). Note the full `ancombc2` path is exercised only by the opt-in live test (a maintainer runs it with ANCOMBC installed).

- [ ] **Step 5: Full suite + gates**

Run: `uv run pytest -q`, then `uv run ty check`, `uv run ruff check .`, `uv run ruff format --check .` (all green).

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/diffab/r/ancombc.R tests/test_diffab_ancombc_rscript.py
git commit -m "feat(diffab): ancombc.R params.json + relevel + rank validation + provenance"
```

---

## Self-Review

**Spec coverage:**
- fix/rand formulas + interactions, `--group` shorthand → Tasks 1-2. ✓
- ANCOM-BC2 controls with native defaults (not force-all-off) → Tasks 1-2 (defaults) + Task 3 (passthrough). ✓
- Column-exists + full-rank validation → Task 3 (`check_formula_columns`, `qr(mm)$rank`). ✓
- Reference levels → Task 1 (`reference` dict) + Task 2 (`--reference col=level`) + Task 3 (`relevel`). ✓
- Provenance (resolved formulas/controls/reference levels/versions) → Task 3 `ancombc_provenance.json`. ✓
- Offline plumbing tests + opt-in/skip live → Tasks 1-3 tests. ✓

**Placeholder scan:** none — full Python, full R, complete tests. The Task-3 test is explicitly gated on Rscript+jsonlite with a documented skip.

**Consistency:** the `params.json` keys are identical across the Global Constraints block, `run_ancombc` (Task 1), and `ancombc.R` (Task 3); `global_test` (Python) ↔ `"global"` (json) is stated and used consistently; the CLI option names map 1:1 to `run_ancombc` kwargs; the R command line `ancombc.R counts metadata params.json output` matches between the wrapper's `run_command` call and the R `args` parsing and the R test's invocation.
