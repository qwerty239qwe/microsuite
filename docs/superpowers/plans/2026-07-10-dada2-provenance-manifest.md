# DADA2 Parameter Provenance Manifest (Round-2 A) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** On every successful `dada2-r` run, write `dada2_denoise_manifest.json` beside the outputs recording the effective DADA2 params (after R-side defaults), `dada2`/R versions, and wrapper run facts (runtime, image, microsuite version, mode, threads, paths, timestamp, command).

**Architecture:** The R backend is the authoritative source of resolved params + versions: it emits a flat JSON to a new `--params-out` path. A pure Python module (`methods/dada2_manifest.py`) reads that, merges wrapper facts, and writes the final manifest. `denoise_dada2_r` wires it in after a successful run, best-effort (warns, never fails the run), identical for local and docker.

**Tech Stack:** Python 3.12, `json`/`warnings`/`datetime`, pytest; R (base only, no `jsonlite`).

## Global Constraints

- Manifest is written on **success only**, **ungated by `validate`** (provenance ≠ validation). A missing/unparseable R params file → `warnings.warn`, never raises.
- The R intermediate `dada2_r_params.json` is removed after the merge; one canonical `dada2_denoise_manifest.json` remains beside `output_stats`.
- The R side is the single source of truth for effective params: hoist the inline `value_after(flag, default)` calls into named `resolved_*` variables, reused in both the dada2 calls and the emitted JSON. Defaults verbatim from the current script: `trim_left*`/`trunc_len*` `"0"`, `max_ee*` `"2"`, `min_overlap` `"12"`, `max_merge_mismatch` `"0"`, `trunc_q` `"2"`, `n_reads_learn` `"1000000"`, `pooling` `"independent"`, `chimera` `"consensus"`, `min_fold_parent_over_abundance` `"1.0"`.
- No new user-facing CLI flag; `--params-out` is an internal wrapper→R detail.
- `run.json`/`events.jsonl` (from `run_command`) are untouched. Methods raise `MicrobiomeSuiteError` (from `microsuite._errors`) for fatal cases.
- `from __future__ import annotations` at the top of every new module.
- **JSON contract** (the flat object R writes and `read_r_params` consumes) — keys, with inactive-branch keys `null`:
  `mode` (`"paired"|"single"`), `trim_left_f`, `trim_left_r`, `trunc_len_f`, `trunc_len_r`, `max_ee_f`, `max_ee_r`, `min_overlap`, `max_merge_mismatch`, `trim_overhang` (paired-only); `trim_left`, `trunc_len`, `max_ee` (single-only); `trunc_q`, `max_n`, `rm_phix`, `pooling_method`, `chimera_method`, `min_fold_parent_over_abundance`, `allow_one_off`, `n_reads_learn` (both); `dada2_version`, `r_version`.

---

### Task 1: `methods/dada2_manifest.py` — read/build/write helpers

**Files:**
- Create: `src/microsuite/methods/dada2_manifest.py`
- Create: `tests/fixtures/dada2_r_params_paired.json` (shared contract fixture, reused by Task 2)
- Test: `tests/test_dada2_manifest.py`

**Interfaces:**
- Produces:
  - `read_r_params(path: Path) -> dict`
  - `build_manifest(r_params: dict, wrapper: dict) -> dict`
  - `write_manifest(manifest: dict, out_dir: Path) -> Path`
  - `MANIFEST_FILENAME = "dada2_denoise_manifest.json"`

- [ ] **Step 1: Write the contract fixture**

Create `tests/fixtures/dada2_r_params_paired.json` (this is exactly what the R backend emits for a paired run; Task 2 reuses it):

```json
{
  "mode": "paired",
  "trim_left_f": 0,
  "trim_left_r": 0,
  "trunc_len_f": 0,
  "trunc_len_r": 0,
  "max_ee_f": 2,
  "max_ee_r": 2,
  "min_overlap": 12,
  "max_merge_mismatch": 0,
  "trim_overhang": false,
  "trim_left": null,
  "trunc_len": null,
  "max_ee": null,
  "trunc_q": 2,
  "max_n": 0,
  "rm_phix": true,
  "pooling_method": "independent",
  "chimera_method": "consensus",
  "min_fold_parent_over_abundance": 1.0,
  "allow_one_off": false,
  "n_reads_learn": 1000000,
  "dada2_version": "1.30.0",
  "r_version": "R version 4.3.2 (2023-10-31)"
}
```

- [ ] **Step 2: Write the failing tests**

```python
# tests/test_dada2_manifest.py
from __future__ import annotations

import json
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.dada2_manifest import (
    MANIFEST_FILENAME,
    build_manifest,
    read_r_params,
    write_manifest,
)

FIXTURE = Path(__file__).parent / "fixtures" / "dada2_r_params_paired.json"


def _wrapper() -> dict:
    return {
        "microsuite_version": "9.9.9",
        "backend": "dada2-r",
        "runtime": "docker",
        "image": "ghcr.io/example/r-dada2:latest",
        "mode": "paired",
        "paired": True,
        "threads": 4,
        "input_dir": "/data/reads",
        "output_table": "/out/table.tsv",
        "output_rep_seqs": "/out/rep.fasta",
        "output_stats": "/out/stats.tsv",
        "output_plot_dir": None,
        "created_at": "2026-07-10T00:00:00+00:00",
        "command": "Rscript dada2_denoise.R --paired",
    }


def test_read_r_params_ok() -> None:
    params = read_r_params(FIXTURE)
    assert params["min_overlap"] == 12
    assert params["dada2_version"] == "1.30.0"
    assert params["trim_left"] is None  # single-only key, null for paired run


def test_read_r_params_missing(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        read_r_params(tmp_path / "nope.json")


def test_read_r_params_malformed(tmp_path: Path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text("{not json", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError):
        read_r_params(bad)


def test_build_manifest_splits_and_drops_nulls() -> None:
    manifest = build_manifest(read_r_params(FIXTURE), _wrapper())
    assert manifest["tool"] == {
        "dada2_version": "1.30.0",
        "r_version": "R version 4.3.2 (2023-10-31)",
    }
    dp = manifest["dada2_params"]
    assert dp["min_overlap"] == 12  # resolved default present, not absent
    assert dp["mode"] == "paired"
    assert "trim_left" not in dp  # single-only null key dropped
    assert "dada2_version" not in dp  # versions live under tool, not params
    assert manifest["run"]["runtime"] == "docker"
    assert manifest["run"]["image"] == "ghcr.io/example/r-dada2:latest"


def test_write_manifest_roundtrip(tmp_path: Path) -> None:
    manifest = build_manifest(read_r_params(FIXTURE), _wrapper())
    path = write_manifest(manifest, tmp_path)
    assert path.name == MANIFEST_FILENAME
    assert json.loads(path.read_text())["dada2_params"]["min_overlap"] == 12
```

- [ ] **Step 3: Run to verify it fails**

Run: `uv run pytest tests/test_dada2_manifest.py -v`
Expected: FAIL (`ModuleNotFoundError: microsuite.methods.dada2_manifest`).

- [ ] **Step 4: Create `src/microsuite/methods/dada2_manifest.py`**

```python
from __future__ import annotations

import json
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError

MANIFEST_FILENAME = "dada2_denoise_manifest.json"

_VERSION_KEYS = ("dada2_version", "r_version")


def read_r_params(path: Path) -> dict:
    if not path.exists():
        raise MicrobiomeSuiteError(f"DADA2 R params file was not written: {path}")
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        raise MicrobiomeSuiteError(f"Could not parse DADA2 R params file {path}: {exc}") from exc


def build_manifest(r_params: dict, wrapper: dict) -> dict:
    params = {
        key: value
        for key, value in r_params.items()
        if key not in _VERSION_KEYS and value is not None
    }
    return {
        "tool": {key: r_params.get(key) for key in _VERSION_KEYS},
        "dada2_params": params,
        "run": dict(wrapper),
    }


def write_manifest(manifest: dict, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / MANIFEST_FILENAME
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return path
```

- [ ] **Step 5: Run to verify pass**

Run: `uv run pytest tests/test_dada2_manifest.py -v`
Expected: PASS (all 5).

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/methods/dada2_manifest.py tests/test_dada2_manifest.py tests/fixtures/dada2_r_params_paired.json
git commit -m "feat(dada2): provenance manifest read/build/write helpers"
```

---

### Task 2: Wire the manifest into `denoise_dada2_r` (Python side, fully offline)

**Files:**
- Modify: `src/microsuite/methods/denoise.py` (`_dada2_r_script_args`, `denoise_dada2_r`, new `_emit_dada2_manifest` helper)
- Test: `tests/test_denoise_cluster_methods.py`

**Interfaces:**
- Consumes: `dada2_manifest` (Task 1); the Task-1 fixture `tests/fixtures/dada2_r_params_paired.json`.
- Produces: `_dada2_r_script_args(..., params_out: str | None = None)`; `_emit_dada2_manifest(output_stats: Path, run_facts: dict) -> None`.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_denoise_cluster_methods.py`)

```python
def _stub_run_writing_r_params(r_params_text, stats_text):
    """subprocess.run stub: writes the R params file the wrapper points --params-out
    at, plus stats/table/rep-seqs so validation passes."""
    import subprocess
    from pathlib import Path

    def fake(command, **kw):
        if "--params-out" in command:
            Path(command[command.index("--params-out") + 1]).write_text(
                r_params_text, encoding="utf-8"
            )
        if "--output-stats" in command:
            Path(command[command.index("--output-stats") + 1]).write_text(
                stats_text, encoding="utf-8"
            )
        for flag, content in (
            ("--output-table", "\tsampleP\nASV1\t5\n"),
            ("--output-rep-seqs", ">ASV1\nACGT\n"),
        ):
            if flag in command:
                Path(command[command.index(flag) + 1]).write_text(content, encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    return fake


def test_denoise_dada2_r_writes_manifest(tmp_path, monkeypatch) -> None:
    import json

    from microsuite.methods.denoise import denoise

    demux = tmp_path / "reads"
    demux.mkdir()
    (demux / "sampleP.fastq.gz").write_text("x")
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    r_params = (Path(__file__).parent / "fixtures" / "dada2_r_params_paired.json").read_text()
    healthy = (
        "\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\n"
        "sampleP\t1000\t950\t940\t930\t900\t880\n"
    )
    monkeypatch.setattr("subprocess.run", _stub_run_writing_r_params(r_params, healthy))
    denoise(
        backend="dada2-r",
        demux=demux,
        output_table=tmp_path / "table.tsv",
        output_rep_seqs=tmp_path / "rep.fasta",
        output_stats=tmp_path / "stats.tsv",
        mode="paired",
        threads=1,
        force=True,
    )
    manifest_path = tmp_path / "dada2_denoise_manifest.json"
    assert manifest_path.exists()
    manifest = json.loads(manifest_path.read_text())
    assert manifest["tool"]["dada2_version"] == "1.30.0"
    assert manifest["dada2_params"]["min_overlap"] == 12
    assert manifest["run"]["backend"] == "dada2-r"
    assert manifest["run"]["mode"] == "paired"
    # intermediate R file removed
    assert not (tmp_path / "dada2_r_params.json").exists()


def test_denoise_dada2_r_missing_params_warns_not_fatal(tmp_path, monkeypatch) -> None:
    from microsuite.methods.denoise import denoise

    demux = tmp_path / "reads"
    demux.mkdir()
    (demux / "sampleP.fastq.gz").write_text("x")
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    healthy = "\tinput\tfiltered\tnonchim\nsampleP\t1000\t950\t880\n"

    def fake_no_params(command, **kw):  # writes stats/table/rep-seqs but NO --params-out file
        import subprocess
        from pathlib import Path

        if "--output-stats" in command:
            Path(command[command.index("--output-stats") + 1]).write_text(healthy, encoding="utf-8")
        for flag, content in (
            ("--output-table", "\tsampleP\nASV1\t5\n"),
            ("--output-rep-seqs", ">ASV1\nACGT\n"),
        ):
            if flag in command:
                Path(command[command.index(flag) + 1]).write_text(content, encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_no_params)
    with pytest.warns(UserWarning, match="provenance manifest"):
        denoise(
            backend="dada2-r",
            demux=demux,
            output_table=tmp_path / "table.tsv",
            output_rep_seqs=tmp_path / "rep.fasta",
            output_stats=tmp_path / "stats.tsv",
            mode="single",
            threads=1,
            force=True,
        )
    assert not (tmp_path / "dada2_denoise_manifest.json").exists()
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_denoise_cluster_methods.py -k "manifest or missing_params" -v`
Expected: FAIL (no manifest written; `--params-out` not passed).

- [ ] **Step 3: Add `params_out` to `_dada2_r_script_args`**

In `_dada2_r_script_args` (denoise.py), add keyword-only param `params_out: str | None = None` to the signature, and after the `--output-stats` block (or anywhere before `return args`) append:

```python
    if params_out is not None:
        args.extend(["--params-out", params_out])
```

- [ ] **Step 4: Add the `_emit_dada2_manifest` helper**

Add near `_dada2_log_params` in denoise.py (module level). Ensure the imports at the top of denoise.py include: `from datetime import UTC, datetime`, `from microsuite import __version__ as _MICROSUITE_VERSION`, and `from microsuite.methods import dada2_manifest` (add whichever are missing).

```python
def _emit_dada2_manifest(output_stats: Path, run_facts: dict) -> None:
    """Merge the R-emitted resolved params with wrapper facts into
    dada2_denoise_manifest.json beside the stats file. Best-effort: a missing or
    unparseable R params file warns and never fails an otherwise-successful run."""
    r_params_path = output_stats.parent / "dada2_r_params.json"
    try:
        r_params = dada2_manifest.read_r_params(r_params_path)
        manifest = dada2_manifest.build_manifest(r_params, run_facts)
        dada2_manifest.write_manifest(manifest, output_stats.parent)
    except MicrobiomeSuiteError as exc:
        warnings.warn(f"Could not write DADA2 provenance manifest: {exc}", stacklevel=2)
    finally:
        r_params_path.unlink(missing_ok=True)
```

- [ ] **Step 5: Build the wrapper facts and pass `--params-out` in both branches**

In `denoise_dada2_r`, before the `if runtime == "docker":` branch (right after `script_res = files(...)`), build the branch-independent facts:

```python
    manifest_facts = {
        "microsuite_version": _MICROSUITE_VERSION,
        "backend": "dada2-r",
        "runtime": runtime,
        "image": resolve_dada2_image(image) if runtime == "docker" else None,
        "mode": "paired" if paired else "single",
        "paired": paired,
        "threads": threads,
        "input_dir": str(input_dir),
        "output_table": str(output_table),
        "output_rep_seqs": str(output_rep_seqs),
        "output_stats": str(output_stats),
        "output_plot_dir": str(output_plot_dir) if output_plot_dir is not None else None,
    }
```

**Docker branch:** add `params_out=mapper.to_container(output_stats.parent / "dada2_r_params.json")` to the `_dada2_r_script_args(...)` call. After the `if validate:` QC block (still inside the `with as_file` scope), add:

```python
            _emit_dada2_manifest(
                output_stats,
                {
                    **manifest_facts,
                    "created_at": datetime.now(UTC).isoformat(),
                    "command": " ".join(command),
                },
            )
```

**Local branch:** add `params_out=str(output_stats.parent / "dada2_r_params.json")` to the `_dada2_r_script_args(...)` call. After the `if validate:` QC block, add the identical `_emit_dada2_manifest(...)` call (with that branch's `command`).

Both calls sit **outside** the `if validate:` gate (provenance is ungated).

- [ ] **Step 6: Run to verify pass + no regression**

Run: `uv run pytest tests/test_denoise_cluster_methods.py -v`
Expected: PASS (new manifest tests + all existing dada2-r tests). Existing argv/validation tests pass `validate=False` or write no params file; the missing-params path only warns, so it won't break them. If an existing test now unexpectedly warns about the manifest, confirm its stub doesn't write a `--params-out` file — that is the expected best-effort warning, and the test can be left as-is unless it uses `filterwarnings("error")`; do NOT weaken `_emit_dada2_manifest`.

- [ ] **Step 7: Commit**

```bash
git add src/microsuite/methods/denoise.py tests/test_denoise_cluster_methods.py
git commit -m "feat(dada2): write provenance manifest after dada2-r run (best-effort)"
```

---

### Task 3: R backend emits resolved params + versions

**Files:**
- Modify: `src/microsuite/methods/r/dada2_denoise.R`
- Modify: `tests/integration/test_dada2_naming_contract_live.py` (assert the manifest on a real run)

**Interfaces:**
- Consumes: nothing new (invoked by Task 2's wiring via `--params-out`).
- Produces: writes the JSON contract (see Global Constraints) to the `--params-out` path on success.

- [ ] **Step 1: Hoist the inline defaults to `resolved_*` variables**

In `dada2_denoise.R`, alongside the existing top-level resolutions (after `n_reads_learn <- ...`, ~line 58), add:

```r
params_out <- value_after("--params-out")

resolved_trim_left_f <- as.integer(value_after("--trim-left-f", "0"))
resolved_trim_left_r <- as.integer(value_after("--trim-left-r", "0"))
resolved_trunc_len_f <- as.integer(value_after("--trunc-len-f", "0"))
resolved_trunc_len_r <- as.integer(value_after("--trunc-len-r", "0"))
resolved_max_ee_f <- as.numeric(value_after("--max-ee-f", "2"))
resolved_max_ee_r <- as.numeric(value_after("--max-ee-r", "2"))
resolved_min_overlap <- as.integer(value_after("--min-overlap", "12"))
resolved_max_merge_mismatch <- as.integer(value_after("--max-merge-mismatch", "0"))
resolved_trim_overhang <- has_flag("--trim-overhang")
resolved_trim_left <- as.integer(value_after("--trim-left", "0"))
resolved_trunc_len <- as.integer(value_after("--trunc-len", "0"))
resolved_max_ee <- as.numeric(value_after("--max-ee", "2"))
```

Replace the inline calls in the **paired** `filterAndTrim` (currently `trimLeft`/`truncLen`/`maxEE` using `value_after(...)`) with the `resolved_*` vars:

```r
    trimLeft = c(resolved_trim_left_f, resolved_trim_left_r),
    truncLen = c(resolved_trunc_len_f, resolved_trunc_len_r),
    maxEE = c(resolved_max_ee_f, resolved_max_ee_r),
```

Replace the **paired** `mergePairs` inline calls:

```r
    minOverlap = resolved_min_overlap,
    maxMismatch = resolved_max_merge_mismatch,
    trimOverhang = resolved_trim_overhang
```

Replace the **single** `filterAndTrim` inline calls:

```r
    trimLeft = resolved_trim_left,
    truncLen = resolved_trunc_len,
    maxEE = resolved_max_ee,
```

- [ ] **Step 2: Add the dependency-free JSON writer + emit on success**

Add these helpers near the other function definitions (e.g. after `count_reads`):

```r
json_scalar <- function(v) {
  if (is.null(v) || (length(v) == 1 && is.na(v))) return("null")
  if (is.logical(v)) return(if (isTRUE(v)) "true" else "false")
  if (is.numeric(v)) return(format(v, scientific = FALSE, trim = TRUE))
  paste0("\"", gsub("\"", "\\\\\"", as.character(v)), "\"")
}

write_params_json <- function(path, params) {
  parts <- vapply(
    names(params),
    function(k) paste0("  \"", k, "\": ", json_scalar(params[[k]])),
    character(1)
  )
  writeLines(c("{", paste(parts, collapse = ",\n"), "}"), path)
}
```

At the very end of the script (after the outputs and retention plot are written, ~line 221), append:

```r
if (!is.null(params_out)) {
  write_params_json(params_out, list(
    mode = if (paired) "paired" else "single",
    trim_left_f = if (paired) resolved_trim_left_f else NA,
    trim_left_r = if (paired) resolved_trim_left_r else NA,
    trunc_len_f = if (paired) resolved_trunc_len_f else NA,
    trunc_len_r = if (paired) resolved_trunc_len_r else NA,
    max_ee_f = if (paired) resolved_max_ee_f else NA,
    max_ee_r = if (paired) resolved_max_ee_r else NA,
    min_overlap = if (paired) resolved_min_overlap else NA,
    max_merge_mismatch = if (paired) resolved_max_merge_mismatch else NA,
    trim_overhang = if (paired) resolved_trim_overhang else NA,
    trim_left = if (!paired) resolved_trim_left else NA,
    trunc_len = if (!paired) resolved_trunc_len else NA,
    max_ee = if (!paired) resolved_max_ee else NA,
    trunc_q = trunc_q,
    max_n = max_n,
    rm_phix = rm_phix,
    pooling_method = pooling_method,
    chimera_method = chimera_method,
    min_fold_parent_over_abundance = min_fold_parent_over_abundance,
    allow_one_off = allow_one_off,
    n_reads_learn = n_reads_learn,
    dada2_version = as.character(packageVersion("dada2")),
    r_version = R.version.string
  ))
}
```

This flat object matches the Task-1 fixture `tests/fixtures/dada2_r_params_paired.json` exactly (keys, types, `null` for inactive-branch keys). Keep them in sync — the fixture is the contract.

- [ ] **Step 3: Verify the R script still parses (skips without Rscript)**

Run: `command -v Rscript && Rscript -e 'invisible(parse("src/microsuite/methods/r/dada2_denoise.R")); cat("parse OK\n")' || echo "no Rscript — skipped"`
Expected: `parse OK` where R is installed, otherwise the skip message. (Normal CI has no R; the real behavioral guard is Task 2's offline test + the live e2e below.)

- [ ] **Step 4: Extend the opt-in live e2e to assert the manifest**

In `tests/integration/test_dada2_naming_contract_live.py`, in the paired real-dada2 case (after the run succeeds and the ASV table is asserted), add:

```python
    import json

    manifest_path = output_stats.parent / "dada2_denoise_manifest.json"
    assert manifest_path.exists(), "dada2-r run must write a provenance manifest"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["tool"]["dada2_version"]  # real dada2 version, non-empty
    assert manifest["dada2_params"]["min_overlap"] == 12  # resolved default
    assert manifest["run"]["backend"] == "dada2-r"
```

(Adapt `output_stats` to the variable name already used in that test for the stats path. This case is skipped unless `MICROSUITE_RUN_EXTERNAL_INTEGRATION=1` and real dada2 is importable.)

- [ ] **Step 5: Run the offline suite (confirms no regression from the R change path)**

Run: `uv run pytest tests/test_denoise_cluster_methods.py tests/test_dada2_manifest.py -q`
Expected: PASS. (The R change is exercised for real only by the opt-in live test; offline tests use the subprocess stub.)

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/methods/r/dada2_denoise.R tests/integration/test_dada2_naming_contract_live.py
git commit -m "feat(dada2): R backend emits resolved params + versions for provenance manifest"
```

---

## Self-Review

**Spec coverage:**
- R emits resolved effective params + `dada2`/R versions via `--params-out`, single source of truth via `resolved_*` vars → Task 3. ✓
- Pure `dada2_manifest.py` (`read_r_params`/`build_manifest`/`write_manifest`) → Task 1. ✓
- Wiring in `denoise_dada2_r`, both branches, ungated by `validate`, R intermediate removed → Task 2 Steps 4-5. ✓
- Best-effort (warn, never fail) on missing/unparseable params → `_emit_dada2_manifest` (Task 2) + `test_denoise_dada2_r_missing_params_warns_not_fatal`. ✓
- Manifest beside `output_stats`; `run.json` untouched → Task 1 `MANIFEST_FILENAME` + Task 2 (no `run_command` change). ✓
- Local + docker parity (params-out in the rw-mounted output dir) → Task 2 Step 5 (`mapper.to_container(...)` vs `str(...)`). ✓
- Effective values equal what R used (e.g. `min_overlap` = 12 when unset) → contract fixture + `test_build_manifest_splits_and_drops_nulls`. ✓

**Placeholder scan:** none — full module, full tests, exact R edits with the verbatim resolved-default literals, and precise call-site placement. The one adaptation flagged (`output_stats` variable name in the live e2e) is called out explicitly.

**Consistency:** the JSON contract keys/types are identical across the Global Constraints block, the Task-1 fixture, `build_manifest`, and the Task-3 R emitter. `read_r_params`/`build_manifest`/`write_manifest`/`MANIFEST_FILENAME` names match between Task 1's module, its tests, and Task 2's `_emit_dada2_manifest`. `params_out` naming aligns across `_dada2_r_script_args`, both wiring call sites, and the R `--params-out` flag.
