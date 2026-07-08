# Denoise Output Validation (P2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add fail-by-default post-run output validation to `microsuite denoise`: a reusable integrity helper (missing / empty / invalid-gzip) wired through `denoise._run` for all backends, plus a `dada2-r` ASV-sample-set contract check, with a `--no-validate` escape.

**Architecture:** `runtime/validation.py` provides `validate_output_file`/`validate_outputs`. A `validate` flag threads `denoise()` → each backend function → `_run`, which calls `validate_outputs(outputs)` after a successful `run_command`. The `dada2-r` path additionally cross-checks the ASV table's sample columns against `_expected_sample_ids(input_dir)`.

**Tech Stack:** Python 3.12, `gzip`, pytest. All tests offline (no real tool execution).

## Global Constraints

- Fail-by-default: a failed check raises `MicrobiomeSuiteError`; `--no-validate` (CLI) / `validate=False` (API) disables all post-run validation.
- Empty (0-byte) output message must mention the incomplete-run / unsynced cloud-placeholder possibility (#9); invalid-gzip message names the file (#8).
- Validation is runtime-agnostic (identical for `--runtime local` and `--runtime docker`).
- `_expected_sample_ids` uses the same PE/SE grouping the R backend uses (`R[12]`, `read[12]`, `[12]`, optional `_001`; unmatched files are single-end samples named by their FASTQ stem).
- The `dada2-r` ASV table is written by R `write.table(asv_table, ..., sep="\t", col.names=NA)` → the header line is `<empty>\tsample1\tsample2\t...` (leading empty cell from `col.names=NA`); sample IDs are the header cells after the first.
- `from __future__ import annotations` at the top of every new module. Existing denoise argv tests that mock a successful `subprocess.run` without creating output files must pass `validate=False` (they test command construction, not outputs).

---

### Task 1: `runtime/validation.py` integrity helper

**Files:**
- Create: `src/microsuite/runtime/validation.py`
- Test: `tests/test_runtime_validation.py`

**Interfaces:**
- Produces:
  - `validate_output_file(path: Path, *, allow_empty: bool = False) -> None` — raises `MicrobiomeSuiteError` if missing / empty / (`.gz` and not valid gzip).
  - `validate_outputs(outputs: dict[str, str], *, allow_empty: bool = False) -> None` — validates each value path.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_runtime_validation.py
from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.runtime.validation import validate_output_file, validate_outputs


def test_missing_output_raises(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError, match="not created"):
        validate_output_file(tmp_path / "nope.tsv")


def test_empty_output_raises_with_placeholder_hint(tmp_path: Path) -> None:
    p = tmp_path / "empty.tsv"
    p.write_bytes(b"")
    with pytest.raises(MicrobiomeSuiteError, match="empty"):
        validate_output_file(p)
    # message mentions the cloud-placeholder possibility
    try:
        validate_output_file(p)
    except MicrobiomeSuiteError as exc:
        assert "placeholder" in str(exc).lower()


def test_nonempty_tsv_passes(tmp_path: Path) -> None:
    p = tmp_path / "t.tsv"
    p.write_text("a\tb\n1\t2\n", encoding="utf-8")
    validate_output_file(p)  # no raise


def test_valid_gzip_passes(tmp_path: Path) -> None:
    p = tmp_path / "reads.fastq.gz"
    p.write_bytes(gzip.compress(b"@r\nACGT\n+\nIIII\n"))
    validate_output_file(p)  # no raise


def test_invalid_gzip_raises(tmp_path: Path) -> None:
    p = tmp_path / "bad.fastq.gz"
    p.write_bytes(b"this is not gzip but has content")
    with pytest.raises(MicrobiomeSuiteError, match="gzip"):
        validate_output_file(p)


def test_allow_empty(tmp_path: Path) -> None:
    p = tmp_path / "empty.txt"
    p.write_bytes(b"")
    validate_output_file(p, allow_empty=True)  # no raise


def test_validate_outputs_raises_on_first_bad(tmp_path: Path) -> None:
    good = tmp_path / "g.tsv"
    good.write_text("x\n", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError):
        validate_outputs({"good": str(good), "missing": str(tmp_path / "no.tsv")})
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_runtime_validation.py -v`
Expected: FAIL (`ModuleNotFoundError: microsuite.runtime.validation`).

- [ ] **Step 3: Create `src/microsuite/runtime/validation.py`**

```python
from __future__ import annotations

import gzip
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError


def validate_output_file(path: Path, *, allow_empty: bool = False) -> None:
    if not path.exists():
        raise MicrobiomeSuiteError(f"Expected output was not created: {path}")
    if not allow_empty and path.stat().st_size == 0:
        raise MicrobiomeSuiteError(
            f"Output is empty: {path}. This can mean an incomplete run or an "
            "unsynced cloud-storage placeholder file."
        )
    if path.name.endswith(".gz"):
        try:
            with gzip.open(path, "rb") as handle:
                handle.read(65536)
        except (gzip.BadGzipFile, EOFError, OSError) as exc:
            raise MicrobiomeSuiteError(
                f"Output is not a valid gzip file: {path}."
            ) from exc


def validate_outputs(outputs: dict[str, str], *, allow_empty: bool = False) -> None:
    for path_str in outputs.values():
        validate_output_file(Path(path_str), allow_empty=allow_empty)
```

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest tests/test_runtime_validation.py -v`
Expected: PASS (7 passed).

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/runtime/validation.py tests/test_runtime_validation.py
git commit -m "feat(runtime): output-integrity validation (missing/empty/invalid-gzip)"
```

---

### Task 2: Thread `validate` through denoise `_run` + CLI

**Files:**
- Modify: `src/microsuite/methods/denoise.py` (`_run` + the three backend functions + `denoise()`)
- Modify: `src/microsuite/cli/method_features_cmd.py` (`--no-validate`)
- Test: `tests/test_denoise_cluster_methods.py` (new generic-validation tests + fix argv tests)

**Interfaces:**
- Consumes: `validate_outputs` (Task 1).
- Produces: `_run(..., validate: bool = True)`; `denoise_qiime2_dada2/_qiime2_deblur/_dada2_r(..., validate: bool = True)`; `denoise(..., validate: bool = True)`; CLI `--no-validate`.

- [ ] **Step 1: Write the failing test** (append to `tests/test_denoise_cluster_methods.py`)

```python
def test_denoise_validates_missing_output(tmp_path, monkeypatch) -> None:
    import subprocess
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import denoise

    input_dir = tmp_path / "reads"
    input_dir.mkdir()
    (input_dir / "s.fastq.gz").write_text("x")
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    # subprocess "succeeds" but writes no output files
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kw: subprocess.CompletedProcess(command, 0, "", ""),
    )
    with pytest.raises(MicrobiomeSuiteError, match="not created|empty"):
        denoise(
            backend="dada2-r", demux=input_dir,
            output_table=tmp_path / "table.tsv", output_rep_seqs=tmp_path / "rep.fasta",
            output_stats=tmp_path / "stats.tsv", mode="single", threads=1, force=True,
        )


def test_denoise_no_validate_skips(tmp_path, monkeypatch) -> None:
    import subprocess
    from microsuite.methods.denoise import denoise

    input_dir = tmp_path / "reads"
    input_dir.mkdir()
    (input_dir / "s.fastq.gz").write_text("x")
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kw: subprocess.CompletedProcess(command, 0, "", ""),
    )
    denoise(  # no raise despite missing outputs
        backend="dada2-r", demux=input_dir,
        output_table=tmp_path / "table.tsv", output_rep_seqs=tmp_path / "rep.fasta",
        output_stats=tmp_path / "stats.tsv", mode="single", threads=1, force=True,
        validate=False,
    )
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_denoise_cluster_methods.py -k "validates_missing or no_validate_skips" -v`
Expected: FAIL (`denoise()` has no `validate` param; missing outputs not detected).

- [ ] **Step 3: Add `validate` to `_run` and validate after success**

In `denoise.py`, import at top: `from microsuite.runtime.validation import validate_outputs`. Change `_run`'s signature to add `validate: bool = True` (keyword-only, after `params`), and after the `run_command(...)` call add:
```python
    if validate and outputs:
        validate_outputs(outputs)
```

- [ ] **Step 4: Thread `validate` through the three backend functions**

Add `validate: bool = True` (keyword-only) to `denoise_qiime2_dada2`, `denoise_qiime2_deblur`, and `denoise_dada2_r`, and pass `validate=validate` to every `_run(...)` call inside them (there are two `_run` calls in `denoise_qiime2_dada2`, one in `denoise_qiime2_deblur`, and two in `denoise_dada2_r` — the local and docker branches).

- [ ] **Step 5: Thread `validate` through `denoise()`**

Add `validate: bool = True` (keyword-only) to `denoise()`. In each backend dispatch branch (`qiime2-dada2`, `qiime2-deblur`, `dada2-r`), pass `validate=validate` to the backend function call.

- [ ] **Step 6: Wire the CLI**

In `method_features_cmd.py` `denoise_cmd`, add:
```python
        no_validate: Annotated[
            bool, typer.Option("--no-validate", help="Skip post-run output validation.")
        ] = False,
```
and add `validate=not no_validate,` to the `denoise(...)` call.

- [ ] **Step 7: Fix existing argv tests broken by default validation**

Run: `uv run pytest tests/test_denoise_cluster_methods.py -v`
Several existing "builds_command" tests mock a successful `subprocess.run` but never create output files, so default validation now raises. For EACH such failing test (they assert the constructed command, not outputs — e.g. `test_denoise_qiime2_dada2_*_builds_command`, `test_denoise_qiime2_deblur_builds_command`, `test_denoise_dada2_r_builds_rscript_command`, and the docker argv tests), add `validate=False` to its `denoise(...)` call. Do NOT add `validate=False` to tests that assert an error is raised *before* execution (missing-Rscript, ordering, rejected-param guards) — those never reach validation.

- [ ] **Step 8: Run to verify pass**

Run: `uv run pytest tests/test_denoise_cluster_methods.py -v`
Expected: PASS (new validation tests + all argv/guard tests, the argv ones now with `validate=False`).

- [ ] **Step 9: Commit**

```bash
git add src/microsuite/methods/denoise.py src/microsuite/cli/method_features_cmd.py tests/test_denoise_cluster_methods.py
git commit -m "feat(denoise): fail-by-default output validation via _run + --no-validate"
```

---

### Task 3: `dada2-r` ASV-sample contract check

**Files:**
- Modify: `src/microsuite/methods/denoise.py` (`_expected_sample_ids`, `_validate_dada2_asv_samples`; call from `denoise_dada2_r`)
- Test: `tests/test_denoise_cluster_methods.py`

**Interfaces:**
- Consumes: nothing new.
- Produces: `_expected_sample_ids(input_dir: Path) -> set[str]`; `_validate_dada2_asv_samples(output_table: Path, input_dir: Path) -> None`.

- [ ] **Step 1: Write the failing tests** (append)

```python
def test_expected_sample_ids_pe_and_se(tmp_path) -> None:
    from microsuite.methods.denoise import _expected_sample_ids
    d = tmp_path / "reads"
    d.mkdir()
    for n in ("a_1.fastq.gz", "a_2.fastq.gz", "b_R1.fastq.gz", "b_R2.fastq.gz", "c.fastq.gz"):
        (d / n).write_text("x")
    assert _expected_sample_ids(d) == {"a", "b", "c"}


def test_validate_asv_samples_ok(tmp_path) -> None:
    from microsuite.methods.denoise import _validate_dada2_asv_samples
    d = tmp_path / "reads"
    d.mkdir()
    (d / "a.fastq.gz").write_text("x")
    (d / "b.fastq.gz").write_text("x")
    table = tmp_path / "asv.tsv"
    # write.table(col.names=NA): leading empty header cell, then samples; rows are ASV ids
    table.write_text("\ta\tb\nASV1\t5\t3\n", encoding="utf-8")
    _validate_dada2_asv_samples(table, d)  # no raise


def test_validate_asv_samples_rejects_fastq_artifact(tmp_path) -> None:
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import _validate_dada2_asv_samples
    d = tmp_path / "reads"
    d.mkdir()
    (d / "a.fastq.gz").write_text("x")
    table = tmp_path / "asv.tsv"
    table.write_text("\ta.R1.filtered.fastq.gz\nASV1\t5\n", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError, match="filtered|fastq"):
        _validate_dada2_asv_samples(table, d)


def test_validate_asv_samples_rejects_mismatch(tmp_path) -> None:
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import _validate_dada2_asv_samples
    d = tmp_path / "reads"
    d.mkdir()
    (d / "a.fastq.gz").write_text("x")
    (d / "b.fastq.gz").write_text("x")
    table = tmp_path / "asv.tsv"
    table.write_text("\ta\tZ\nASV1\t5\t3\n", encoding="utf-8")  # Z not an input sample
    with pytest.raises(MicrobiomeSuiteError, match="sample"):
        _validate_dada2_asv_samples(table, d)
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_denoise_cluster_methods.py -k "expected_sample_ids or validate_asv" -v`
Expected: FAIL (helpers not defined).

- [ ] **Step 3: Add the helpers to `denoise.py`**

```python
import re

_FASTQ_EXTS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
_READ_PATTERNS = [
    re.compile(r"^(?P<sample>.+?)[._-]R(?P<read>[12])(?:[._-]001)?$"),
    re.compile(r"^(?P<sample>.+?)[._-]read(?P<read>[12])(?:[._-]001)?$", re.IGNORECASE),
    re.compile(r"^(?P<sample>.+?)[._-](?P<read>[12])(?:[._-]001)?$"),
]
_FASTQ_ARTIFACTS = (".fastq", ".fq", ".filtered")


def _fastq_stem(name: str) -> str:
    for ext in _FASTQ_EXTS:
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def _expected_sample_ids(input_dir: Path) -> set[str]:
    samples: set[str] = set()
    for path in sorted(input_dir.iterdir()):
        if not (path.is_file() and path.name.endswith(_FASTQ_EXTS)):
            continue
        stem = _fastq_stem(path.name)
        match = next((m for m in (p.match(stem) for p in _READ_PATTERNS) if m), None)
        samples.add(match.group("sample") if match else stem)
    return samples


def _validate_dada2_asv_samples(output_table: Path, input_dir: Path) -> None:
    lines = output_table.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise MicrobiomeSuiteError(f"ASV table is empty: {output_table}")
    # write.table(col.names=NA): header is "<empty>\tsample1\tsample2..."
    columns = lines[0].split("\t")[1:]
    if not columns:
        raise MicrobiomeSuiteError(f"ASV table has no sample columns: {output_table}")
    if len(set(columns)) != len(columns):
        raise MicrobiomeSuiteError(f"ASV table has duplicate sample columns: {output_table}")
    for col in columns:
        if not col or any(token in col for token in _FASTQ_ARTIFACTS):
            raise MicrobiomeSuiteError(
                f"ASV table sample column looks like a raw/filtered FASTQ name, not a "
                f"sample id: {col!r} in {output_table}"
            )
    expected = _expected_sample_ids(input_dir)
    actual = set(columns)
    if actual != expected:
        raise MicrobiomeSuiteError(
            f"ASV table sample columns do not match the input samples. "
            f"Unexpected: {sorted(actual - expected)}; missing: {sorted(expected - actual)}."
        )
```

- [ ] **Step 4: Call it from `denoise_dada2_r`**

In `denoise_dada2_r`, after the `_run(...)` returns in BOTH the docker and local branches (or once after the branch, if the function returns there), and only when `validate`, call:
```python
        if validate:
            _validate_dada2_asv_samples(output_table, input_dir)
```
Place it so it runs after a successful `_run` in each runtime path (the ASV table exists on the host for both local and docker). Keep it inside the `validate` gate so `--no-validate` skips it.

- [ ] **Step 5: Run to verify pass + full denoise suite**

Run: `uv run pytest tests/test_denoise_cluster_methods.py -v`
Expected: PASS (new semantic tests + everything from Task 2).

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/methods/denoise.py tests/test_denoise_cluster_methods.py
git commit -m "feat(dada2): validate ASV sample columns against input sample set"
```

---

## Self-Review

**Spec coverage:**
- `validate_output_file`/`validate_outputs` (missing/empty-placeholder/invalid-gzip) → Task 1. ✓
- Wired through `_run` for all denoise backends via `validate` threading → Task 2. ✓
- `--no-validate` / `validate=False` escape → Task 2 (CLI + API). ✓
- dada2-r ASV-sample contract (`_expected_sample_ids` + `_validate_dada2_asv_samples`, FASTQ-artifact + set-match) → Task 3. ✓
- Runtime-agnostic (local + docker) → Task 3 Step 4 calls after `_run` in both paths. ✓
- Existing argv tests adjusted for default validation → Task 2 Step 7. ✓

**Placeholder scan:** none — full helper code, wiring instructions with exact call-site counts, and full tests provided. Task 2 Steps 4/7 reference existing code by exact function names and the observable failure (validation raise) rather than reproducing unrelated bodies.

**Consistency:** `validate` name is used identically across `_run`, the three backend functions, `denoise()`, and the CLI (`--no-validate` → `validate=not no_validate`). `_expected_sample_ids`/`_validate_dada2_asv_samples` signatures match between Task 3's definition, tests, and the `denoise_dada2_r` call site. The ASV header format (`col.names=NA` leading empty cell) is consistent between the helper (`split("\t")[1:]`) and the test fixtures.
