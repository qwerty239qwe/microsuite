# DADA2 QC Warnings & Summaries (Round-2 B) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add DADA2 run-quality feedback to the `dada2-r` backend: a post-run retention summary (JSON+TSV) with warnings on low retention, a pre-run paired overlap check via `--amplicon-length`, and `--strict-qc` to make warnings fatal.

**Architecture:** A pure-Python `methods/dada2_qc.py` parses the denoising-stats TSV and computes retention/overlap. `denoise_dada2_r` runs the overlap check pre-flight and writes the summary + emits `warnings.warn` post-run (gated by `validate`); `--strict-qc` raises instead. No R-script changes; works for local and docker.

**Tech Stack:** Python 3.12, `gzip`/`json`/`warnings`, pytest (`pytest.warns`). All tests offline.

## Global Constraints

- Pure Python; no changes to `dada2_denoise.R`. Runs identically for `--runtime local` and `docker`.
- Stats TSV format: R `write.table(track, col.names=NA)` → header `<empty>\tinput\tfiltered\t...`; rows are samples. Paired columns: `input, filtered, denoised_f, denoised_r, merged, nonchim`; single: `input, filtered, denoised, nonchim`. Presence of a `merged` column ⇒ paired.
- Warnings via `warnings.warn` (non-fatal); `--strict-qc` (default off) turns them into `MicrobiomeSuiteError`.
- QC is gated by `validate` (so `--no-validate` skips it), consistent with P2.
- Retention thresholds (documented heuristics): `filtered/input < 0.5`, `merged/input < 0.5`, `nonchim/input < 0.4`.
- Overlap: `retained_x = trunc_len_x if trunc_len_x > 0 else read_len_x`; warn when `retained_f + retained_r - amplicon_length < min_overlap` (min_overlap = `--min-overlap` value or DADA2 default 12).
- `from __future__ import annotations` at the top of every new module.

---

### Task 1: `methods/dada2_qc.py` — retention/overlap helpers

**Files:**
- Create: `src/microsuite/methods/dada2_qc.py`
- Test: `tests/test_dada2_qc.py`

**Interfaces:**
- Produces:
  - `summarize_dada2_stats(stats_path: Path) -> dict`
  - `write_qc_summary(summary: dict, out_dir: Path) -> tuple[Path, Path]`
  - `retention_warnings(summary: dict) -> list[str]`
  - `check_overlap(*, trunc_len_f: int, trunc_len_r: int, read_len_f: int, read_len_r: int, amplicon_length: int, min_overlap: int) -> str | None`
  - `first_read_length(fastq: Path) -> int`
  - `RETENTION_THRESHOLDS: dict[str, float]`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_dada2_qc.py
from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.dada2_qc import (
    check_overlap,
    first_read_length,
    retention_warnings,
    summarize_dada2_stats,
    write_qc_summary,
)


def _paired_stats(tmp_path: Path) -> Path:
    p = tmp_path / "stats.tsv"
    p.write_text(
        "\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\n"
        "sA\t1000\t500\t480\t470\t400\t300\n"
        "sB\t2000\t1800\t1700\t1690\t1600\t1500\n",
        encoding="utf-8",
    )
    return p


def _single_stats(tmp_path: Path) -> Path:
    p = tmp_path / "stats.tsv"
    p.write_text(
        "\tinput\tfiltered\tdenoised\tnonchim\n"
        "sA\t1000\t900\t880\t850\n",
        encoding="utf-8",
    )
    return p


def test_summarize_paired(tmp_path: Path) -> None:
    s = summarize_dada2_stats(_paired_stats(tmp_path))
    assert s["paired"] is True
    assert s["overall"]["input"] == 3000
    assert s["overall"]["filtered_frac"] == pytest.approx((500 + 1800) / 3000, abs=1e-3)
    assert s["overall"]["nonchim_frac"] == pytest.approx((300 + 1500) / 3000, abs=1e-3)
    assert s["overall"]["merged_frac"] is not None
    assert s["per_sample"]["sA"]["nonchim_frac"] == pytest.approx(0.3, abs=1e-3)
    assert "->" in s["bottleneck"]


def test_summarize_single_has_no_merged(tmp_path: Path) -> None:
    s = summarize_dada2_stats(_single_stats(tmp_path))
    assert s["paired"] is False
    assert s["overall"]["merged_frac"] is None
    assert s["per_sample"]["sA"]["merged_frac"] is None


def test_write_qc_summary_roundtrip(tmp_path: Path) -> None:
    s = summarize_dada2_stats(_paired_stats(tmp_path))
    json_path, tsv_path = write_qc_summary(s, tmp_path / "qc")
    assert json_path.exists() and tsv_path.exists()
    assert json.loads(json_path.read_text())["overall"]["input"] == 3000
    assert "sample_id" in tsv_path.read_text().splitlines()[0]


def test_retention_warnings_flags_low(tmp_path: Path) -> None:
    s = summarize_dada2_stats(_paired_stats(tmp_path))  # nonchim 0.6, filtered ~0.77, merged ~0.67
    # lower it: craft a low-retention summary directly
    s["overall"] = {"input": 1000, "filtered_frac": 0.3, "merged_frac": 0.3, "nonchim_frac": 0.2}
    s["bottleneck"] = "input->filtered"
    msgs = retention_warnings(s)
    assert msgs and any("filtered/input" in m for m in msgs)
    assert any("nonchim/input" in m for m in msgs)


def test_retention_warnings_healthy_is_empty(tmp_path: Path) -> None:
    s = summarize_dada2_stats(_paired_stats(tmp_path))
    s["overall"] = {"input": 1000, "filtered_frac": 0.9, "merged_frac": 0.85, "nonchim_frac": 0.8}
    assert retention_warnings(s) == []


def test_check_overlap() -> None:
    # 150 + 150 - 250 = 50 >= 12 -> ok
    assert check_overlap(trunc_len_f=0, trunc_len_r=0, read_len_f=150, read_len_r=150,
                         amplicon_length=250, min_overlap=12) is None
    # 150 + 150 - 295 = 5 < 12 -> warn
    msg = check_overlap(trunc_len_f=0, trunc_len_r=0, read_len_f=150, read_len_r=150,
                        amplicon_length=295, min_overlap=12)
    assert msg is not None and "overlap" in msg.lower()
    # truncLen overrides read length: 120 + 120 - 250 = -10 < 12 -> warn
    assert check_overlap(trunc_len_f=120, trunc_len_r=120, read_len_f=150, read_len_r=150,
                         amplicon_length=250, min_overlap=12) is not None
    # generous truncLen: 200 + 200 - 250 = 150 >= 12 -> ok
    assert check_overlap(trunc_len_f=200, trunc_len_r=200, read_len_f=150, read_len_r=150,
                         amplicon_length=250, min_overlap=12) is None
```

```python
def test_first_read_length(tmp_path: Path) -> None:
    plain = tmp_path / "r.fastq"
    plain.write_text("@x\nACGTACGT\n+\nIIIIIIII\n", encoding="utf-8")
    assert first_read_length(plain) == 8
    gz = tmp_path / "r.fastq.gz"
    gz.write_bytes(gzip.compress(b"@x\nACGT\n+\nIIII\n"))
    assert first_read_length(gz) == 4


def test_summarize_empty_raises(tmp_path: Path) -> None:
    p = tmp_path / "empty.tsv"
    p.write_text("", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError):
        summarize_dada2_stats(p)
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_dada2_qc.py -v`
Expected: FAIL (`ModuleNotFoundError: microsuite.methods.dada2_qc`).

- [ ] **Step 3: Create `src/microsuite/methods/dada2_qc.py`**

```python
from __future__ import annotations

import gzip
import json
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError

RETENTION_THRESHOLDS: dict[str, float] = {"filtered": 0.5, "merged": 0.5, "nonchim": 0.4}


def _read_stats(stats_path: Path) -> tuple[list[str], dict[str, dict[str, int]]]:
    lines = [ln for ln in stats_path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not lines:
        raise MicrobiomeSuiteError(f"Denoising stats table is empty: {stats_path}")
    metrics = lines[0].split("\t")[1:]  # col.names=NA -> leading empty cell
    rows: dict[str, dict[str, int]] = {}
    for line in lines[1:]:
        cells = line.split("\t")
        rows[cells[0]] = {m: int(float(v)) for m, v in zip(metrics, cells[1:])}
    return metrics, rows


def _frac(vals: dict[str, int], key: str) -> float:
    inp = vals.get("input", 0)
    return round(vals.get(key, 0) / inp, 4) if inp > 0 else 0.0


def _bottleneck(totals: dict[str, int], paired: bool) -> str:
    steps = ["input", "filtered", "merged", "nonchim"] if paired else ["input", "filtered", "denoised", "nonchim"]
    worst_label, worst_ratio = steps[0] + "->" + steps[1], 2.0
    for prev, cur in zip(steps, steps[1:]):
        p = totals.get(prev, 0)
        ratio = (totals.get(cur, 0) / p) if p > 0 else 1.0
        if ratio < worst_ratio:
            worst_ratio, worst_label = ratio, f"{prev}->{cur}"
    return worst_label


def summarize_dada2_stats(stats_path: Path) -> dict:
    metrics, rows = _read_stats(stats_path)
    paired = "merged" in metrics
    per_sample: dict[str, dict] = {}
    totals: dict[str, int] = {m: 0 for m in metrics}
    for sample, vals in rows.items():
        per_sample[sample] = {
            "input": vals.get("input", 0),
            "filtered_frac": _frac(vals, "filtered"),
            "merged_frac": _frac(vals, "merged") if paired else None,
            "nonchim_frac": _frac(vals, "nonchim"),
        }
        for m in metrics:
            totals[m] += vals.get(m, 0)
    overall = {
        "input": totals.get("input", 0),
        "filtered_frac": _frac(totals, "filtered"),
        "merged_frac": _frac(totals, "merged") if paired else None,
        "nonchim_frac": _frac(totals, "nonchim"),
    }
    return {"per_sample": per_sample, "overall": overall, "bottleneck": _bottleneck(totals, paired), "paired": paired}


def write_qc_summary(summary: dict, out_dir: Path) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "dada2_qc_summary.json"
    tsv_path = out_dir / "dada2_qc_summary.tsv"
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    header = ["sample_id", "input", "filtered_frac", "merged_frac", "nonchim_frac"]
    out = ["\t".join(header)]
    for sample, entry in sorted(summary["per_sample"].items()):
        merged = "" if entry["merged_frac"] is None else entry["merged_frac"]
        out.append("\t".join(str(x) for x in [sample, entry["input"], entry["filtered_frac"], merged, entry["nonchim_frac"]]))
    tsv_path.write_text("\n".join(out) + "\n", encoding="utf-8")
    return json_path, tsv_path


def retention_warnings(summary: dict) -> list[str]:
    overall = summary["overall"]
    checks = [("filtered", "filtered/input")]
    if summary["paired"]:
        checks.append(("merged", "merged/input"))
    checks.append(("nonchim", "nonchim/input"))
    messages: list[str] = []
    for key, label in checks:
        frac = overall.get(f"{key}_frac")
        if frac is None:
            continue
        if frac < RETENTION_THRESHOLDS[key]:
            messages.append(
                f"Low DADA2 retention: {label} = {frac:.0%} (below {RETENTION_THRESHOLDS[key]:.0%}). "
                f"Bottleneck step: {summary['bottleneck']}. Check quality profiles, maxEE, truncLen, "
                "and (paired) the expected read overlap."
            )
    return messages


def check_overlap(
    *,
    trunc_len_f: int,
    trunc_len_r: int,
    read_len_f: int,
    read_len_r: int,
    amplicon_length: int,
    min_overlap: int,
) -> str | None:
    retained_f = trunc_len_f if trunc_len_f > 0 else read_len_f
    retained_r = trunc_len_r if trunc_len_r > 0 else read_len_r
    overlap = retained_f + retained_r - amplicon_length
    if overlap < min_overlap:
        return (
            f"Insufficient paired overlap: retained_f({retained_f}) + retained_r({retained_r}) "
            f"- amplicon({amplicon_length}) = {overlap} < min_overlap({min_overlap}). "
            "Merging is likely to fail; increase truncLen or verify the amplicon length."
        )
    return None


def first_read_length(fastq: Path) -> int:
    opener = gzip.open if fastq.name.endswith(".gz") else open
    with opener(fastq, "rt", encoding="utf-8") as handle:
        handle.readline()  # @header
        seq = handle.readline().strip()
    if not seq:
        raise MicrobiomeSuiteError(f"Could not read a sequence record from {fastq}")
    return len(seq)
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_dada2_qc.py -v`
Expected: PASS (all).

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/methods/dada2_qc.py tests/test_dada2_qc.py
git commit -m "feat(dada2): retention/overlap QC helpers (summary, warnings, overlap check)"
```

---

### Task 2: Wire QC into `dada2-r` + CLI

**Files:**
- Modify: `src/microsuite/methods/denoise.py` (`denoise_dada2_r`, `denoise`)
- Modify: `src/microsuite/cli/method_features_cmd.py` (`--amplicon-length`, `--strict-qc`)
- Test: `tests/test_denoise_cluster_methods.py`

**Interfaces:**
- Consumes: `dada2_qc` (Task 1).
- Produces: `denoise_dada2_r(..., amplicon_length: int | None = None, strict_qc: bool = False)`; `denoise(..., amplicon_length=None, strict_qc=False)`; CLI `--amplicon-length`/`--strict-qc`.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_denoise_cluster_methods.py`)

```python
def _stub_run_writing_stats(stats_text):
    import subprocess
    def fake(command, **kw):
        # find --output-stats and write the fixture there
        if "--output-stats" in command:
            from pathlib import Path
            Path(command[command.index("--output-stats") + 1]).write_text(stats_text, encoding="utf-8")
        # write table + rep-seqs so P2 integrity + ASV checks don't trip first
        for flag, content in (("--output-table", "\tsampleP\nASV1\t5\n"), ("--output-rep-seqs", ">ASV1\nACGT\n")):
            if flag in command:
                from pathlib import Path
                Path(command[command.index(flag) + 1]).write_text(content, encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")
    return fake


def test_denoise_dada2_r_low_retention_warns(tmp_path, monkeypatch) -> None:
    from microsuite.methods.denoise import denoise
    demux = tmp_path / "reads"
    demux.mkdir()
    (demux / "sampleP.fastq.gz").write_text("x")
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    low = "\tinput\tfiltered\tnonchim\nsampleP\t1000\t300\t200\n"
    monkeypatch.setattr("subprocess.run", _stub_run_writing_stats(low))
    with pytest.warns(UserWarning, match="Low DADA2 retention"):
        denoise(
            backend="dada2-r", demux=demux,
            output_table=tmp_path / "table.tsv", output_rep_seqs=tmp_path / "rep.fasta",
            output_stats=tmp_path / "stats.tsv", mode="single", threads=1, force=True,
        )
    assert (tmp_path / "dada2_qc_summary.json").exists()


def test_denoise_dada2_r_strict_qc_raises(tmp_path, monkeypatch) -> None:
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import denoise
    demux = tmp_path / "reads"
    demux.mkdir()
    (demux / "sampleP.fastq.gz").write_text("x")
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    low = "\tinput\tfiltered\tnonchim\nsampleP\t1000\t300\t200\n"
    monkeypatch.setattr("subprocess.run", _stub_run_writing_stats(low))
    with pytest.raises(MicrobiomeSuiteError, match="retention"):
        denoise(
            backend="dada2-r", demux=demux,
            output_table=tmp_path / "table.tsv", output_rep_seqs=tmp_path / "rep.fasta",
            output_stats=tmp_path / "stats.tsv", mode="single", threads=1, force=True,
            strict_qc=True,
        )
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_denoise_cluster_methods.py -k "low_retention or strict_qc" -v`
Expected: FAIL (`denoise()` has no `strict_qc`/QC behavior).

- [ ] **Step 3: Add QC to `denoise_dada2_r`**

Add `amplicon_length: int | None = None` and `strict_qc: bool = False` (keyword-only) to `denoise_dada2_r`'s signature. Add the import at the top of `denoise.py`:
```python
import warnings
from microsuite.methods import dada2_qc
```

**Pre-run overlap check** — place it right after the `_prepare_outputs(...)`/plot-dir block and before the `runtime`-branch execution, only for paired + `amplicon_length` set:
```python
    if amplicon_length is not None and paired:
        peek = next((p for p in sorted(input_dir.iterdir())
                     if p.is_file() and p.name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))), None)
        read_len = dada2_qc.first_read_length(peek) if peek is not None else 0
        overlap_msg = dada2_qc.check_overlap(
            trunc_len_f=tuning.trunc_len_f, trunc_len_r=tuning.trunc_len_r,
            read_len_f=read_len, read_len_r=read_len,
            amplicon_length=amplicon_length, min_overlap=tuning.min_overlap or 12,
        )
        if overlap_msg is not None:
            if strict_qc:
                raise MicrobiomeSuiteError(overlap_msg)
            warnings.warn(overlap_msg, stacklevel=2)
```

**Post-run retention summary+warn** — after the successful `_run(...)` returns, in the SAME `if validate:` region as P2's `_validate_dada2_asv_samples(...)` call (both runtime branches converge there or each has it; add it right after that ASV check):
```python
        if validate:
            _validate_dada2_asv_samples(output_table, input_dir, paired=paired)
            summary = dada2_qc.summarize_dada2_stats(output_stats)
            dada2_qc.write_qc_summary(summary, output_stats.parent)
            for message in dada2_qc.retention_warnings(summary):
                if strict_qc:
                    raise MicrobiomeSuiteError(message)
                warnings.warn(message, stacklevel=2)
```
(If the docker and local branches each call the ASV check separately, add the QC block to each; keep it inside the existing `if validate:` gate so `--no-validate` skips it.)

- [ ] **Step 4: Thread through `denoise()`**

Add `amplicon_length: int | None = None` and `strict_qc: bool = False` (keyword-only) to `denoise()`. In the `dada2-r` dispatch branch, pass `amplicon_length=amplicon_length, strict_qc=strict_qc`. Add a guard near the top: if `amplicon_length is not None` and `backend != "dada2-r"`, raise `MicrobiomeSuiteError("--amplicon-length only applies to --backend dada2-r.")`.

- [ ] **Step 5: Wire the CLI**

In `method_features_cmd.py` `denoise_cmd`, add:
```python
        amplicon_length: Annotated[
            int | None, typer.Option("--amplicon-length", help="Expected paired amplicon length (dada2-r overlap check).")
        ] = None,
        strict_qc: Annotated[
            bool, typer.Option("--strict-qc", help="Make DADA2 QC warnings (low retention / overlap) fatal.")
        ] = False,
```
and add `amplicon_length=amplicon_length, strict_qc=strict_qc,` to the `denoise(...)` call.

- [ ] **Step 6: Run to verify pass + no regression**

Run: `uv run pytest tests/test_denoise_cluster_methods.py -v`
Expected: PASS (new QC tests + existing tests). If an existing dada2-r argv/validation test now emits a `warnings.warn` that pytest escalates, it means that test's stub writes a real stats table with low retention — those existing tests either pass `validate=False` (argv tests, from P2) so QC is skipped, or don't write a stats file at all. If a test unexpectedly warns, confirm it carries `validate=False`; do NOT weaken the QC logic.

- [ ] **Step 7: Commit**

```bash
git add src/microsuite/methods/denoise.py src/microsuite/cli/method_features_cmd.py tests/test_denoise_cluster_methods.py
git commit -m "feat(dada2): retention/overlap QC warnings + --amplicon-length/--strict-qc"
```

---

## Self-Review

**Spec coverage:**
- `dada2_qc.py` (summarize/write/warnings/overlap/first_read_length) → Task 1. ✓
- Post-run retention summary + warnings, gated by `validate` → Task 2 Step 3. ✓
- Pre-run overlap check via `--amplicon-length` (paired) → Task 2 Step 3. ✓
- `--strict-qc` makes warnings fatal → Task 2 (both checks). ✓
- CLI options + `denoise()` threading + non-dada2-r guard → Task 2 Steps 4/5. ✓
- Runtime-agnostic (stats/ASV on host either way); QC in the `validate` region → Task 2 Step 3. ✓

**Placeholder scan:** none — full helper code, full tests, and exact wiring with call-site placement provided. The one flagged test correction (the `check_overlap` truncLen assertion) is called out explicitly to fix before running, not left ambiguous.

**Consistency:** `summarize_dada2_stats`/`write_qc_summary`/`retention_warnings`/`check_overlap`/`first_read_length` names and signatures match between Task 1's module, its tests, and Task 2's call sites; `amplicon_length`/`strict_qc` names align across CLI → `denoise()` → `denoise_dada2_r`; the stats header parse (`split("\t")[1:]`) matches P2 and the R `col.names=NA` format.
