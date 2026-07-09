# DADA2 Naming-Contract Test (P4) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Lock the dada2-r wrapper⇄backend sample-naming contract (#10): fix the lowercase-`_r1` `IGNORECASE` gap, pin `_expected_sample_ids` to a shared corpus with a CI unit test, and add an opt-in test that runs the real `dada2-r` backend and asserts the ASV columns equal the intended sample IDs.

**Architecture:** A shared corpus (`tests/naming_contract_cases.py`) drives both a CI-runnable Python unit test and an opt-in real-dada2 end-to-end test, so the two can't drift. The `IGNORECASE` fix makes lowercase read suffixes strip like the R backend.

**Tech Stack:** Python 3.12, pytest, `gzip`; the opt-in test invokes the real R DADA2 backend (skips without `Rscript`).

## Global Constraints

- The R backend is the source of truth: single mode keeps the FULL FASTQ stem; paired mode strips the read suffix (R1/R2, read1/2, `_1`/`_2`, forward/reverse, `_001`, case-insensitive).
- `_READ_PATTERNS[0]` gains `re.IGNORECASE` (and `[2]` for uniformity) so lowercase `_r1/_r2` paired matches the R backend.
- The corpus is the single source of expected IDs; both tests import it.
- The opt-in e2e is gated by `MICROSUITE_RUN_EXTERNAL_INTEGRATION=1` and skips unless `Rscript` is on PATH; it never runs in the default suite.
- `from __future__ import annotations` at the top of every new module.

---

### Task 1: IGNORECASE fix + shared corpus + CI unit test

**Files:**
- Modify: `src/microsuite/methods/denoise.py` (`_READ_PATTERNS`)
- Create: `tests/naming_contract_cases.py`
- Create: `tests/test_naming_contract.py`

**Interfaces:**
- Consumes: `microsuite.methods.denoise._expected_sample_ids(input_dir, *, paired: bool) -> set[str]`.
- Produces: `NamingCase` dataclass and `CASES: tuple[NamingCase, ...]` in `tests/naming_contract_cases.py`.

- [ ] **Step 1: Write the failing test**

```python
# tests/naming_contract_cases.py
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class NamingCase:
    label: str
    filenames: tuple[str, ...]
    paired: bool
    expected: frozenset[str]


CASES: tuple[NamingCase, ...] = (
    NamingCase("pe_R1R2", ("sampleA_R1.fastq.gz", "sampleA_R2.fastq.gz"), True, frozenset({"sampleA"})),
    NamingCase("pe_R1R2_001", ("sampleA_R1_001.fastq.gz", "sampleA_R2_001.fastq.gz"), True, frozenset({"sampleA"})),
    NamingCase("pe_lower_r1r2", ("sampleA_r1.fastq.gz", "sampleA_r2.fastq.gz"), True, frozenset({"sampleA"})),
    NamingCase("pe_read1read2", ("sampleA_read1.fastq.gz", "sampleA_read2.fastq.gz"), True, frozenset({"sampleA"})),
    NamingCase("pe_1_2", ("sampleA_1.fastq.gz", "sampleA_2.fastq.gz"), True, frozenset({"sampleA"})),
    NamingCase("pe_forward_reverse", ("sampleA_forward.fastq.gz", "sampleA_reverse.fastq.gz"), True, frozenset({"sampleA"})),
    NamingCase(
        "pe_multi",
        ("s1_R1.fastq.gz", "s1_R2.fastq.gz", "s2_R1.fastq.gz", "s2_R2.fastq.gz"),
        True,
        frozenset({"s1", "s2"}),
    ),
    NamingCase("se_plain", ("sampleA.fastq.gz",), False, frozenset({"sampleA"})),
    NamingCase("se_keeps_suffix", ("sampleA_R1.fastq.gz",), False, frozenset({"sampleA_R1"})),
)
```

```python
# tests/test_naming_contract.py
from __future__ import annotations

from pathlib import Path

import pytest

from microsuite.methods.denoise import _expected_sample_ids

from tests.naming_contract_cases import CASES, NamingCase


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.label)
def test_expected_sample_ids_matches_contract(case: NamingCase, tmp_path: Path) -> None:
    for name in case.filenames:
        (tmp_path / name).write_text("x", encoding="utf-8")
    result = _expected_sample_ids(tmp_path, paired=case.paired)
    assert result == set(case.expected)
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_naming_contract.py -v`
Expected: FAIL on the `pe_lower_r1r2` case (lowercase `_r1/_r2` not stripped → `{"sampleA_r1", "sampleA_r2"}` ≠ `{"sampleA"}`), and possibly an import-path error if `tests` isn't importable — if the `from tests.naming_contract_cases import ...` fails, change it to a direct import matching how other tests import sibling helpers (check an existing test that imports a `tests/` helper; if none, use `import sys; sys.path...`-free `from naming_contract_cases import ...` since pytest adds the test dir to `sys.path`). Prefer `from naming_contract_cases import CASES, NamingCase` if that resolves.

- [ ] **Step 3: Apply the IGNORECASE fix**

In `src/microsuite/methods/denoise.py`, change `_READ_PATTERNS` so items 0 and 2 are case-insensitive:
```python
_READ_PATTERNS = [
    re.compile(r"^(?P<sample>.+?)[._-]R(?P<read>[12])(?:[._-]001)?$", re.IGNORECASE),
    re.compile(r"^(?P<sample>.+?)[._-]read(?P<read>[12])(?:[._-]001)?$", re.IGNORECASE),
    re.compile(r"^(?P<sample>.+?)[._-](?P<read>[12])(?:[._-]001)?$", re.IGNORECASE),
    re.compile(r"^(?P<sample>.+?)[._-]?(?:forward|reverse)$", re.IGNORECASE),
]
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_naming_contract.py -v`
Expected: PASS (all `CASES`, including `pe_lower_r1r2`).

- [ ] **Step 5: Confirm no regression in the existing dada2 tests**

Run: `uv run pytest tests/test_denoise_cluster_methods.py -q`
Expected: PASS (the IGNORECASE change only broadens matching; existing ASV-sample tests still pass).

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/methods/denoise.py tests/naming_contract_cases.py tests/test_naming_contract.py
git commit -m "test(dada2): pin sample-naming contract + fix lowercase read-suffix matching"
```

---

### Task 2: Opt-in real-dada2 end-to-end contract test

**Files:**
- Create: `tests/integration/test_dada2_naming_contract_live.py`

**Interfaces:**
- Consumes: `microsuite.methods.denoise.denoise`, a real `Rscript` + `dada2` (skipped otherwise).
- Produces: an opt-in test proving the real backend's ASV columns equal the intended sample IDs.

- [ ] **Step 1: Write the test**

```python
# tests/integration/test_dada2_naming_contract_live.py
from __future__ import annotations

import gzip
import os
import random
import shutil
from pathlib import Path

import pytest

from microsuite.methods.denoise import denoise

pytestmark = pytest.mark.skipif(
    os.environ.get("MICROSUITE_RUN_EXTERNAL_INTEGRATION") != "1",
    reason="set MICROSUITE_RUN_EXTERNAL_INTEGRATION=1 to run external-tool integration tests",
)

_BASES = "ACGT"
_REVCOMP = str.maketrans("ACGT", "TGCA")


def _read(tpl: str, rng: random.Random) -> tuple[str, str]:
    seq, qual = [], []
    for pos, base in enumerate(tpl):
        q = max(20, min(40, int(38 - pos * 0.13) + rng.randint(-3, 3)))
        if rng.random() < 10 ** (-q / 10.0):
            base = rng.choice([b for b in _BASES if b != base])
        seq.append(base)
        qual.append(chr(q + 33))
    return "".join(seq), "".join(qual)


def _write_fastq_gz(path: Path, records: list[tuple[str, str]]) -> None:
    body = "".join(f"@r{i + 1}\n{seq}\n+\n{qual}\n" for i, (seq, qual) in enumerate(records))
    path.write_bytes(gzip.compress(body.encode()))


def _single_end(path: Path, seed: int, n: int = 5000, length: int = 150) -> None:
    rng = random.Random(seed)
    templates = ["".join(rng.choice(_BASES) for _ in range(length)) for _ in range(3)]
    _write_fastq_gz(path, [_read(rng.choice(templates), rng) for _ in range(n)])


def _paired_end(r1: Path, r2: Path, seed: int, n: int = 5000, amplicon: int = 250, length: int = 150) -> None:
    rng = random.Random(seed)
    templates = ["".join(rng.choice(_BASES) for _ in range(amplicon)) for _ in range(3)]
    fwd, rev = [], []
    for _ in range(n):
        tpl = rng.choice(templates)
        fwd.append(_read(tpl[:length], rng))
        rc = tpl[amplicon - length:].translate(_REVCOMP)[::-1]
        rev.append(_read(rc, rng))
    _write_fastq_gz(r1, fwd)
    _write_fastq_gz(r2, rev)


def _asv_columns(table: Path) -> set[str]:
    header = table.read_text(encoding="utf-8").splitlines()[0]
    return set(header.split("\t")[1:])  # col.names=NA -> leading empty cell


def _run(demux: Path, out: Path, *, mode: str, **kw) -> Path:
    table = out / "table.tsv"
    denoise(
        backend="dada2-r", demux=demux,
        output_table=table, output_rep_seqs=out / "rep.fasta", output_stats=out / "stats.tsv",
        mode=mode, threads=2, force=True, trunc_len=0, max_ee=2.0, **kw,
    )
    return table


def test_single_end_asv_columns_equal_sample_id(tmp_path: Path) -> None:
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is not installed on PATH")
    demux = tmp_path / "reads"
    demux.mkdir()
    _single_end(demux / "sampleS.fastq.gz", seed=7)
    table = _run(demux, tmp_path / "out", mode="single")
    assert _asv_columns(table) == {"sampleS"}


def test_paired_end_asv_columns_strip_read_suffix(tmp_path: Path) -> None:
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is not installed on PATH")
    demux = tmp_path / "reads"
    demux.mkdir()
    _paired_end(demux / "sampleP_R1.fastq.gz", demux / "sampleP_R2.fastq.gz", seed=11)
    table = _run(
        demux, tmp_path / "out", mode="paired",
        trunc_len_f=0, trunc_len_r=0, max_ee_f=2.0, max_ee_r=2.0,
    )
    assert _asv_columns(table) == {"sampleP"}
```

- [ ] **Step 2: Verify it skips by default**

Run: `uv run pytest tests/integration/test_dada2_naming_contract_live.py -q`
Expected: 2 skipped (env var unset).

- [ ] **Step 3: Verify skip-when-no-Rscript and (if available) a real run**

Run: `MICROSUITE_RUN_EXTERNAL_INTEGRATION=1 uv run pytest tests/integration/test_dada2_naming_contract_live.py -q`
Expected: if `Rscript` (with `dada2`) is NOT installed → 2 skipped ("Rscript is not installed"); if it IS installed → both pass. Record which occurred. If a real paired run yields an **empty** ASV table (learnErrors/mergePairs failed on the fixture), do NOT weaken the assertion — first increase `n`/adjust the `amplicon`/overlap so a real paired run produces `>ASV1`; only if a reliably-learnable paired fixture proves infeasible after genuine tuning, reduce the paired test to single-end coverage with an explicit `pytest.skip("paired learnable fixture infeasible: <reason>")` and note it in the report (never silently drop the assertion).

- [ ] **Step 4: Lint**

Run: `uv run ruff check tests/integration/test_dada2_naming_contract_live.py` and `uv run ruff format --check tests/integration/test_dada2_naming_contract_live.py`
Expected: clean.

- [ ] **Step 5: Commit**

```bash
git add tests/integration/test_dada2_naming_contract_live.py
git commit -m "test(dada2): opt-in end-to-end naming-contract test (real dada2-r)"
```

---

## Self-Review

**Spec coverage:**
- IGNORECASE fix (lowercase `_r1/_r2`) → Task 1 Step 3 + the `pe_lower_r1r2` corpus case. ✓
- Shared corpus (`tests/naming_contract_cases.py`) imported by both tests → Task 1 + Task 2 (the e2e asserts the same intended IDs). ✓
- CI-runnable unit test pinning `_expected_sample_ids` across conventions → Task 1. ✓
- Opt-in real-dada2 e2e (single + paired), reuses the docker.yml learnable generator, asserts ASV columns == intended IDs, skips cleanly → Task 2. ✓
- Explicit, non-silent paired fallback → Task 2 Step 3. ✓

**Placeholder scan:** none — full corpus, full tests, the exact IGNORECASE patch, and the learnable generator (ported from `.github/workflows/docker.yml`) are provided. The import-path note in Task 1 Step 2 is a resolution instruction (pick the form that imports), not a missing value.

**Consistency:** `_expected_sample_ids(input_dir, *, paired)`, `NamingCase(label, filenames, paired, expected)`, and `denoise(backend="dada2-r", mode=..., trunc_len=0, max_ee=2.0, ...)` names/signatures match the current code and across both tasks. The ASV header parse (`split("\t")[1:]`) matches P2 and the R `col.names=NA` format.
