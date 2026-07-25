# mothur Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `cluster --backend mothur` and `tax_classify --backend mothur`, plus a `microsuite workflow mothur` that composes them into the mothur MiSeq SOP from a FASTQ directory.

**Architecture:** One module, `methods/mothur.py`, owns all knowledge of how to talk to mothur: locate the binary, run a single command, parse its `Output File Names:` stdout block, and select a named output. The two backends are straight-line sequences of calls to that primitive; the workflow is an orchestrator that adds only `make.contigs`.

**Tech Stack:** Python 3.11+, typer CLI, pytest, existing `microsuite.runtime.runner.run_command`, external `mothur` 1.48.x.

## Global Constraints

- Python target `py311`; the suite runs on 3.11 and 3.12.
- ruff line-length **100**, lint rules `E, F, I, UP, B`. `uv run ruff check .`, `uv run ruff format --check .`, and `uv run ty check` must pass.
- **No new runtime dependencies.** mothur is an external binary, not a Python package.
- Every module starts with `from __future__ import annotations`.
- All failures raise `microsuite._errors.MicrobiomeSuiteError`, never bare exceptions.
- The whole test suite must pass **without mothur installed**. Every test mocks `subprocess.run` or operates on recorded text.
- Spec: `docs/superpowers/specs/2026-07-25-mothur-workflow-design.md`.

### Deviations from the spec, already agreed

1. The spec names the primitive `_run_mothur`. It is called from `cluster.py`, `tax_classify.py`, and `workflows/mothur_sop.py`, so it is public: **`run_mothur`**.
2. The spec says the parser returns outputs "keyed by extension". `make.contigs` emits **two** `.fasta` files (`.trim.contigs.fasta` and `.scrap.contigs.fasta`), so a dict silently drops one. The parser returns an **ordered list**, and `select_output()` picks by suffix, raising on 0 or >1 match.

### Amendments (agreed at pre-flight, 2026-07-25)

These override the task text below wherever they conflict.

1. **Container and fixture capture run first, as Task 0.** The original plan built the parser against hand-written fixtures and replaced them with real captures in Task 7. Docker is available, so the capture is promoted to Task 0 and the parser is written against **real mothur 1.48.2 output** from the start. This removes both the rework risk and the "tests assert against invented data" review finding. Task 7 keeps documentation only.
2. **`_reject_options` moves to `methods/_dispatch.py`.** Task 5's original text mandated copying it verbatim from `trim.py:382-389`. `methods/_dispatch.py` already exists as the shared home for cross-method dispatch helpers (it holds `require_backend`), so Task 5 lifts the function there and imports it in both `trim.py` and `tax_classify.py` instead of duplicating it.

---

### Task 0: mothur container and captured stdout fixtures

Everything downstream trusts the parser, and the parser is only as good as the output it was written against. This task produces **ground truth**: real stdout from real mothur.

**Files:**
- Create: `containers/mothur/Dockerfile`
- Modify: `.github/workflows/docker.yml` (build matrix)
- Modify: `tests/test_container_skeletons.py:14-35`
- Create: `tests/fixtures/mothur/unique_seqs.txt`
- Create: `tests/fixtures/mothur/error_exit_zero.txt`
- Create: `tests/fixtures/mothur/make_contigs.txt`

**Interfaces:**
- Produces: three fixture files containing **verbatim, unedited** mothur stdout, consumed by Task 1's parser tests.

- [ ] **Step 1: Write the failing container test**

In `tests/test_container_skeletons.py`, add to the `expected` dict (after the `microsuite-dada2` entry, line 34):

```python
        "mothur": ["mothur"],
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_container_skeletons.py -v`
Expected: FAIL — `AssertionError: mothur`

- [ ] **Step 3: Write the Dockerfile**

`containers/mothur/Dockerfile`:

```dockerfile
FROM condaforge/miniforge3:24.9.2-0

LABEL org.opencontainers.image.title="mothur"
LABEL org.opencontainers.image.description="mothur amplicon OTU clustering and classification backend"

# Expected commands: mothur
RUN conda install -y -c bioconda -c conda-forge mothur=1.48.2 && \
    conda clean -afy

ENTRYPOINT ["mothur"]
```

- [ ] **Step 4: Add to the docker build matrix**

In `.github/workflows/docker.yml`, add after the `microsuite-dada2` entry (around line 119). Copy the exact key set used by the `vsearch` entry at lines 71-74:

```yaml
          - image: mothur
            dockerfile: containers/mothur/Dockerfile
            context: .
```

- [ ] **Step 5: Run to verify the container test passes**

Run: `uv run pytest tests/test_container_skeletons.py -v`
Expected: PASS

- [ ] **Step 6: Build the image**

```bash
docker build -t microsuite/mothur:local containers/mothur
docker run --rm microsuite/mothur:local "#help()" | head -20
```

Expected: mothur's version banner, reporting 1.48.2.

- [ ] **Step 7: Capture `unique.seqs` output**

Create a scratch input and run the command, saving stdout **verbatim**:

```bash
mkdir -p scratch
printf '>a_sampleA\nACGTACGTACGT\n>b_sampleA\nACGTACGTACGT\n>c_sampleB\nTTTTACGTACGT\n' > scratch/test.fasta
docker run --rm -v "$PWD/scratch:/data" microsuite/mothur:local \
  "#set.dir(output=/data); unique.seqs(fasta=/data/test.fasta, format=count)" \
  > tests/fixtures/mothur/unique_seqs.txt
cat tests/fixtures/mothur/unique_seqs.txt
```

- [ ] **Step 8: Capture the exit-0 error case**

This is the most important fixture: it proves mothur returns 0 on failure.

```bash
docker run --rm -v "$PWD/scratch:/data" microsuite/mothur:local \
  "#set.dir(output=/data); align.seqs(fasta=/data/test.fasta, reference=/data/missing.fasta)" \
  > tests/fixtures/mothur/error_exit_zero.txt
echo "exit code: $?"
grep -c '\[ERROR\]' tests/fixtures/mothur/error_exit_zero.txt
```

Expected: exit code `0`, at least one `[ERROR]` line. **If the exit code is non-zero**, mothur 1.48.2 has changed this behaviour — record that in the report, because it weakens the case for `check_mothur_errors` and Task 1 must be told.

- [ ] **Step 9: Capture `make.contigs` output**

This fixture proves the two-`.fasta` ambiguity that `select_output` exists to catch.

```bash
printf '@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n' > scratch/sampleA_R1.fastq
printf '@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n' > scratch/sampleA_R2.fastq
printf 'sampleA\t/data/sampleA_R1.fastq\t/data/sampleA_R2.fastq\n' > scratch/stability.files
docker run --rm -v "$PWD/scratch:/data" microsuite/mothur:local \
  "#set.dir(output=/data); make.contigs(file=/data/stability.files)" \
  > tests/fixtures/mothur/make_contigs.txt
grep -A6 'Output File Names' tests/fixtures/mothur/make_contigs.txt
```

Confirm the block lists both a `.trim.contigs.fasta` and a `.scrap.contigs.fasta`. If it lists only one, note it — Task 1's ambiguity test needs a different source.

- [ ] **Step 10: Record the observed format**

Write to the report file, verbatim from the captured files:
- The exact `Output File Names:` header line, **including any trailing whitespace** (check with `grep -n 'Output File Names' -A1 tests/fixtures/mothur/*.txt | cat -A | head`).
- What terminates the block: blank line, `mothur >`, EOF, or something else.
- The exact `[ERROR]` line prefix.
- mothur's exit code on the failure case.
- Whether `make.contigs` produced two `.fasta` entries.

Task 1 is written directly from these observations, so precision here is the whole point of the task.

- [ ] **Step 11: Clean up and commit**

```bash
rm -rf scratch
uv run pytest tests/test_container_skeletons.py -v
git add containers/mothur .github/workflows/docker.yml tests/test_container_skeletons.py tests/fixtures/mothur
git commit -m "feat(mothur): add container and capture real stdout fixtures"
```

Do **not** hand-edit the captured fixtures. They are evidence; edited evidence is worthless.

---

## File Structure

| Path | Responsibility |
|---|---|
| `containers/mothur/Dockerfile` | **Create (Task 0).** bioconda mothur 1.48.2. Source of the captured fixtures. |
| `tests/fixtures/mothur/*.txt` | **Create (Task 0).** Verbatim captured mothur stdout. Evidence — never hand-edited. |
| `src/microsuite/methods/_dispatch.py` | **Modify.** Gains `reject_options`, lifted out of `trim.py`. |
| `src/microsuite/methods/mothur.py` | **Create.** Binary discovery, `run_mothur`, stdout parsing, output selection, emptiness guard. The only module that knows mothur's CLI and output conventions. |
| `src/microsuite/methods/cluster.py` | **Modify.** Add `cluster_mothur()`, extend `SUPPORTED_BACKENDS`, add the `.shared` transpose. |
| `src/microsuite/methods/tax_classify.py` | **Modify.** Add `tax_classify_mothur()`, extend `SUPPORTED_METHODS`, add option rejection. |
| `src/microsuite/cli/method_features_cmd.py` | **Modify.** mothur options on `cluster`. |
| `src/microsuite/cli/method_taxonomy_cmd.py` | **Modify.** `--taxonomy-reference` / `--taxonomy-map` on `tax_classify`. |
| `src/microsuite/workflows/mothur_sop.py` | **Create.** Stability file, `make.contigs`, orchestration. |
| `src/microsuite/workflows/catalog.py` | **Modify.** `WorkflowSpec` entry. |
| `src/microsuite/cli/workflow_cmd.py` | **Modify.** `microsuite workflow mothur` command. |
| `tests/test_mothur_parser.py` | **Create.** Parser and selector, on captured stdout. |
| `tests/test_mothur_backends.py` | **Create.** Command construction, dispatch, option rejection, transpose. |

---

### Task 1: mothur stdout parser and output selector

Pure text handling. No subprocess, no mothur. This is the load-bearing component: everything downstream trusts it to say what files a step produced.

**Files:**
- Create: `src/microsuite/methods/mothur.py`
- Create: `tests/test_mothur_parser.py`
- Read only (do not edit): `tests/fixtures/mothur/*.txt` from Task 0

**Interfaces:**
- Consumes: `microsuite._errors.MicrobiomeSuiteError`; the three captured fixtures from Task 0
- Produces:
  - `MOTHUR_ERROR_MARKER: str`
  - `parse_mothur_outputs(stdout: str) -> list[Path]`
  - `select_output(outputs: list[Path], suffix: str, *, step: str, exclude: tuple[str, ...] = ()) -> Path`
  - `check_mothur_errors(stdout: str, *, step: str) -> None`

- [ ] **Step 1: Read the captured fixtures and the Task 0 report**

The fixtures already exist — Task 0 captured them from real mothur 1.48.2. Read all three plus Task 0's recorded observations before writing anything:

```bash
cat -A tests/fixtures/mothur/unique_seqs.txt | head -30
cat tests/fixtures/mothur/error_exit_zero.txt
cat tests/fixtures/mothur/make_contigs.txt
```

`cat -A` matters: the `Output File Names:` header may carry trailing whitespace, and the parser must tolerate whatever is actually there.

**The fixtures are ground truth. Never edit them to suit the parser** — if a test fails, the parser is wrong. The test assertions below use placeholder filenames; replace them with the **actual** filenames in the captured output.

- [ ] **Step 2: Write the failing tests**

`tests/test_mothur_parser.py`:

```python
from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.mothur import (
    check_mothur_errors,
    parse_mothur_outputs,
    select_output,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "mothur"


def _fixture(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


def test_parse_returns_output_paths_in_order() -> None:
    outputs = parse_mothur_outputs(_fixture("unique_seqs.txt"))

    assert [p.name for p in outputs] == [
        "stability.trim.contigs.good.count_table",
        "stability.trim.contigs.good.unique.fasta",
    ]


def test_parse_returns_empty_when_block_absent() -> None:
    assert parse_mothur_outputs("mothur > quit()\n") == []


def test_select_output_matches_by_suffix() -> None:
    outputs = parse_mothur_outputs(_fixture("unique_seqs.txt"))

    assert select_output(outputs, ".fasta", step="unique.seqs").name.endswith(".unique.fasta")
    assert select_output(outputs, ".count_table", step="unique.seqs").name.endswith(".count_table")


def test_select_output_raises_when_missing_and_lists_what_was_produced() -> None:
    outputs = parse_mothur_outputs(_fixture("unique_seqs.txt"))

    with pytest.raises(MicrobiomeSuiteError) as excinfo:
        select_output(outputs, ".shared", step="unique.seqs")

    message = str(excinfo.value)
    assert "unique.seqs" in message
    assert ".shared" in message
    # The message must name what mothur actually produced, or the user is blind.
    assert "stability.trim.contigs.good.unique.fasta" in message


def test_select_output_raises_on_ambiguous_match() -> None:
    # make.contigs emits both .trim.contigs.fasta and .scrap.contigs.fasta.
    # Silently taking the first would hand the scrap reads downstream.
    outputs = parse_mothur_outputs(_fixture("make_contigs.txt"))

    with pytest.raises(MicrobiomeSuiteError, match="ambiguous"):
        select_output(outputs, ".fasta", step="make.contigs")


def test_select_output_exclude_resolves_ambiguity() -> None:
    outputs = parse_mothur_outputs(_fixture("make_contigs.txt"))

    chosen = select_output(outputs, ".fasta", step="make.contigs", exclude=("scrap",))

    assert chosen.name == "stability.trim.contigs.fasta"


def test_select_output_raises_when_no_outputs_at_all() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="produced no output files"):
        select_output([], ".fasta", step="screen.seqs")


def test_check_errors_raises_on_error_marker() -> None:
    # mothur exits 0 here. Trusting the return code would sail straight past it.
    with pytest.raises(MicrobiomeSuiteError, match="does not exist"):
        check_mothur_errors(_fixture("error_exit_zero.txt"), step="align.seqs")


def test_check_errors_passes_clean_output() -> None:
    check_mothur_errors(_fixture("unique_seqs.txt"), step="unique.seqs")


def test_check_errors_ignores_warnings() -> None:
    check_mothur_errors("[WARNING]: blank sequence removed\n", step="screen.seqs")
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `uv run pytest tests/test_mothur_parser.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'microsuite.methods.mothur'`

- [ ] **Step 4: Write the implementation**

`src/microsuite/methods/mothur.py`:

```python
"""Wrapper around the external ``mothur`` command.

mothur derives its own output filenames by appending a tag per command, so a
caller cannot predict them. Every command instead prints an ``Output File
Names:`` block; this module reads that block rather than re-deriving mothur's
naming convention, which changes between releases.
"""

from __future__ import annotations

from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError

MOTHUR_ERROR_MARKER = "[ERROR]"

_OUTPUT_HEADER = "Output File Names:"


def parse_mothur_outputs(stdout: str) -> list[Path]:
    """Return the paths listed under mothur's final ``Output File Names:`` block."""
    lines = stdout.splitlines()
    start = None
    for index, line in enumerate(lines):
        if line.strip().startswith(_OUTPUT_HEADER):
            start = index + 1
    if start is None:
        return []

    outputs: list[Path] = []
    for line in lines[start:]:
        candidate = line.strip()
        if not candidate:
            break
        if candidate.startswith("mothur >") or candidate.startswith("["):
            break
        outputs.append(Path(candidate))
    return outputs


def select_output(
    outputs: list[Path],
    suffix: str,
    *,
    step: str,
    exclude: tuple[str, ...] = (),
) -> Path:
    """Return the single output whose filename ends with ``suffix``.

    Raises when nothing matches or when more than one does. Both are ambiguity,
    and guessing produces a plausible-looking wrong result rather than a crash.
    """
    if not outputs:
        raise MicrobiomeSuiteError(f"mothur step '{step}' produced no output files.")

    matches = [
        path
        for path in outputs
        if path.name.endswith(suffix) and not any(token in path.name for token in exclude)
    ]
    produced = ", ".join(path.name for path in outputs)
    if not matches:
        raise MicrobiomeSuiteError(
            f"mothur step '{step}' produced no '{suffix}' output. Produced: {produced}"
        )
    if len(matches) > 1:
        ambiguous = ", ".join(path.name for path in matches)
        raise MicrobiomeSuiteError(
            f"mothur step '{step}' produced an ambiguous '{suffix}' match: {ambiguous}"
        )
    return matches[0]


def check_mothur_errors(stdout: str, *, step: str) -> None:
    """Raise if mothur reported an error.

    mothur returns exit code 0 even when a command fails, so the return code
    cannot be trusted. Without this scan a failed step is silently skipped and
    the next step consumes a stale file from an earlier run, producing a
    well-formed but wrong result.
    """
    errors = [
        line.strip()
        for line in stdout.splitlines()
        if line.strip().startswith(MOTHUR_ERROR_MARKER)
    ]
    if errors:
        detail = " ".join(errors)
        raise MicrobiomeSuiteError(f"mothur step '{step}' failed: {detail}")
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `uv run pytest tests/test_mothur_parser.py -v`
Expected: PASS — 10 passed

- [ ] **Step 6: Lint and format**

Run: `uv run ruff check src/microsuite/methods/mothur.py tests/test_mothur_parser.py && uv run ruff format src/microsuite/methods/mothur.py tests/test_mothur_parser.py && uv run ty check`
Expected: no errors

- [ ] **Step 7: Commit**

```bash
git add src/microsuite/methods/mothur.py tests/test_mothur_parser.py
git commit -m "feat(mothur): parse Output File Names and detect exit-0 errors"
```

> **Fixture provenance:** the three samples are verbatim stdout captured from mothur 1.48.2 in Task 0. They are evidence, not illustration. If a parser test fails, fix the parser — never the fixture.

---

### Task 2: The `run_mothur` primitive

**Files:**
- Modify: `src/microsuite/methods/mothur.py`
- Create: `tests/test_mothur_backends.py`

**Interfaces:**
- Consumes: `parse_mothur_outputs`, `select_output`, `check_mothur_errors` (Task 1); `microsuite.runtime.runner.run_command`, `CommandLog`
- Produces:
  - `find_mothur() -> str`
  - `format_mothur_command(command: str, params: dict[str, str]) -> str`
  - `run_mothur(command: str, params: dict[str, str], *, work_dir: Path, run_dir: Path | None = None, timeout: float | None = None) -> list[Path]`
  - `ensure_non_empty_fasta(path: Path, *, step: str) -> Path`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_mothur_backends.py` (create the file with this content):

```python
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.mothur import (
    ensure_non_empty_fasta,
    find_mothur,
    format_mothur_command,
    run_mothur,
)

CLEAN_STDOUT = """mothur > unique.seqs(fasta=in.fasta, format=count)

Output File Names: 
out.count_table
out.unique.fasta

"""


def _fake_which(name: str) -> str | None:
    return "/usr/bin/mothur" if name == "mothur" else None


def test_format_mothur_command_sets_output_dir_first() -> None:
    text = format_mothur_command(
        "unique.seqs", {"fasta": "in.fasta", "format": "count"}
    )

    assert text == "#unique.seqs(fasta=in.fasta, format=count)"


def test_find_mothur_raises_with_install_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("shutil.which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match="container"):
        find_mothur()


def test_run_mothur_builds_command_with_set_dir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr("shutil.which", _fake_which)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, CLEAN_STDOUT, "")

    monkeypatch.setattr("subprocess.run", fake_run)

    outputs = run_mothur(
        "unique.seqs", {"fasta": "in.fasta"}, work_dir=tmp_path
    )

    assert calls == [
        [
            "/usr/bin/mothur",
            f"#set.dir(output={tmp_path}); unique.seqs(fasta=in.fasta)",
        ]
    ]
    assert [p.name for p in outputs] == ["out.count_table", "out.unique.fasta"]


def test_run_mothur_raises_on_exit_zero_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("shutil.which", _fake_which)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(command, 0, "[ERROR]: it broke.\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    with pytest.raises(MicrobiomeSuiteError, match="it broke"):
        run_mothur("align.seqs", {"fasta": "in.fasta"}, work_dir=tmp_path)


def test_ensure_non_empty_fasta_raises_on_empty_file(tmp_path: Path) -> None:
    empty = tmp_path / "empty.fasta"
    empty.write_text("", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="screen.seqs"):
        ensure_non_empty_fasta(empty, step="screen.seqs")


def test_ensure_non_empty_fasta_accepts_a_record(tmp_path: Path) -> None:
    populated = tmp_path / "ok.fasta"
    populated.write_text(">seq1\nACGT\n", encoding="utf-8")

    assert ensure_non_empty_fasta(populated, step="screen.seqs") == populated
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest tests/test_mothur_backends.py -v`
Expected: FAIL — `ImportError: cannot import name 'find_mothur'`

- [ ] **Step 3: Write the implementation**

Add to `src/microsuite/methods/mothur.py` — new imports at the top, then the functions:

```python
import shutil

from microsuite.runtime.runner import CommandLog, run_command
```

```python
def find_mothur() -> str:
    """Return the mothur executable path, or raise with an install hint."""
    mothur = shutil.which("mothur")
    if mothur is None:
        raise MicrobiomeSuiteError(
            "mothur requires the external 'mothur' command. "
            "Install mothur or use the microsuite/mothur container and rerun this command."
        )
    return mothur


def format_mothur_command(command: str, params: dict[str, str]) -> str:
    """Render a single mothur command in mothur's '#command(k=v, ...)' syntax."""
    rendered = ", ".join(f"{key}={value}" for key, value in params.items())
    return f"#{command}({rendered})"


def run_mothur(
    command: str,
    params: dict[str, str],
    *,
    work_dir: Path,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> list[Path]:
    """Run one mothur command and return the files it reports producing.

    ``set.dir(output=...)`` is not cosmetic: without it mothur writes its
    intermediates and its ``mothur.<timestamp>.logfile`` into the process
    working directory.
    """
    mothur = find_mothur()
    work_dir.mkdir(parents=True, exist_ok=True)
    script = f"#set.dir(output={work_dir}); {format_mothur_command(command, params)[1:]}"

    result = run_command(
        [mothur, script],
        f"mothur step '{command}' failed.",
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(task="mothur", backend="mothur", params={"step": command, **params}),
    )
    check_mothur_errors(result.stdout or "", step=command)
    return parse_mothur_outputs(result.stdout or "")


def ensure_non_empty_fasta(path: Path, *, step: str) -> Path:
    """Raise if a FASTA has no records.

    mothur writes an empty FASTA and continues when a filter removes every
    sequence, so the failure would otherwise surface as an empty feature table
    many steps later.
    """
    if not path.exists() or not any(
        line.startswith(">") for line in path.read_text(encoding="utf-8").splitlines()
    ):
        raise MicrobiomeSuiteError(
            f"mothur step '{step}' removed every sequence. "
            "Relax the screening parameters and rerun."
        )
    return path
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest tests/test_mothur_backends.py -v`
Expected: PASS — 6 passed

- [ ] **Step 5: Lint, format, full suite**

Run: `uv run ruff check . && uv run ruff format --check . && uv run ty check && uv run pytest -q`
Expected: no new failures beyond the 11 known Windows-platform failures (`test_metadata_models`, `test_metadata_redact`, `test_run_fastp_multiqc_script`, `test_runtime_container`, `test_system_doctor`). On Linux CI all should pass.

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/methods/mothur.py tests/test_mothur_backends.py
git commit -m "feat(mothur): add run_mothur primitive with error and emptiness guards"
```

---

### Task 3: `.shared` to feature-major TSV transpose

mothur's `.shared` is sample-major; microsuite's table contract is feature-major. Pure function, no mothur.

**Files:**
- Modify: `src/microsuite/methods/cluster.py`
- Modify: `tests/test_mothur_backends.py`

**Interfaces:**
- Produces: `write_otu_table_from_shared(shared: Path, output: Path) -> None` in `microsuite.methods.cluster`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_mothur_backends.py`:

```python
from microsuite.methods.cluster import write_otu_table_from_shared


def test_write_otu_table_from_shared_transposes_to_feature_major(tmp_path: Path) -> None:
    shared = tmp_path / "final.opti_mcc.shared"
    shared.write_text(
        "label\tGroup\tnumOtus\tOtu0001\tOtu0002\n"
        "0.03\tsampleA\t2\t5\t3\n"
        "0.03\tsampleB\t2\t0\t7\n",
        encoding="utf-8",
    )
    output = tmp_path / "table.tsv"

    write_otu_table_from_shared(shared, output)

    assert output.read_text(encoding="utf-8") == (
        "feature-id\tsampleA\tsampleB\n"
        "Otu0001\t5\t0\n"
        "Otu0002\t3\t7\n"
    )


def test_write_otu_table_from_shared_keeps_all_zero_samples(tmp_path: Path) -> None:
    # A sample that survived filtering but shares no OTUs must stay as a column,
    # or downstream sample counts silently disagree with the metadata.
    shared = tmp_path / "final.opti_mcc.shared"
    shared.write_text(
        "label\tGroup\tnumOtus\tOtu0001\n"
        "0.03\tsampleA\t1\t9\n"
        "0.03\tsampleB\t1\t0\n",
        encoding="utf-8",
    )
    output = tmp_path / "table.tsv"

    write_otu_table_from_shared(shared, output)

    assert output.read_text(encoding="utf-8").splitlines()[0] == "feature-id\tsampleA\tsampleB"


def test_write_otu_table_from_shared_rejects_empty_file(tmp_path: Path) -> None:
    shared = tmp_path / "empty.shared"
    shared.write_text("", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="no rows"):
        write_otu_table_from_shared(shared, tmp_path / "table.tsv")
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_mothur_backends.py -k shared -v`
Expected: FAIL — `ImportError: cannot import name 'write_otu_table_from_shared'`

- [ ] **Step 3: Implement**

Append to `src/microsuite/methods/cluster.py`:

```python
def write_otu_table_from_shared(shared: Path, output: Path) -> None:
    """Convert mothur's sample-major .shared into microsuite's feature-major TSV.

    .shared columns are: label, Group, numOtus, then one column per OTU.
    """
    rows = [
        line.rstrip("\n").split("\t")
        for line in shared.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(rows) < 2:
        raise MicrobiomeSuiteError(f"mothur .shared file has no rows: {shared}")

    header, *records = rows
    otus = header[3:]
    samples = [record[1] for record in records]

    with output.open("w", encoding="utf-8", newline="") as handle:
        handle.write("feature-id\t" + "\t".join(samples) + "\n")
        for index, otu in enumerate(otus):
            counts = [record[3 + index] for record in records]
            handle.write(otu + "\t" + "\t".join(counts) + "\n")
```

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest tests/test_mothur_backends.py -k shared -v`
Expected: PASS — 3 passed

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/methods/cluster.py tests/test_mothur_backends.py
git commit -m "feat(mothur): transpose .shared into feature-major table TSV"
```

---

### Task 4: `cluster --backend mothur`

**Files:**
- Modify: `src/microsuite/methods/cluster.py`
- Modify: `src/microsuite/cli/method_features_cmd.py:334-398`
- Modify: `tests/test_mothur_backends.py`

**Interfaces:**
- Consumes: `run_mothur`, `select_output`, `ensure_non_empty_fasta` (Tasks 1–2); `write_otu_table_from_shared` (Task 3)
- Produces: `cluster_mothur(*, rep_seqs, output_table, output_rep_seqs, reference_alignment, identity, maxambig, maxhomop, pre_cluster_diffs, force, run_dir, timeout) -> None`; `SUPPORTED_BACKENDS` gains `"mothur"`

**Pipeline.** Eleven `run_mothur` calls, straight-line:

| # | mothur command | Key params | Output taken |
|---|---|---|---|
| 1 | `unique.seqs` | `fasta`, `format=count` | `.fasta`, `.count_table` |
| 2 | `align.seqs` | `fasta`, `reference` | `.align` |
| 3 | `screen.seqs` | `fasta`, `count`, `optimize=start-end`, `criteria=90`, `maxambig`, `maxhomop` | `.fasta`, `.count_table` |
| 4 | `filter.seqs` | `fasta`, `vertical=T`, `trump=.` | `.fasta` |
| 5 | `unique.seqs` | `fasta`, `count` | `.fasta`, `.count_table` |
| 6 | `pre.cluster` | `fasta`, `count`, `diffs` | `.fasta`, `.count_table` |
| 7 | `chimera.vsearch` | `fasta`, `count`, `dereplicate=t` | `.fasta`, `.count_table` |
| 8 | `dist.seqs` | `fasta`, `cutoff` | `.dist` |
| 9 | `cluster` | `column`, `count`, `method=opti`, `cutoff` | `.list` |
| 10 | `make.shared` | `list`, `count`, `label` | `.shared` |
| 11 | `get.oturep` | `column`, `list`, `count`, `fasta`, `label`, `method=abundance` | `.fasta` |

`optimize=start-end, criteria=90` is used rather than literal `start=`/`end=` coordinates because those depend on which reference alignment the user supplied; hardcoding them would silently discard everything against a different reference.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_mothur_backends.py`:

```python
from microsuite.methods.cluster import SUPPORTED_BACKENDS, cluster


def _mothur_stdout(*names: str) -> str:
    listed = "\n".join(names)
    return f"mothur > step\n\nOutput File Names: \n{listed}\n\n"


def _sop_stdouts(tmp_path: Path) -> list[str]:
    """One canned stdout per SOP step, in order."""
    base = str(tmp_path / "seqs")
    return [
        _mothur_stdout(f"{base}.unique.fasta", f"{base}.count_table"),
        _mothur_stdout(f"{base}.unique.align"),
        _mothur_stdout(f"{base}.good.fasta", f"{base}.good.count_table"),
        _mothur_stdout(f"{base}.filter.fasta"),
        _mothur_stdout(f"{base}.filter.unique.fasta", f"{base}.filter.count_table"),
        _mothur_stdout(f"{base}.precluster.fasta", f"{base}.precluster.count_table"),
        _mothur_stdout(f"{base}.pick.fasta", f"{base}.pick.count_table"),
        _mothur_stdout(f"{base}.dist"),
        _mothur_stdout(f"{base}.opti_mcc.list"),
        _mothur_stdout(f"{base}.opti_mcc.shared"),
        _mothur_stdout(f"{base}.rep.fasta"),
    ]


def test_mothur_is_a_supported_cluster_backend() -> None:
    assert "mothur" in SUPPORTED_BACKENDS


def test_cluster_mothur_runs_the_sop_in_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a_1\nACGT\n", encoding="utf-8")
    reference = tmp_path / "silva.align"
    reference.write_text(">ref\nAC-GT\n", encoding="utf-8")
    monkeypatch.setattr("shutil.which", _fake_which)

    scripts: list[str] = []
    stdouts = iter(_sop_stdouts(tmp_path))

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        scripts.append(command[1])
        return subprocess.CompletedProcess(command, 0, next(stdouts), "")

    monkeypatch.setattr("subprocess.run", fake_run)

    # get.oturep and make.shared outputs are read back, so create them.
    (tmp_path / "seqs.opti_mcc.shared").write_text(
        "label\tGroup\tnumOtus\tOtu0001\n0.03\tsampleA\t1\t4\n", encoding="utf-8"
    )
    (tmp_path / "seqs.rep.fasta").write_text(">Otu0001\nACGT\n", encoding="utf-8")
    for name in ("seqs.unique.fasta", "seqs.good.fasta", "seqs.pick.fasta"):
        (tmp_path / name).write_text(">a_1\nACGT\n", encoding="utf-8")

    cluster(
        backend="mothur",
        rep_seqs=seqs,
        output_table=tmp_path / "table.tsv",
        output_rep_seqs=tmp_path / "rep.fasta",
        reference_alignment=reference,
        identity=0.97,
    )

    invoked = [script.split("; ", 1)[1].split("(", 1)[0] for script in scripts]
    assert invoked == [
        "unique.seqs",
        "align.seqs",
        "screen.seqs",
        "filter.seqs",
        "unique.seqs",
        "pre.cluster",
        "chimera.vsearch",
        "dist.seqs",
        "cluster",
        "make.shared",
        "get.oturep",
    ]


def test_cluster_mothur_converts_identity_to_distance_cutoff(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a_1\nACGT\n", encoding="utf-8")
    reference = tmp_path / "silva.align"
    reference.write_text(">ref\nAC-GT\n", encoding="utf-8")
    monkeypatch.setattr("shutil.which", _fake_which)

    scripts: list[str] = []
    stdouts = iter(_sop_stdouts(tmp_path))

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        scripts.append(command[1])
        return subprocess.CompletedProcess(command, 0, next(stdouts), "")

    monkeypatch.setattr("subprocess.run", fake_run)
    (tmp_path / "seqs.opti_mcc.shared").write_text(
        "label\tGroup\tnumOtus\tOtu0001\n0.03\tsampleA\t1\t4\n", encoding="utf-8"
    )
    (tmp_path / "seqs.rep.fasta").write_text(">Otu0001\nACGT\n", encoding="utf-8")
    for name in ("seqs.unique.fasta", "seqs.good.fasta", "seqs.pick.fasta"):
        (tmp_path / name).write_text(">a_1\nACGT\n", encoding="utf-8")

    cluster(
        backend="mothur",
        rep_seqs=seqs,
        output_table=tmp_path / "table.tsv",
        output_rep_seqs=tmp_path / "rep.fasta",
        reference_alignment=reference,
        identity=0.97,
    )

    dist_script = next(s for s in scripts if "dist.seqs(" in s)
    assert "cutoff=0.03" in dist_script


def test_cluster_mothur_requires_reference_alignment(tmp_path: Path) -> None:
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a_1\nACGT\n", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="--reference-alignment"):
        cluster(
            backend="mothur",
            rep_seqs=seqs,
            output_table=tmp_path / "table.tsv",
            output_rep_seqs=tmp_path / "rep.fasta",
            identity=0.97,
        )


def test_cluster_mothur_validates_reference_before_running_mothur(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # A bad reference path must fail immediately, not six steps in.
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a_1\nACGT\n", encoding="utf-8")
    monkeypatch.setattr("shutil.which", _fake_which)

    def explode(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        raise AssertionError("mothur must not run when the reference is missing")

    monkeypatch.setattr("subprocess.run", explode)

    with pytest.raises(MicrobiomeSuiteError, match="does not exist"):
        cluster(
            backend="mothur",
            rep_seqs=seqs,
            output_table=tmp_path / "table.tsv",
            output_rep_seqs=tmp_path / "rep.fasta",
            reference_alignment=tmp_path / "absent.align",
            identity=0.97,
        )
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_mothur_backends.py -k cluster_mothur -v`
Expected: FAIL — `TypeError: cluster() got an unexpected keyword argument 'reference_alignment'`

- [ ] **Step 3: Extend the dispatcher**

In `src/microsuite/methods/cluster.py`, change line 11 and the `cluster()` signature and body:

```python
SUPPORTED_BACKENDS = ("vsearch", "usearch", "qiime2-vsearch", "mothur")
```

Add these parameters to `cluster()` after `sample_field: int = 0,`:

```python
    reference_alignment: Path | None = None,
    maxambig: int = 0,
    maxhomop: int = 8,
    pre_cluster_diffs: int = 2,
```

And add this branch immediately after `backend = require_backend(...)`:

```python
    if backend == "mothur":
        if reference_alignment is None:
            raise MicrobiomeSuiteError(
                "--reference-alignment is required for --backend mothur."
            )
        cluster_mothur(
            rep_seqs=rep_seqs,
            output_table=output_table,
            output_rep_seqs=output_rep_seqs,
            reference_alignment=reference_alignment,
            identity=identity,
            maxambig=maxambig,
            maxhomop=maxhomop,
            pre_cluster_diffs=pre_cluster_diffs,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
```

- [ ] **Step 4: Implement `cluster_mothur`**

Add to `src/microsuite/methods/cluster.py`, with the import at the top:

```python
from microsuite.methods.mothur import ensure_non_empty_fasta, run_mothur, select_output
```

```python
def cluster_mothur(
    *,
    rep_seqs: Path,
    output_table: Path,
    output_rep_seqs: Path,
    reference_alignment: Path,
    identity: float,
    maxambig: int,
    maxhomop: int,
    pre_cluster_diffs: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    """Run the clustering half of mothur's MiSeq SOP over a FASTA."""
    if not 0 < identity <= 1:
        raise MicrobiomeSuiteError("--identity must be greater than 0 and less than or equal to 1.")

    ensure_input(rep_seqs)
    ensure_input(reference_alignment)
    prepare_output(output_table, force=force)
    prepare_output(output_rep_seqs, force=force)
    shared_sidecar = output_table.with_suffix(".shared")
    prepare_output(shared_sidecar, force=force)

    cutoff = round(1.0 - identity, 4)
    work_dir = output_table.parent / f"{output_table.stem}.mothur"
    work_dir.mkdir(parents=True, exist_ok=True)

    def step(name: str, params: dict[str, str]) -> list[Path]:
        step_run_dir = None if run_dir is None else run_dir / name.replace(".", "_")
        return run_mothur(
            name, params, work_dir=work_dir, run_dir=step_run_dir, timeout=timeout
        )

    out = step("unique.seqs", {"fasta": str(rep_seqs), "format": "count"})
    fasta = ensure_non_empty_fasta(
        select_output(out, ".fasta", step="unique.seqs"), step="unique.seqs"
    )
    count = select_output(out, ".count_table", step="unique.seqs")

    out = step("align.seqs", {"fasta": str(fasta), "reference": str(reference_alignment)})
    aligned = select_output(out, ".align", step="align.seqs")

    out = step(
        "screen.seqs",
        {
            "fasta": str(aligned),
            "count": str(count),
            "optimize": "start-end",
            "criteria": "90",
            "maxambig": str(maxambig),
            "maxhomop": str(maxhomop),
        },
    )
    fasta = ensure_non_empty_fasta(
        select_output(out, ".fasta", step="screen.seqs"), step="screen.seqs"
    )
    count = select_output(out, ".count_table", step="screen.seqs")

    out = step("filter.seqs", {"fasta": str(fasta), "vertical": "T", "trump": "."})
    fasta = select_output(out, ".fasta", step="filter.seqs")

    out = step("unique.seqs", {"fasta": str(fasta), "count": str(count)})
    fasta = select_output(out, ".fasta", step="unique.seqs")
    count = select_output(out, ".count_table", step="unique.seqs")

    out = step(
        "pre.cluster",
        {"fasta": str(fasta), "count": str(count), "diffs": str(pre_cluster_diffs)},
    )
    fasta = select_output(out, ".fasta", step="pre.cluster")
    count = select_output(out, ".count_table", step="pre.cluster")

    out = step(
        "chimera.vsearch",
        {"fasta": str(fasta), "count": str(count), "dereplicate": "t"},
    )
    fasta = ensure_non_empty_fasta(
        select_output(out, ".fasta", step="chimera.vsearch"), step="chimera.vsearch"
    )
    count = select_output(out, ".count_table", step="chimera.vsearch")

    out = step("dist.seqs", {"fasta": str(fasta), "cutoff": str(cutoff)})
    column = select_output(out, ".dist", step="dist.seqs")

    out = step(
        "cluster",
        {
            "column": str(column),
            "count": str(count),
            "method": "opti",
            "cutoff": str(cutoff),
        },
    )
    otu_list = select_output(out, ".list", step="cluster")

    out = step(
        "make.shared",
        {"list": str(otu_list), "count": str(count), "label": str(cutoff)},
    )
    shared = select_output(out, ".shared", step="make.shared")

    out = step(
        "get.oturep",
        {
            "column": str(column),
            "list": str(otu_list),
            "count": str(count),
            "fasta": str(fasta),
            "label": str(cutoff),
            "method": "abundance",
        },
    )
    representatives = select_output(out, ".fasta", step="get.oturep")

    write_otu_table_from_shared(shared, output_table)
    shutil.copyfile(shared, shared_sidecar)
    shutil.copyfile(representatives, output_rep_seqs)
```

- [ ] **Step 5: Run to verify it passes**

Run: `uv run pytest tests/test_mothur_backends.py -v`
Expected: PASS — all tests

- [ ] **Step 6: Wire the CLI**

In `src/microsuite/cli/method_features_cmd.py`, add these options to `cluster_cmd` after `sample_field` (line 375):

```python
        reference_alignment: Annotated[
            Path | None,
            typer.Option(
                "--reference-alignment",
                help="Aligned reference FASTA (mothur backend). Required for --backend mothur.",
            ),
        ] = None,
        maxambig: Annotated[
            int,
            typer.Option("--maxambig", min=0, help="Max ambiguous bases (mothur screen.seqs)."),
        ] = 0,
        maxhomop: Annotated[
            int,
            typer.Option("--maxhomop", min=1, help="Max homopolymer length (mothur screen.seqs)."),
        ] = 8,
        pre_cluster_diffs: Annotated[
            int,
            typer.Option("--pre-cluster-diffs", min=0, help="mothur pre.cluster diffs."),
        ] = 2,
```

And pass them through in the `cluster(...)` call:

```python
            reference_alignment=reference_alignment,
            maxambig=maxambig,
            maxhomop=maxhomop,
            pre_cluster_diffs=pre_cluster_diffs,
```

- [ ] **Step 7: Verify the CLI help**

Run: `uv run microsuite cluster --help`
Expected: `--reference-alignment`, `--maxambig`, `--maxhomop`, `--pre-cluster-diffs` all listed.

- [ ] **Step 8: Lint, format, full suite, commit**

```bash
uv run ruff check . && uv run ruff format --check . && uv run ty check && uv run pytest -q
git add src/microsuite/methods/cluster.py src/microsuite/cli/method_features_cmd.py tests/test_mothur_backends.py
git commit -m "feat(cluster): add mothur OTU backend running the MiSeq SOP"
```

---

### Task 5: `tax_classify --backend mothur`

**Files:**
- Modify: `src/microsuite/methods/tax_classify.py:12-101`
- Modify: `src/microsuite/cli/method_taxonomy_cmd.py:28-95`
- Modify: `tests/test_mothur_backends.py`

**Interfaces:**
- Consumes: `run_mothur`, `select_output` (Tasks 1–2)
- Produces: `tax_classify_mothur(*, rep_seqs, output, taxonomy_reference, taxonomy_map, otu_list, count_table, force, run_dir, timeout) -> None`; `SUPPORTED_METHODS` gains `"mothur"`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_mothur_backends.py`:

```python
from microsuite.methods.tax_classify import SUPPORTED_METHODS, tax_classify


def test_mothur_is_a_supported_taxonomy_backend() -> None:
    assert "mothur" in SUPPORTED_METHODS


def test_tax_classify_mothur_builds_classify_seqs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a\nACGT\n", encoding="utf-8")
    ref = tmp_path / "trainset.fasta"
    ref.write_text(">r\nACGT\n", encoding="utf-8")
    tax = tmp_path / "trainset.tax"
    tax.write_text("r\tBacteria;Firmicutes;\n", encoding="utf-8")
    monkeypatch.setattr("shutil.which", _fake_which)

    scripts: list[str] = []
    produced = tmp_path / "seqs.wang.taxonomy"
    produced.write_text("a\tBacteria(100);\n", encoding="utf-8")

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        scripts.append(command[1])
        return subprocess.CompletedProcess(command, 0, _mothur_stdout(str(produced)), "")

    monkeypatch.setattr("subprocess.run", fake_run)

    tax_classify(
        backend="mothur",
        rep_seqs=seqs,
        output=tmp_path / "taxonomy.tsv",
        taxonomy_reference=ref,
        taxonomy_map=tax,
    )

    assert "classify.seqs(" in scripts[0]
    assert f"reference={ref}" in scripts[0]
    assert f"taxonomy={tax}" in scripts[0]


def test_tax_classify_mothur_rejects_classifier(tmp_path: Path) -> None:
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a\nACGT\n", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="--classifier"):
        tax_classify(
            backend="mothur",
            rep_seqs=seqs,
            output=tmp_path / "taxonomy.tsv",
            classifier=tmp_path / "classifier.qza",
            taxonomy_reference=tmp_path / "ref.fasta",
            taxonomy_map=tmp_path / "ref.tax",
        )


def test_tax_classify_mothur_requires_both_reference_files(tmp_path: Path) -> None:
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a\nACGT\n", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="--taxonomy-map"):
        tax_classify(
            backend="mothur",
            rep_seqs=seqs,
            output=tmp_path / "taxonomy.tsv",
            taxonomy_reference=tmp_path / "ref.fasta",
        )


def test_non_mothur_backend_rejects_taxonomy_reference(tmp_path: Path) -> None:
    # Silently ignoring this would classify against the wrong database and
    # return a well-formed, wrong taxonomy table.
    seqs = tmp_path / "seqs.qza"
    seqs.write_text("", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="--taxonomy-reference"):
        tax_classify(
            backend="kraken2",
            rep_seqs=seqs,
            output=tmp_path / "report.txt",
            classifier=tmp_path / "db",
            taxonomy_reference=tmp_path / "ref.fasta",
        )
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_mothur_backends.py -k tax_classify -v`
Expected: FAIL — `TypeError: tax_classify() got an unexpected keyword argument 'taxonomy_reference'`

- [ ] **Step 3: Extend the dispatcher**

In `src/microsuite/methods/tax_classify.py`, change line 12:

```python
SUPPORTED_METHODS = ("qiime2", "kraken2", "bracken", "metaphlan", "emu", "mothur")
```

Add to the `tax_classify()` signature after `classifier`:

```python
    taxonomy_reference: Path | None = None,
    taxonomy_map: Path | None = None,
    otu_list: Path | None = None,
    count_table: Path | None = None,
```

Insert this immediately after `backend = require_backend(...)` on line 34, **before** `resolve_classifier` is applied:

```python
    if backend == "mothur":
        _reject_options(
            "mothur",
            {"--classifier": classifier, "--input-type": None},
        )
        if taxonomy_reference is None or taxonomy_map is None:
            raise MicrobiomeSuiteError(
                "--taxonomy-reference and --taxonomy-map are both required "
                "for --backend mothur."
            )
    else:
        _reject_options(
            backend,
            {
                "--taxonomy-reference": taxonomy_reference,
                "--taxonomy-map": taxonomy_map,
                "--otu-list": otu_list,
                "--count-table": count_table,
            },
        )
```

Move the `resolve_classifier` block (lines 29-32) to **after** this check so a rejected `--classifier` fails before any refdb lookup.

**Amendment 2 applies here.** Do not copy `_reject_options`. Move it to the shared dispatch module.

Cut the function from `src/microsuite/methods/trim.py:382-389` and add it to `src/microsuite/methods/_dispatch.py`, renamed to public since it now crosses modules:

```python
def reject_options(backend: str, options: dict[str, object | None]) -> None:
    """Raise if any option that the chosen backend does not support was supplied.

    Silently ignoring an unsupported option is worse than failing: the command
    succeeds and returns a result computed without it.
    """
    rejected = [
        option
        for option, value in options.items()
        if value is not None and value is not False and value != []
    ]
    if rejected:
        raise MicrobiomeSuiteError(f"{', '.join(rejected)} not supported by --backend {backend}.")
```

In `trim.py`, delete the local definition and import it, keeping the existing call sites working by aliasing at the import:

```python
from microsuite.methods._dispatch import reject_options as _reject_options
```

In `tax_classify.py`, import it directly:

```python
from microsuite.methods._dispatch import reject_options, require_backend
```

Use `reject_options(...)` (not `_reject_options`) in the new `tax_classify` code below. Run `uv run pytest tests/test_methods.py -v` after the move to confirm the existing trim tests still pass.

Add the dispatch branch inside the `stage_execution` block, alongside the others:

```python
        if backend == "mothur":
            tax_classify_mothur(
                rep_seqs=rep_seqs,
                output=output,
                taxonomy_reference=taxonomy_reference,
                taxonomy_map=taxonomy_map,
                otu_list=otu_list,
                count_table=count_table,
                force=force,
                run_dir=run_dir,
                timeout=timeout,
            )
            return
```

- [ ] **Step 4: Implement `tax_classify_mothur`**

Add to `src/microsuite/methods/tax_classify.py`, importing at the top:

```python
from microsuite.methods.mothur import run_mothur, select_output
```

```python
def tax_classify_mothur(
    *,
    rep_seqs: Path,
    output: Path,
    taxonomy_reference: Path,
    taxonomy_map: Path,
    otu_list: Path | None,
    count_table: Path | None,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    """Classify sequences with mothur's naive Bayes classifier.

    With an ``otu_list`` from the mothur cluster backend, the per-sequence
    assignments are consensus-collapsed per OTU via classify.otu.
    """
    ensure_input(rep_seqs)
    ensure_input(taxonomy_reference)
    ensure_input(taxonomy_map)
    prepare_output(output, force=force)

    work_dir = output.parent / f"{output.stem}.mothur"
    work_dir.mkdir(parents=True, exist_ok=True)

    def step(name: str, params: dict[str, str]) -> list[Path]:
        step_run_dir = None if run_dir is None else run_dir / name.replace(".", "_")
        return run_mothur(
            name, params, work_dir=work_dir, run_dir=step_run_dir, timeout=timeout
        )

    classify_params = {
        "fasta": str(rep_seqs),
        "reference": str(taxonomy_reference),
        "taxonomy": str(taxonomy_map),
    }
    if count_table is not None:
        ensure_input(count_table)
        classify_params["count"] = str(count_table)

    out = step("classify.seqs", classify_params)
    assignments = select_output(out, ".taxonomy", step="classify.seqs")

    if otu_list is not None:
        ensure_input(otu_list)
        otu_params = {"list": str(otu_list), "taxonomy": str(assignments)}
        if count_table is not None:
            otu_params["count"] = str(count_table)
        out = step("classify.otu", otu_params)
        assignments = select_output(out, ".cons.taxonomy", step="classify.otu")

    shutil.copyfile(assignments, output)
```

- [ ] **Step 5: Run to verify it passes**

Run: `uv run pytest tests/test_mothur_backends.py -v`
Expected: PASS — all tests

- [ ] **Step 6: Wire the CLI**

In `src/microsuite/cli/method_taxonomy_cmd.py`, add after the `classifier` option (line 52):

```python
        taxonomy_reference: Annotated[
            Path | None,
            typer.Option(
                "--taxonomy-reference",
                help="Reference sequence FASTA. mothur backend only.",
            ),
        ] = None,
        taxonomy_map: Annotated[
            Path | None,
            typer.Option(
                "--taxonomy-map",
                help="Reference taxonomy file, 'id<TAB>lineage'. mothur backend only.",
            ),
        ] = None,
        otu_list: Annotated[
            Path | None,
            typer.Option(
                "--otu-list",
                help="mothur .list file; enables per-OTU consensus taxonomy.",
            ),
        ] = None,
        count_table: Annotated[
            Path | None,
            typer.Option("--count-table", help="mothur .count_table file."),
        ] = None,
```

And in the `tax_classify(...)` call:

```python
            taxonomy_reference=taxonomy_reference,
            taxonomy_map=taxonomy_map,
            otu_list=otu_list,
            count_table=count_table,
```

- [ ] **Step 7: Verify and commit**

```bash
uv run microsuite tax_classify --help
uv run ruff check . && uv run ruff format --check . && uv run ty check && uv run pytest -q
git add src/microsuite/methods/tax_classify.py src/microsuite/cli/method_taxonomy_cmd.py tests/test_mothur_backends.py
git commit -m "feat(tax_classify): add mothur classifier backend with paired reference files"
```

---

### Task 6: `microsuite workflow mothur`

**Files:**
- Create: `src/microsuite/workflows/mothur_sop.py`
- Modify: `src/microsuite/workflows/catalog.py:15-39`
- Modify: `src/microsuite/cli/workflow_cmd.py`
- Create: `tests/test_mothur_workflow.py`

**Interfaces:**
- Consumes: `run_mothur`, `select_output` (Tasks 1–2); `cluster` (Task 4); `tax_classify` (Task 5)
- Produces: `write_stability_file(reads_dir: Path, output: Path) -> Path`; `run_mothur_sop(*, reads_dir, output_dir, reference_alignment, taxonomy_reference, taxonomy_map, identity=0.97, force=False, timeout=None) -> None`

- [ ] **Step 1: Write the failing tests**

`tests/test_mothur_workflow.py`:

```python
from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.workflows.catalog import WORKFLOWS
from microsuite.workflows.mothur_sop import write_stability_file


def test_mothur_workflow_is_in_the_catalog() -> None:
    assert any(workflow.name == "mothur" for workflow in WORKFLOWS)


def test_write_stability_file_pairs_r1_and_r2(tmp_path: Path) -> None:
    reads = tmp_path / "reads"
    reads.mkdir()
    for name in ("sampleA_R1.fastq.gz", "sampleA_R2.fastq.gz",
                 "sampleB_R1.fastq.gz", "sampleB_R2.fastq.gz"):
        (reads / name).write_text("", encoding="utf-8")

    stability = write_stability_file(reads, tmp_path / "stability.files")

    lines = [line.split("\t") for line in stability.read_text(encoding="utf-8").splitlines()]
    assert [line[0] for line in lines] == ["sampleA", "sampleB"]
    assert lines[0][1].endswith("sampleA_R1.fastq.gz")
    assert lines[0][2].endswith("sampleA_R2.fastq.gz")


def test_write_stability_file_rejects_unpaired_reads(tmp_path: Path) -> None:
    # A dropped mate silently halves a sample's depth if it is not caught here.
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "sampleA_R1.fastq.gz").write_text("", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="sampleA"):
        write_stability_file(reads, tmp_path / "stability.files")


def test_write_stability_file_rejects_empty_directory(tmp_path: Path) -> None:
    reads = tmp_path / "reads"
    reads.mkdir()

    with pytest.raises(MicrobiomeSuiteError, match="No paired FASTQ"):
        write_stability_file(reads, tmp_path / "stability.files")
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_mothur_workflow.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'microsuite.workflows.mothur_sop'`

- [ ] **Step 3: Implement**

`src/microsuite/workflows/mothur_sop.py`:

```python
"""End-to-end mothur MiSeq SOP: FASTQ directory to OTU table and taxonomy."""

from __future__ import annotations

import re
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.cluster import cluster
from microsuite.methods.mothur import run_mothur, select_output
from microsuite.methods.tax_classify import tax_classify

_MATE = re.compile(r"^(?P<sample>.+?)_(?:R)?(?P<mate>[12])(?:_001)?\.f(?:ast)?q(?:\.gz)?$")


def write_stability_file(reads_dir: Path, output: Path) -> Path:
    """Write mothur's stability file: sample name, R1 path, R2 path per line."""
    pairs: dict[str, dict[str, Path]] = {}
    for path in sorted(reads_dir.iterdir()):
        match = _MATE.match(path.name)
        if match is None:
            continue
        pairs.setdefault(match.group("sample"), {})[match.group("mate")] = path

    if not pairs:
        raise MicrobiomeSuiteError(f"No paired FASTQ files found in {reads_dir}.")

    unpaired = sorted(sample for sample, mates in pairs.items() if len(mates) != 2)
    if unpaired:
        raise MicrobiomeSuiteError(
            f"Samples missing a mate file: {', '.join(unpaired)}. "
            "mothur's make.contigs requires both reads of every pair."
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        for sample in sorted(pairs):
            mates = pairs[sample]
            handle.write(f"{sample}\t{mates['1']}\t{mates['2']}\n")
    return output


def run_mothur_sop(
    *,
    reads_dir: Path,
    output_dir: Path,
    reference_alignment: Path,
    taxonomy_reference: Path,
    taxonomy_map: Path,
    identity: float = 0.97,
    force: bool = False,
    timeout: float | None = None,
) -> None:
    """Run make.contigs, then the mothur cluster and taxonomy backends."""
    if not reads_dir.is_dir():
        raise MicrobiomeSuiteError(f"Reads directory does not exist: {reads_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    run_dir = output_dir / "logs"

    work_dir = output_dir / "contigs"
    stability = write_stability_file(reads_dir, work_dir / "stability.files")

    outputs = run_mothur(
        "make.contigs",
        {"file": str(stability)},
        work_dir=work_dir,
        run_dir=run_dir / "make_contigs",
        timeout=timeout,
    )
    # make.contigs emits both .trim.contigs.fasta and .scrap.contigs.fasta;
    # the scrap file holds the reads that failed assembly.
    contigs = select_output(outputs, ".fasta", step="make.contigs", exclude=("scrap",))

    output_table = output_dir / "table.tsv"
    output_rep_seqs = output_dir / "rep-seqs.fasta"
    cluster(
        backend="mothur",
        rep_seqs=contigs,
        output_table=output_table,
        output_rep_seqs=output_rep_seqs,
        reference_alignment=reference_alignment,
        identity=identity,
        force=force,
        run_dir=run_dir / "cluster",
        timeout=timeout,
    )

    tax_classify(
        backend="mothur",
        rep_seqs=output_rep_seqs,
        output=output_dir / "taxonomy.tsv",
        taxonomy_reference=taxonomy_reference,
        taxonomy_map=taxonomy_map,
        force=force,
        run_dir=run_dir / "taxonomy",
        timeout=timeout,
    )
```

- [ ] **Step 4: Add the catalog entry**

Append to `WORKFLOWS` in `src/microsuite/workflows/catalog.py`:

```python
    WorkflowSpec(
        name="mothur",
        summary="Run the mothur MiSeq SOP from paired FASTQ files to an OTU table.",
        inputs=(
            "Paired FASTQ directory, an aligned reference FASTA, and a "
            "mothur trainset (reference FASTA + taxonomy file)."
        ),
        outputs="table.tsv, table.shared, rep-seqs.fasta, taxonomy.tsv, per-step logs.",
        status="ready",
    ),
```

- [ ] **Step 5: Add the CLI command**

In `src/microsuite/cli/workflow_cmd.py`, add the import and command:

```python
from microsuite.workflows.mothur_sop import run_mothur_sop
```

```python
@app.command("mothur")
def mothur_workflow(
    output: Annotated[Path, typer.Option("--out", "-o", help="Output run directory.")],
    reads_dir: Annotated[Path, typer.Option("--reads-dir", help="Paired FASTQ directory.")],
    reference_alignment: Annotated[
        Path, typer.Option("--reference-alignment", help="Aligned reference FASTA.")
    ],
    taxonomy_reference: Annotated[
        Path, typer.Option("--taxonomy-reference", help="Trainset reference FASTA.")
    ],
    taxonomy_map: Annotated[
        Path, typer.Option("--taxonomy-map", help="Trainset taxonomy file.")
    ],
    identity: Annotated[
        float, typer.Option("--identity", min=0.0, max=1.0, help="Clustering identity.")
    ] = 0.97,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
    timeout: Annotated[
        float | None, typer.Option("--timeout", help="Per-step timeout in seconds.")
    ] = None,
) -> None:
    run_mothur_sop(
        reads_dir=reads_dir,
        output_dir=output,
        reference_alignment=reference_alignment,
        taxonomy_reference=taxonomy_reference,
        taxonomy_map=taxonomy_map,
        identity=identity,
        force=force,
        timeout=timeout,
    )
```

- [ ] **Step 6: Run to verify it passes**

Run: `uv run pytest tests/test_mothur_workflow.py -v && uv run microsuite workflow list`
Expected: 4 passed; `mothur` appears in the workflow listing.

- [ ] **Step 7: Commit**

```bash
uv run ruff check . && uv run ruff format --check . && uv run ty check && uv run pytest -q
git add src/microsuite/workflows/mothur_sop.py src/microsuite/workflows/catalog.py src/microsuite/cli/workflow_cmd.py tests/test_mothur_workflow.py
git commit -m "feat(workflow): add end-to-end mothur MiSeq SOP workflow"
```

---

### Task 7: Documentation

The container and fixtures landed in Task 0. This task is documentation only.

**Files:**
- Modify: `docs/methods.md`
- Create: `docs/mothur.md`

- [ ] **Step 1: Update the method reference**

In `docs/methods.md`:

Add to the *Feature Table Generation* table (after line 77):

```markdown
| Generate OTU-style table | `mothur` | FASTA sequences with sample IDs in labels | OTU count table | table TSV, representative FASTA, sidecar `.shared` | mothur MiSeq SOP clustering; requires an aligned reference. |
```

Add to *Denoising And Clustering Backends* (after line 89):

```markdown
| `mothur` | mothur 1.48.2 | ready | `microsuite cluster --backend mothur --reference-alignment silva.v4.fasta` | `cluster(backend="mothur", rep_seqs=..., reference_alignment=...)` | [mothur](../containers/mothur/Dockerfile) or external `mothur` | Alignment-based OTU clustering with OptiClust; needs a user-supplied aligned reference and is slower than VSEARCH. | mothur MiSeq SOP OTU clustering. |
```

Add to *Taxonomy And Phylogeny* (after line 113):

```markdown
| `mothur` | mothur 1.48.2 | ready | `microsuite tax_classify --backend mothur --taxonomy-reference trainset.fasta --taxonomy-map trainset.tax` | `tax_classify(backend="mothur", rep_seqs=..., taxonomy_reference=..., taxonomy_map=...)` | [mothur](../containers/mothur/Dockerfile) or external `mothur` | Naive Bayes classification with optional per-OTU consensus; takes a reference FASTA + taxonomy pair rather than `--classifier`. | mothur classify.seqs / classify.otu. |
```

Add to the *Backend Validation Status* table (after line 32):

```markdown
| mothur | ready | Unit-tested wrapper + user environment | Command construction, stdout parsing, and option rejection are covered by Python tests; mothur itself and its reference data are user supplied. |
```

- [ ] **Step 2: Write the reference-data guide**

`docs/mothur.md` must cover, with working commands:
- Where to obtain the SILVA reference alignment for `--reference-alignment` (mothur.org's SILVA reference files) and how to trim it to a region with `pcr.seqs`.
- Where to obtain a trainset for `--taxonomy-reference` / `--taxonomy-map` (mothur.org's RDP or SILVA trainsets).
- A worked `microsuite workflow mothur` invocation with all five required paths.
- A note that these are user-supplied by design, cross-referencing `docs/superpowers/specs/2026-07-25-mothur-workflow-design.md`.

- [ ] **Step 3: Final verification and commit**

```bash
uv run ruff check . && uv run ruff format --check . && uv run ty check && uv run pytest -q
git add docs/methods.md docs/mothur.md
git commit -m "docs(mothur): document backends, validation status, and reference data"
```

---

## Self-Review

**Spec coverage:**

| Spec deliverable | Task |
|---|---|
| `methods/mothur.py` primitive + parser | 1, 2 |
| `cluster --backend mothur` | 4 |
| `tax_classify --backend mothur` | 5 |
| `microsuite workflow mothur` | 6 |
| `containers/mothur/Dockerfile` | 0 |
| Error contract: exit-0 `[ERROR]` scan | 1 (`check_mothur_errors`), 2 (wired into `run_mothur`) |
| Error contract: empty/absent output block | 1 (`select_output` raises "produced no output files") |
| Error contract: missing extension | 1 (`select_output` lists what was produced) |
| Error contract: reference validated before step 1 | 4 (`ensure_input` before the first `run_mothur`) |
| Error contract: step removes every sequence | 2 (`ensure_non_empty_fasta`), applied in 4 |
| `--identity` → cutoff conversion | 4 |
| `.shared` sidecar at `<output-table stem>.shared` | 4 |
| Option rejection, both directions | 5 |
| Testing: parser, command construction, rejection, transpose | 1, 2, 3, 4, 5 |
| Docs: methods.md, api-cli.md, mothur.md | 7 |
| CI smoke deferred | Not implemented — correct, per the spec's Deferred section |

`docs/api-cli.md` is listed in the spec's Documentation section but has no explicit step. It is regenerated content; Task 7 Step 1 covers `docs/methods.md`, and the api-cli entry follows the same edit. **Add it to Task 7 Step 1 if `docs/api-cli.md` turns out to be hand-maintained** — check whether a generator exists before editing by hand.

**Type consistency:** `run_mothur` returns `list[Path]` in Tasks 2, 4, 5, 6. `select_output(outputs, suffix, *, step, exclude)` keeps the same signature throughout. `write_otu_table_from_shared(shared, output)` matches between Tasks 3 and 4. `SUPPORTED_BACKENDS` (cluster) and `SUPPORTED_METHODS` (tax_classify) are the existing names, unchanged. `reject_options` is the public name in `_dispatch.py`; `trim.py` aliases it to its existing private name at import.

**Known risk, stated plainly:** the exact mothur parameter names in Task 4's pipeline table are written against mothur 1.48.2's documented interface but have **not** been executed. Task 0 verifies only `unique.seqs`, `align.seqs`, and `make.contigs`. The remaining eight SOP commands are first exercised for real by a user, not by this plan — the CI smoke test that would have caught parameter drift is deferred per the spec. Expect field reports.
