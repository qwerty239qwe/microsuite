# Portable `run_fastp_multiqc.sh` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a portable, parallel `scripts/run_fastp_multiqc.sh` that trims every FASTQ sample in one accession directory with fastp (via the microsuite CLI) and summarizes with MultiQC, following the repo's PATH-tools and opt-in-`--force` conventions.

**Architecture:** One self-contained bash script: parse args → build a `sample_id/layout/read1/read2` manifest with an inline `python3` block → preflight-check `microsuite`/`fastp`/`multiqc` on PATH → run per-sample `microsuite trim --backend fastp` concurrently via `xargs -P` (one sample_id per invocation, worker looks the row up in the manifest) → run `microsuite qc --backend multiqc` once. Tests target the script's observable contracts (manifest correctness, missing-tool guard, syntax) without running fastp/multiqc.

**Tech Stack:** Bash (portable, no GNU-only flags), inline `python3` for manifest detection, `xargs -P` for cross-sample parallelism, pytest + `subprocess` for tests, optional `shellcheck`.

## Global Constraints

- Single new file `scripts/run_fastp_multiqc.sh`, `#!/usr/bin/env bash`, `set -euo pipefail`.
- No `apt`/`dpkg`/`/tmp`-uv-env/`/mnt/c`/`LD_LIBRARY_PATH`/`rsync --no-perms`/`cp --no-preserve`.
- Required tools on PATH: `microsuite`, `fastp`, `multiqc`; missing any ⇒ actionable non-zero error citing `docs/installation.md`.
- `--force` is opt-in (default off); absent it, do NOT pass `--force` to microsuite.
- `--jobs N` (default 1), `--threads T` (default 4); `--jobs 1` == prior one-at-a-time behavior. MultiQC runs once, after all trims.
- Portable on macOS (BSD userland) and Linux (GNU): use `xargs -P N -n 1`; no `wait -n`, no `xargs -d`, no BSD/GNU-divergent flags.
- Manifest columns exactly `sample_id\tlayout\tread1\tread2`; `layout` is `PE` or `SE`; SE rows have an empty read2. `sample_id` must contain no whitespace (the builder rejects it) so it is safe to pass through `xargs`.
- New python test file starts with `from __future__ import annotations`.

---

### Task 1: The script + its tests

**Files:**
- Create: `scripts/run_fastp_multiqc.sh`
- Test: `tests/test_run_fastp_multiqc_script.py`

**Interfaces:**
- Consumes: the installed `microsuite` CLI (`trim --backend fastp`, `qc --backend multiqc`) — not modified here.
- Produces: an executable script with the CLI surface in Global Constraints; a manifest at `${OUTPUT_ROOT}/${ACCESSION}/fastp_outputs/fastq_manifest.tsv`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_run_fastp_multiqc_script.py
from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_fastp_multiqc.sh"


def _tiny_gz(path: Path) -> None:
    path.write_bytes(gzip.compress(b"@r\nACGT\n+\nIIII\n"))


def test_script_exists_and_syntax_ok() -> None:
    assert SCRIPT.exists(), SCRIPT
    result = subprocess.run(["bash", "-n", str(SCRIPT)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_shellcheck_clean() -> None:
    if shutil.which("shellcheck") is None:
        pytest.skip("shellcheck not installed")
    result = subprocess.run(["shellcheck", str(SCRIPT)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


def test_manifest_only_groups_pe_and_se(tmp_path: Path) -> None:
    fq = tmp_path / "acc"
    fq.mkdir()
    for name in ("x_R1.fastq.gz", "x_R2.fastq.gz", "y.fastq.gz"):
        _tiny_gz(fq / name)
    out = tmp_path / "results"
    result = subprocess.run(
        ["bash", str(SCRIPT), "acc", "--input-dir", str(fq),
         "--output-root", str(out), "--manifest-only"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    manifest = out / "acc" / "fastp_outputs" / "fastq_manifest.tsv"
    rows = [line.split("\t") for line in manifest.read_text().splitlines() if line]
    assert rows[0] == ["sample_id", "layout", "read1", "read2"]
    data = {row[0]: row for row in rows[1:]}
    assert data["x"][1] == "PE"
    assert data["x"][2].endswith("x_R1.fastq.gz")
    assert data["x"][3].endswith("x_R2.fastq.gz")
    assert data["y"][1] == "SE"
    assert data["y"][3] == ""


def test_missing_fastp_tool_guard(tmp_path: Path) -> None:
    fq = tmp_path / "acc"
    fq.mkdir()
    _tiny_gz(fq / "y.fastq.gz")
    # Build a PATH that has python3 + bash + stub microsuite/multiqc but NO fastp.
    shim = tmp_path / "bin"
    shim.mkdir()
    for tool in ("microsuite", "multiqc"):
        p = shim / tool
        p.write_text("#!/usr/bin/env bash\nexit 0\n")
        p.chmod(0o755)
    real_dirs = {str(Path(shutil.which("bash")).parent), str(Path(shutil.which("python3")).parent)}
    path = os.pathsep.join([str(shim), *real_dirs])
    env = dict(os.environ, PATH=path)
    # Guard the guard: ensure fastp really is unreachable on this constructed PATH.
    assert shutil.which("fastp", path=path) is None
    out = tmp_path / "results"
    result = subprocess.run(
        ["bash", str(SCRIPT), "acc", "--input-dir", str(fq), "--output-root", str(out)],
        capture_output=True, text=True, env=env,
    )
    assert result.returncode != 0
    assert "fastp" in result.stderr
```

- [ ] **Step 2: Run to verify they fail**

Run: `uv run pytest tests/test_run_fastp_multiqc_script.py -v`
Expected: FAIL (script does not exist yet).

- [ ] **Step 3: Create `scripts/run_fastp_multiqc.sh`**

```bash
#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_fastp_multiqc.sh ACCESSION [options]

Trim every FASTQ sample in one accession directory with fastp (via microsuite),
then summarize the fastp reports with MultiQC.

Options:
  --input-root DIR     Root containing accession directories (default: ./data)
  --input-dir DIR      FASTQ directory to use directly (overrides --input-root)
  --output-root DIR    Root for outputs (default: ./results)
  --jobs N             Samples to trim concurrently (default: 1)
  --threads T          fastp threads per sample (default: 4)
  --force              Overwrite existing outputs (default: off)
  --manifest-only      Detect FASTQ layout, write the manifest, then exit
  --help               Show this help

Total cores used is approximately jobs x threads. --jobs 1 (default) trims one
sample at a time. Requires microsuite, fastp, and multiqc on PATH
(see docs/installation.md).
EOF
}

INPUT_ROOT="./data"
INPUT_DIR=""
OUTPUT_ROOT="./results"
JOBS=1
THREADS=4
FORCE=0
MANIFEST_ONLY=0
ACCESSION=""

while (($#)); do
  case "$1" in
    --input-root) INPUT_ROOT="$2"; shift 2 ;;
    --input-dir) INPUT_DIR="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --manifest-only) MANIFEST_ONLY=1; shift ;;
    --help | -h) usage; exit 0 ;;
    --*) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    *)
      if [[ -n "${ACCESSION}" ]]; then
        echo "Unexpected argument: $1" >&2; usage >&2; exit 2
      fi
      ACCESSION="$1"; shift ;;
  esac
done

if [[ -z "${ACCESSION}" ]]; then
  echo "ACCESSION is required." >&2; usage >&2; exit 2
fi

if [[ -z "${INPUT_DIR}" ]]; then
  INPUT_DIR="${INPUT_ROOT}/${ACCESSION}"
fi
if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "Input directory does not exist: ${INPUT_DIR}" >&2; exit 1
fi

OUT_DIR="${OUTPUT_ROOT}/${ACCESSION}/fastp_outputs"
TRIMMED_DIR="${OUT_DIR}/trimmed_fastq"
FASTP_REPORT_DIR="${OUT_DIR}/fastp_reports"
MULTIQC_DIR="${OUT_DIR}/multiqc"
RUN_DIR="${OUT_DIR}/run_logs"
LOG_DIR="${OUT_DIR}/logs"
MANIFEST="${OUT_DIR}/fastq_manifest.tsv"

mkdir -p "${TRIMMED_DIR}" "${FASTP_REPORT_DIR}" "${MULTIQC_DIR}" "${RUN_DIR}" "${LOG_DIR}"

python3 - "${INPUT_DIR}" "${MANIFEST}" <<'PY'
from __future__ import annotations

import re
import sys
from pathlib import Path

input_dir = Path(sys.argv[1])
manifest = Path(sys.argv[2])

extensions = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
fastqs = sorted(
    p for p in input_dir.iterdir() if p.is_file() and p.name.endswith(extensions)
)
if not fastqs:
    raise SystemExit(f"No FASTQ files found in {input_dir}")

patterns = [
    re.compile(r"^(?P<sample>.+?)(?P<sep>[._-])R(?P<read>[12])(?:[._-]001)?$"),
    re.compile(r"^(?P<sample>.+?)(?P<sep>[._-])read(?P<read>[12])(?:[._-]001)?$", re.IGNORECASE),
    re.compile(r"^(?P<sample>.+?)(?P<sep>[._-])(?P<read>[12])(?:[._-]001)?$"),
]


def stem(path: Path) -> str:
    name = path.name
    for ext in extensions:
        if name.endswith(ext):
            return name[: -len(ext)]
    return path.stem


groups: dict[str, dict[str, Path]] = {}
singletons: list[Path] = []
for path in fastqs:
    s = stem(path)
    match = next((m for m in (p.match(s) for p in patterns) if m), None)
    if match is None:
        singletons.append(path)
        continue
    groups.setdefault(match.group("sample"), {})[match.group("read")] = path

rows: list[tuple[str, str, Path, Path | None]] = []
for sample, reads in sorted(groups.items()):
    if "1" in reads and "2" in reads:
        rows.append((sample, "PE", reads["1"], reads["2"]))
    elif "1" in reads:
        rows.append((sample, "SE", reads["1"], None))
    else:
        raise SystemExit(f"Found R2 without R1 for sample {sample}: {reads['2']}")
for path in singletons:
    rows.append((stem(path), "SE", path, None))

if not rows:
    raise SystemExit(f"No usable FASTQ samples found in {input_dir}")

for sample, _, _, _ in rows:
    if re.search(r"\s", sample):
        raise SystemExit(f"Sample id contains whitespace, unsupported: {sample!r}")

manifest.parent.mkdir(parents=True, exist_ok=True)
with manifest.open("w", encoding="utf-8") as handle:
    handle.write("sample_id\tlayout\tread1\tread2\n")
    for sample, layout, read1, read2 in sorted(rows):
        handle.write(f"{sample}\t{layout}\t{read1}\t{read2 or ''}\n")
print(f"Wrote manifest for {len(rows)} samples: {manifest}")
PY

if [[ "${MANIFEST_ONLY}" == "1" ]]; then
  echo "Manifest: ${MANIFEST}"
  exit 0
fi

for tool in microsuite fastp multiqc; do
  if ! command -v "${tool}" >/dev/null 2>&1; then
    echo "Required tool not found on PATH: ${tool}." >&2
    echo "Install it (see docs/installation.md; 'trim --backend fastp' needs fastp) or use the microsuite containers." >&2
    exit 1
  fi
done

run_one_sample_by_id() {
  local sample="$1"
  local line layout read1 read2
  line="$(awk -F '\t' -v s="${sample}" '$1 == s { print; exit }' "${MANIFEST}")"
  IFS=$'\t' read -r _ layout read1 read2 <<<"${line}"

  local force_args=()
  if [[ "${FORCE}" == "1" ]]; then force_args=(--force); fi

  if [[ "${layout}" == "PE" ]]; then
    microsuite trim --backend fastp \
      --read1 "${read1}" --read2 "${read2}" \
      --output1 "${TRIMMED_DIR}/${sample}_1.fastq.gz" \
      --output2 "${TRIMMED_DIR}/${sample}_2.fastq.gz" \
      --html "${FASTP_REPORT_DIR}/${sample}.fastp.html" \
      --json-report "${FASTP_REPORT_DIR}/${sample}.fastp.json" \
      --threads "${THREADS}" --run-dir "${RUN_DIR}/fastp/${sample}" \
      "${force_args[@]}" \
      >"${LOG_DIR}/${sample}.fastp.stdout.log" \
      2>"${LOG_DIR}/${sample}.fastp.stderr.log"
  else
    microsuite trim --backend fastp \
      --read1 "${read1}" \
      --output1 "${TRIMMED_DIR}/${sample}.fastq.gz" \
      --html "${FASTP_REPORT_DIR}/${sample}.fastp.html" \
      --json-report "${FASTP_REPORT_DIR}/${sample}.fastp.json" \
      --threads "${THREADS}" --run-dir "${RUN_DIR}/fastp/${sample}" \
      "${force_args[@]}" \
      >"${LOG_DIR}/${sample}.fastp.stdout.log" \
      2>"${LOG_DIR}/${sample}.fastp.stderr.log"
  fi
}
export -f run_one_sample_by_id
export MANIFEST TRIMMED_DIR FASTP_REPORT_DIR RUN_DIR LOG_DIR THREADS FORCE

echo "Trimming samples (jobs=${JOBS}, threads=${THREADS})"
tail -n +2 "${MANIFEST}" | cut -f1 \
  | xargs -P "${JOBS}" -n 1 bash -c 'run_one_sample_by_id "$1"' _

force_args=()
if [[ "${FORCE}" == "1" ]]; then force_args=(--force); fi

echo "Summarizing with MultiQC"
microsuite qc --backend multiqc \
  --input-dir "${FASTP_REPORT_DIR}" \
  --output-dir "${MULTIQC_DIR}" \
  --run-dir "${RUN_DIR}/multiqc" \
  "${force_args[@]}" \
  >"${LOG_DIR}/multiqc.stdout.log" \
  2>"${LOG_DIR}/multiqc.stderr.log"

echo "Done"
echo "Accession:   ${ACCESSION}"
echo "Input:       ${INPUT_DIR}"
echo "Manifest:    ${MANIFEST}"
echo "Output:      ${OUT_DIR}"
echo "MultiQC:     ${MULTIQC_DIR}/multiqc_report.html"
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_run_fastp_multiqc_script.py -v`
Expected: PASS (`test_shellcheck_clean` skips if `shellcheck` is absent). If `shellcheck` IS present and flags anything, fix the script until clean.

- [ ] **Step 5: Make the script executable + confirm syntax**

Run: `chmod +x scripts/run_fastp_multiqc.sh && bash -n scripts/run_fastp_multiqc.sh && echo OK`
Expected: `OK`

- [ ] **Step 6: Commit**

```bash
git add scripts/run_fastp_multiqc.sh tests/test_run_fastp_multiqc_script.py
git commit -m "feat(scripts): portable parallel run_fastp_multiqc.sh with PATH preflight"
```

---

## Self-Review

**Spec coverage:**
- New `scripts/run_fastp_multiqc.sh`, bash strict mode → Task 1 script. ✓
- No apt/dpkg/tmp-uv/mnt-c/LD_LIBRARY_PATH/rsync/cp-no-preserve → none present in the script. ✓
- PATH preflight for microsuite/fastp/multiqc with actionable error → the `for tool in ...` loop. ✓
- `--force` opt-in, not passed absent the flag → `FORCE=0` default, `force_args` only set when `FORCE==1`. ✓
- `--jobs`/`--threads`, MultiQC once after → `xargs -P "${JOBS}" -n 1`, multiqc after the pipeline. ✓
- Portable (no `wait -n`, `xargs -d`, GNU-only) → uses `xargs -P N -n 1`, awk, POSIX-ish bash. ✓
- Manifest columns + no-whitespace sample_id guard → python builder writes exact header and raises on whitespace. ✓
- Tests: manifest PE/SE, missing-tool guard, syntax + optional shellcheck → the four test functions. ✓

**Placeholder scan:** none — full script and full tests provided.

**Consistency:** manifest column order `sample_id/layout/read1/read2` is written by the python builder and read back by `run_one_sample_by_id` (`read -r _ layout read1 read2`) — aligned. The missing-tool test asserts `fastp` unreachable on the constructed PATH before running, so the guard is what triggers the failure.

**Note on the worker/xargs choice:** only `sample_id` (guaranteed whitespace-free) passes through `xargs`; the worker re-reads the full row from the manifest via `awk`, avoiding any field-splitting hazard on paths containing spaces. `xargs -P N -n 1` is the portable subset available on both BSD (macOS) and GNU.
