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
      "${force_args[@]+"${force_args[@]}"}" \
      >"${LOG_DIR}/${sample}.fastp.stdout.log" \
      2>"${LOG_DIR}/${sample}.fastp.stderr.log"
  else
    microsuite trim --backend fastp \
      --read1 "${read1}" \
      --output1 "${TRIMMED_DIR}/${sample}.fastq.gz" \
      --html "${FASTP_REPORT_DIR}/${sample}.fastp.html" \
      --json-report "${FASTP_REPORT_DIR}/${sample}.fastp.json" \
      --threads "${THREADS}" --run-dir "${RUN_DIR}/fastp/${sample}" \
      "${force_args[@]+"${force_args[@]}"}" \
      >"${LOG_DIR}/${sample}.fastp.stdout.log" \
      2>"${LOG_DIR}/${sample}.fastp.stderr.log"
  fi
}
export -f run_one_sample_by_id
export MANIFEST TRIMMED_DIR FASTP_REPORT_DIR RUN_DIR LOG_DIR THREADS FORCE

echo "Trimming samples (jobs=${JOBS}, threads=${THREADS})"
tail -n +2 "${MANIFEST}" | cut -f1 \
  | xargs -P "${JOBS}" -n 1 bash -c 'set -euo pipefail; run_one_sample_by_id "$1"' _

force_args=()
if [[ "${FORCE}" == "1" ]]; then force_args=(--force); fi

echo "Summarizing with MultiQC"
microsuite qc --backend multiqc \
  --input-dir "${FASTP_REPORT_DIR}" \
  --output-dir "${MULTIQC_DIR}" \
  --run-dir "${RUN_DIR}/multiqc" \
  "${force_args[@]+"${force_args[@]}"}" \
  >"${LOG_DIR}/multiqc.stdout.log" \
  2>"${LOG_DIR}/multiqc.stderr.log"

echo "Done"
echo "Accession:   ${ACCESSION}"
echo "Input:       ${INPUT_DIR}"
echo "Manifest:    ${MANIFEST}"
echo "Output:      ${OUT_DIR}"
echo "MultiQC:     ${MULTIQC_DIR}/multiqc_report.html"
