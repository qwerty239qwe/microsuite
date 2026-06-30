#!/usr/bin/env bash
set -Eeuo pipefail

# microsuite-driven variant of the PRJNA321534 alpha-diversity pipeline.
#
# The point of this pipeline is to exercise the microsuite CLI itself end to end
# (cluster -> import -> diversity -> workflow) on real SRA data, rather than the
# raw vsearch/breakaway plumbing in run_prjna321534.sh. SRA fetching stays raw
# because microsuite has no SRA fetcher.

BIOPROJECT="${BIOPROJECT:-PRJNA321534}"
WORKDIR="${WORKDIR:-/work}"
THREADS="${THREADS:-$(nproc)}"
MAX_RUNS="${MAX_RUNS:-0}"
OTU_ID="${OTU_ID:-0.97}"
SRA_MAX_SIZE="${SRA_MAX_SIZE:-200G}"

MANIFEST_DIR="${WORKDIR}/manifest"
FASTQ_DIR="${WORKDIR}/fastq"
FASTA_DIR="${WORKDIR}/fasta"
RESULTS_DIR="${WORKDIR}/results"
ALPHA_DIR="${RESULTS_DIR}/alpha"
LOG_DIR="${WORKDIR}/logs"
TMP_DIR="${WORKDIR}/tmp"

mkdir -p "${MANIFEST_DIR}" "${FASTQ_DIR}" "${FASTA_DIR}" "${RESULTS_DIR}" \
  "${ALPHA_DIR}" "${LOG_DIR}" "${TMP_DIR}"

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "${LOG_DIR}/run.log"
}

# --- 1. Select SRA runs (raw edirect; microsuite has no SRA fetcher) -----------
log "Fetching SRA RunInfo for ${BIOPROJECT}"
esearch -db sra -query "${BIOPROJECT}" \
  | efetch -format runinfo \
  > "${MANIFEST_DIR}/sra_runinfo.csv"

python3 - "${MANIFEST_DIR}/sra_runinfo.csv" "${MANIFEST_DIR}/runs.txt" "${MAX_RUNS}" <<'PY'
import csv
import sys

runinfo, output, max_runs = sys.argv[1], sys.argv[2], int(sys.argv[3])
with open(runinfo, newline="", encoding="utf-8") as handle:
    sample = handle.read(2048)
    handle.seek(0)
    if sample.lstrip().startswith("<"):
        raise SystemExit("SRA RunInfo fetch returned HTML, not CSV.")
    reader = csv.DictReader(handle)
    if "Run" not in (reader.fieldnames or []):
        raise SystemExit("SRA RunInfo CSV is missing the required Run column.")
    rows = [row for row in reader if row.get("Run")]
runs = [row["Run"] for row in rows]
if max_runs > 0:
    runs = runs[:max_runs]
if not runs:
    raise SystemExit("No SRA runs were found; refusing to launch an empty run.")
with open(output, "w", encoding="utf-8") as handle:
    handle.write("\n".join(runs) + "\n")
PY

RUN_COUNT="$(wc -l < "${MANIFEST_DIR}/runs.txt" | tr -d ' ')"
log "Selected ${RUN_COUNT} SRA runs"

# --- 2. Download each run and emit per-sample FASTA ----------------------------
fastq_to_fasta() {
  local sample="$1"
  local fastq="$2"
  local fasta="$3"
  python3 - "${sample}" "${fastq}" "${fasta}" <<'PY'
import gzip
import sys

sample, fastq, fasta = sys.argv[1:4]
open_in = gzip.open if fastq.endswith(".gz") else open
with open_in(fastq, "rt", encoding="utf-8", errors="replace") as inp, open(
    fasta, "w", encoding="utf-8"
) as out:
    index = 0
    while True:
        header = inp.readline()
        if not header:
            break
        seq = inp.readline().strip()
        inp.readline()
        inp.readline()
        index += 1
        out.write(f">{sample}_{index};sample={sample};\n{seq}\n")
PY
}

while read -r run; do
  [[ -z "${run}" ]] && continue
  sample_fasta="${FASTA_DIR}/${run}.fasta"
  if [[ -s "${sample_fasta}" ]]; then
    log "Skipping ${run}; FASTA already exists"
    continue
  fi

  log "Downloading ${run}"
  prefetch --max-size "${SRA_MAX_SIZE}" -O "${TMP_DIR}" "${run}" 2>&1 | tee -a "${LOG_DIR}/${run}.prefetch.log"
  fasterq-dump --threads "${THREADS}" --split-files -O "${FASTQ_DIR}" "${TMP_DIR}/${run}" \
    2>&1 | tee -a "${LOG_DIR}/${run}.fasterq.log"

  r1="${FASTQ_DIR}/${run}_1.fastq"
  r2="${FASTQ_DIR}/${run}_2.fastq"
  single="${FASTQ_DIR}/${run}.fastq"
  merged="${FASTQ_DIR}/${run}.merged.fastq"

  if [[ -s "${r1}" && -s "${r2}" ]]; then
    log "Merging paired reads for ${run}"
    vsearch --fastq_mergepairs "${r1}" \
      --reverse "${r2}" \
      --fastqout "${merged}" \
      --threads "${THREADS}" \
      2>&1 | tee -a "${LOG_DIR}/${run}.merge.log"
    fastq_to_fasta "${run}" "${merged}" "${sample_fasta}"
  elif [[ -s "${single}" ]]; then
    fastq_to_fasta "${run}" "${single}" "${sample_fasta}"
  else
    log "No FASTQ output found for ${run}; skipping"
  fi
done < "${MANIFEST_DIR}/runs.txt"

log "Building combined FASTA"
cat "${FASTA_DIR}"/*.fasta > "${RESULTS_DIR}/all_samples.fasta"

# --- 3. Sample metadata (one column so import has structure) -------------------
{
  printf 'sample\tbioproject\n'
  while read -r run; do
    [[ -z "${run}" ]] && continue
    printf '%s\t%s\n' "${run}" "${BIOPROJECT}"
  done < "${MANIFEST_DIR}/runs.txt"
} > "${RESULTS_DIR}/metadata.tsv"

# --- 4. Cluster into an OTU table via the microsuite CLI -----------------------
log "Clustering OTUs with: microsuite cluster --backend vsearch"
microsuite cluster \
  --backend vsearch \
  --rep-seqs "${RESULTS_DIR}/all_samples.fasta" \
  --output-table "${RESULTS_DIR}/otu_table.tsv" \
  --output-rep-seqs "${RESULTS_DIR}/otus.fasta" \
  --identity "${OTU_ID}" \
  --sample-delimiter '_' \
  --sample-field 0 \
  --run-dir "${LOG_DIR}/cluster" \
  2>&1 | tee -a "${LOG_DIR}/microsuite.cluster.log"

# --- 5. Import the table into the microsuite AnnData format --------------------
log "Importing OTU table with: microsuite import tsv"
microsuite import tsv "${RESULTS_DIR}/otu_table.tsv" \
  --metadata "${RESULTS_DIR}/metadata.tsv" \
  --output "${RESULTS_DIR}/table.h5ad" \
  2>&1 | tee -a "${LOG_DIR}/microsuite.import.log"

# --- 6. Alpha diversity via the microsuite CLI (R + native metrics) ------------
for metric in breakaway inext shannon observed_features chao1; do
  log "microsuite diversity alpha --metric ${metric}"
  if microsuite diversity alpha "${RESULTS_DIR}/table.h5ad" \
      --metric "${metric}" \
      --output "${ALPHA_DIR}/${metric}.tsv" \
      2>&1 | tee -a "${LOG_DIR}/microsuite.alpha.${metric}.log"; then
    log "alpha ${metric}: ok"
  else
    log "alpha ${metric}: FAILED (continuing)"
  fi
done

# --- 7. Whole-table summary workflow (alpha + beta + ordination + report) ------
log "microsuite workflow table-summary"
microsuite workflow table-summary \
  --out "${RESULTS_DIR}/summary" \
  --table "${RESULTS_DIR}/otu_table.tsv" \
  --metadata "${RESULTS_DIR}/metadata.tsv" \
  2>&1 | tee -a "${LOG_DIR}/microsuite.workflow.log" || log "table-summary: FAILED (continuing)"

cat > "${RESULTS_DIR}/run_summary.json" <<JSON
{
  "bioproject": "${BIOPROJECT}",
  "run_count": ${RUN_COUNT},
  "otu_id": "${OTU_ID}",
  "threads": ${THREADS},
  "driver": "microsuite-cli",
  "completed_at": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
JSON

log "Done. Results are in ${RESULTS_DIR}"
