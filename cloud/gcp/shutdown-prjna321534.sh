#!/usr/bin/env bash
set -Eeuo pipefail

metadata() {
  curl -fsH "Metadata-Flavor: Google" \
    "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1"
}

BUCKET="$(metadata bucket)"
BIOPROJECT="$(metadata bioproject || true)"
BIOPROJECT="${BIOPROJECT:-PRJNA321534}"
WORK_ROOT="/mnt/microsuite-prjna321534"
RESULTS_GCS="gs://${BUCKET}/results/${BIOPROJECT}"

if [[ -d "${WORK_ROOT}/results" ]]; then
  gcloud storage rsync -r "${WORK_ROOT}/results" "${RESULTS_GCS}" || true
fi
if [[ -f /var/log/prjna321534-startup.log ]]; then
  gcloud storage cp /var/log/prjna321534-startup.log "${RESULTS_GCS}/logs/startup.shutdown-copy.log" || true
fi
if [[ -f /var/log/prjna321534-container.log ]]; then
  gcloud storage cp /var/log/prjna321534-container.log "${RESULTS_GCS}/logs/container.shutdown-copy.log" || true
fi
