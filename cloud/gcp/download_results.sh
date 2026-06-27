#!/usr/bin/env bash
set -Eeuo pipefail

: "${BUCKET:?Set BUCKET without gs://}"

BIOPROJECT="${BIOPROJECT:-PRJNA321534}"
LOCAL_DIR="${LOCAL_DIR:-results/${BIOPROJECT}}"

mkdir -p "${LOCAL_DIR}"
gcloud storage rsync -r "gs://${BUCKET}/results/${BIOPROJECT}" "${LOCAL_DIR}"

echo "Downloaded results to ${LOCAL_DIR}"
