#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Fall back to the settings recorded by the last create_spot_vm.sh run so this
# works in a fresh shell without re-exporting BUCKET/BIOPROJECT.
if [[ -z "${BUCKET:-}" && -f "${SCRIPT_DIR}/.last_run.env" ]]; then
  # shellcheck disable=SC1091
  source "${SCRIPT_DIR}/.last_run.env"
fi

: "${BUCKET:?Set BUCKET without gs://}"

BIOPROJECT="${BIOPROJECT:-PRJNA321534}"
LOCAL_DIR="${LOCAL_DIR:-results/${BIOPROJECT}}"

mkdir -p "${LOCAL_DIR}"
gcloud storage rsync -r "gs://${BUCKET}/results/${BIOPROJECT}" "${LOCAL_DIR}"

echo "Downloaded results to ${LOCAL_DIR}"
