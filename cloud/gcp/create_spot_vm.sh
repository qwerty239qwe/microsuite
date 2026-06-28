#!/usr/bin/env bash
set -Eeuo pipefail

: "${PROJECT_ID:?Set PROJECT_ID, e.g. my-gcp-project}"
: "${BUCKET:?Set BUCKET without gs://, e.g. microsuite-prjna321534-$USER}"
: "${IMAGE:?Set IMAGE, e.g. ghcr.io/qwerty239qwe/microsuite-prjna321534-alpha:latest}"

ZONE="${ZONE:-us-central1-a}"
REGION="${REGION:-us-central1}"
VM_NAME="${VM_NAME:-microsuite-prjna321534-spot}"
MACHINE_TYPE="${MACHINE_TYPE:-e2-standard-16}"
BOOT_DISK_GB="${BOOT_DISK_GB:-200}"
BIOPROJECT="${BIOPROJECT:-PRJNA321534}"
MAX_RUNS="${MAX_RUNS:-10}"
THREADS="${THREADS:-16}"
TERMINATION_ACTION="${TERMINATION_ACTION:-STOP}"
MAX_RUN_DURATION="${MAX_RUN_DURATION:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

gcloud config set project "${PROJECT_ID}" >/dev/null

if ! gcloud storage buckets describe "gs://${BUCKET}" >/dev/null 2>&1; then
  gcloud storage buckets create "gs://${BUCKET}" --location="${REGION}"
fi

CREATE_ARGS=(--instance-termination-action="${TERMINATION_ACTION}")
if [[ -n "${MAX_RUN_DURATION}" ]]; then
  CREATE_ARGS+=(--max-run-duration="${MAX_RUN_DURATION}")
fi

gcloud compute instances create "${VM_NAME}" \
  --project="${PROJECT_ID}" \
  --zone="${ZONE}" \
  --machine-type="${MACHINE_TYPE}" \
  --provisioning-model=SPOT \
  "${CREATE_ARGS[@]}" \
  --boot-disk-size="${BOOT_DISK_GB}GB" \
  --boot-disk-type=pd-balanced \
  --image-family=ubuntu-2404-lts-amd64 \
  --image-project=ubuntu-os-cloud \
  --scopes=https://www.googleapis.com/auth/cloud-platform \
  --metadata="image=${IMAGE},bucket=${BUCKET},bioproject=${BIOPROJECT},max-runs=${MAX_RUNS},threads=${THREADS}" \
  --metadata-from-file="startup-script=${SCRIPT_DIR}/startup-prjna321534.sh,shutdown-script=${SCRIPT_DIR}/shutdown-prjna321534.sh" \
  --tags=microsuite-prjna321534

cat <<EOF

Created Spot VM: ${VM_NAME}

Watch logs:
  gcloud compute ssh ${VM_NAME} --zone ${ZONE} --command 'sudo tail -f /var/log/prjna321534-startup.log'

Results will sync to:
  gs://${BUCKET}/results/${BIOPROJECT}

Download when complete:
  BUCKET=${BUCKET} BIOPROJECT=${BIOPROJECT} LOCAL_DIR=results/${BIOPROJECT} cloud/gcp/download_results.sh
EOF
