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
# SPOT is cheapest; switch to STANDARD for long full runs that you do not want
# preempted (the pipeline resumes downloads across a SPOT STOP, but STANDARD
# avoids preemption entirely).
PROVISIONING_MODEL="${PROVISIONING_MODEL:-SPOT}"
# On a SPOT preemption: STOP keeps the boot disk so the run resumes; DELETE
# discards it. Use STOP for full runs so per-sample progress survives.
TERMINATION_ACTION="${TERMINATION_ACTION:-STOP}"
MAX_RUN_DURATION="${MAX_RUN_DURATION:-}"
# When true, the VM deletes itself after results sync so a finished or failed
# run stops billing immediately instead of idling until MAX_RUN_DURATION.
DELETE_ON_FINISH="${DELETE_ON_FINISH:-false}"
# raw = raw esearch/vsearch/breakaway plumbing; microsuite = drive the
# microsuite CLI (cluster/import/diversity/workflow) on the same data.
DRIVER="${DRIVER:-raw}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

gcloud config set project "${PROJECT_ID}" >/dev/null

if ! gcloud storage buckets describe "gs://${BUCKET}" >/dev/null 2>&1; then
  gcloud storage buckets create "gs://${BUCKET}" --location="${REGION}"
fi

# --instance-termination-action is only valid for SPOT, or for any VM that also
# sets --max-run-duration.
CREATE_ARGS=(--provisioning-model="${PROVISIONING_MODEL}")
if [[ -n "${MAX_RUN_DURATION}" ]]; then
  CREATE_ARGS+=(--max-run-duration="${MAX_RUN_DURATION}" --instance-termination-action="${TERMINATION_ACTION}")
elif [[ "${PROVISIONING_MODEL}" == "SPOT" ]]; then
  CREATE_ARGS+=(--instance-termination-action="${TERMINATION_ACTION}")
fi

gcloud compute instances create "${VM_NAME}" \
  --project="${PROJECT_ID}" \
  --zone="${ZONE}" \
  --machine-type="${MACHINE_TYPE}" \
  "${CREATE_ARGS[@]}" \
  --boot-disk-size="${BOOT_DISK_GB}GB" \
  --boot-disk-type=pd-balanced \
  --image-family=ubuntu-2404-lts-amd64 \
  --image-project=ubuntu-os-cloud \
  --scopes=https://www.googleapis.com/auth/cloud-platform \
  --metadata="image=${IMAGE},bucket=${BUCKET},bioproject=${BIOPROJECT},max-runs=${MAX_RUNS},threads=${THREADS},delete-on-finish=${DELETE_ON_FINISH},driver=${DRIVER}" \
  --metadata-from-file="startup-script=${SCRIPT_DIR}/startup-prjna321534.sh,shutdown-script=${SCRIPT_DIR}/shutdown-prjna321534.sh" \
  --tags=microsuite-prjna321534

# Record this run's settings so download_results.sh / delete_spot_vm.sh work in a
# fresh shell without re-exporting every variable.
cat > "${SCRIPT_DIR}/.last_run.env" <<ENV
PROJECT_ID=${PROJECT_ID}
ZONE=${ZONE}
REGION=${REGION}
VM_NAME=${VM_NAME}
BUCKET=${BUCKET}
BIOPROJECT=${BIOPROJECT}
ENV

cat <<EOF

Created Spot VM: ${VM_NAME}

Watch logs:
  gcloud compute ssh ${VM_NAME} --zone ${ZONE} --command 'sudo tail -f /var/log/prjna321534-startup.log'

Results will sync to:
  gs://${BUCKET}/results/${BIOPROJECT}

Download when complete:
  BUCKET=${BUCKET} BIOPROJECT=${BIOPROJECT} LOCAL_DIR=results/${BIOPROJECT} cloud/gcp/download_results.sh
EOF
