#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Fall back to the settings recorded by the last create_spot_vm.sh run so this
# works in a fresh shell without re-exporting PROJECT_ID/ZONE/VM_NAME.
if [[ -z "${PROJECT_ID:-}" && -f "${SCRIPT_DIR}/.last_run.env" ]]; then
  # shellcheck disable=SC1091
  source "${SCRIPT_DIR}/.last_run.env"
fi

: "${PROJECT_ID:?Set PROJECT_ID}"

ZONE="${ZONE:-us-central1-a}"
VM_NAME="${VM_NAME:-microsuite-prjna321534-spot}"

gcloud config set project "${PROJECT_ID}" >/dev/null
gcloud compute instances delete "${VM_NAME}" --zone="${ZONE}"
