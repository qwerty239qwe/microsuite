#!/usr/bin/env bash
set -Eeuo pipefail

: "${PROJECT_ID:?Set PROJECT_ID}"

ZONE="${ZONE:-us-central1-a}"
VM_NAME="${VM_NAME:-microsuite-prjna321534-spot}"

gcloud config set project "${PROJECT_ID}" >/dev/null
gcloud compute instances delete "${VM_NAME}" --zone="${ZONE}"
