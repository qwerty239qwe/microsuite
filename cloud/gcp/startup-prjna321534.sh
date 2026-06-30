#!/usr/bin/env bash
set -Eeuo pipefail

metadata() {
  curl -fsH "Metadata-Flavor: Google" \
    "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1"
}

IMAGE="$(metadata image)"
BUCKET="$(metadata bucket)"
BIOPROJECT="$(metadata bioproject || true)"
BIOPROJECT="${BIOPROJECT:-PRJNA321534}"
MAX_RUNS="$(metadata max-runs || true)"
MAX_RUNS="${MAX_RUNS:-0}"
THREADS="$(metadata threads || true)"
THREADS="${THREADS:-$(nproc)}"
DELETE_ON_FINISH="$(metadata delete-on-finish || true)"
DELETE_ON_FINISH="${DELETE_ON_FINISH:-false}"
# raw = raw esearch/vsearch/breakaway plumbing (run_prjna321534);
# microsuite = drive the microsuite CLI (run_prjna321534_microsuite).
DRIVER="$(metadata driver || true)"
DRIVER="${DRIVER:-raw}"
if [[ "${DRIVER}" == "microsuite" ]]; then
  ENTRY="run_prjna321534_microsuite"
else
  ENTRY="run_prjna321534"
fi
WORK_ROOT="/mnt/microsuite-prjna321534"
RESULTS_GCS="gs://${BUCKET}/results/${BIOPROJECT}"

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a /var/log/prjna321534-startup.log
}

log "Installing Docker if needed"
if ! command -v docker >/dev/null 2>&1; then
  apt-get update
  apt-get install -y --no-install-recommends ca-certificates curl docker.io
fi
systemctl enable --now docker

mkdir -p "${WORK_ROOT}"

log "Pulling ${IMAGE}"
docker pull "${IMAGE}"

log "Running ${BIOPROJECT}; MAX_RUNS=${MAX_RUNS}; THREADS=${THREADS}; DRIVER=${DRIVER}"
set +e
docker run --rm \
  --name prjna321534-alpha \
  --entrypoint "${ENTRY}" \
  -e BIOPROJECT="${BIOPROJECT}" \
  -e MAX_RUNS="${MAX_RUNS}" \
  -e THREADS="${THREADS}" \
  -e WORKDIR=/work \
  -v "${WORK_ROOT}:/work" \
  "${IMAGE}" 2>&1 | tee -a /var/log/prjna321534-container.log
exit_code="${PIPESTATUS[0]}"
set -e

log "Syncing results to ${RESULTS_GCS}"
gcloud storage rsync -r "${WORK_ROOT}/results" "${RESULTS_GCS}" || true
gcloud storage cp /var/log/prjna321534-startup.log "${RESULTS_GCS}/logs/startup.log" || true
gcloud storage cp /var/log/prjna321534-container.log "${RESULTS_GCS}/logs/container.log" || true

if [[ "${exit_code}" -ne 0 ]]; then
  log "Container failed with exit code ${exit_code}"
else
  log "Completed successfully"
fi

# Stop billing immediately once results and logs are in GCS. Runs on both success
# and failure so a silently failing run does not idle until MAX_RUN_DURATION.
if [[ "${DELETE_ON_FINISH}" == "true" ]]; then
  name="$(curl -fsH 'Metadata-Flavor: Google' http://metadata.google.internal/computeMetadata/v1/instance/name)"
  zone="$(curl -fsH 'Metadata-Flavor: Google' http://metadata.google.internal/computeMetadata/v1/instance/zone)"
  zone="${zone##*/}"
  log "delete-on-finish set; deleting ${name} in ${zone}"
  gcloud compute instances delete "${name}" --zone "${zone}" --quiet || true
fi

exit "${exit_code}"
