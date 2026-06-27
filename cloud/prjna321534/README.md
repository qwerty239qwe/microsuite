# PRJNA321534 Cloud Alpha-Diversity Run

This is an optional manual real-data validation kit. It can incur Google Cloud
costs and is not run by default CI.

This kit builds a Docker image in GitHub Actions, pulls it on a Google Cloud
Spot VM, runs a VSEARCH OTU-table workflow for BioProject `PRJNA321534`, then
computes Breakaway and iNEXT alpha-diversity estimates.

## 1. Build and Push the Image

Run the GitHub Actions workflow:

```text
PRJNA321534 Cloud Image
```

The workflow is `workflow_dispatch` only. It will not run automatically on
normal pushes.

It pushes:

```text
ghcr.io/<github-owner>/microsuite-prjna321534-alpha:latest
ghcr.io/<github-owner>/microsuite-prjna321534-alpha:<commit-sha>
```

Use the `<commit-sha>` tag for real runs so the pilot and full run use the
same image. Use `latest` only for quick experiments. If the repository or
package is private, make sure the VM can pull from GHCR. For private images,
run `docker login ghcr.io` on the VM or use a public image.

## 2. Pilot Run on a Spot VM

Run a 10-run pilot first:

```bash
export PROJECT_ID="my-gcp-project"
export BUCKET="microsuite-prjna321534-$USER"
export IMAGE="ghcr.io/qwerty239qwe/microsuite-prjna321534-alpha:<commit-sha>"
export ZONE="us-central1-a"
export REGION="us-central1"
export MACHINE_TYPE="e2-standard-16"
export BOOT_DISK_GB="200"
export MAX_RUNS="10"
export THREADS="16"

bash cloud/gcp/create_spot_vm.sh
```

This creates billable Google Cloud resources. The default pilot disk is 200 GB;
increase it only after the pilot tells you the real storage requirement.

Watch logs:

```bash
gcloud compute ssh microsuite-prjna321534-spot \
  --zone "$ZONE" \
  --command 'sudo tail -f /var/log/prjna321534-startup.log'
```

Download results:

```bash
export LOCAL_DIR="results/PRJNA321534-pilot"
bash cloud/gcp/download_results.sh
```

Delete the VM when done:

```bash
bash cloud/gcp/delete_spot_vm.sh
```

## 3. Full Run

After the pilot succeeds, run the full BioProject:

```bash
export MAX_RUNS="0"
export VM_NAME="microsuite-prjna321534-full"
export BOOT_DISK_GB="2000"
bash cloud/gcp/create_spot_vm.sh
```

`MAX_RUNS=0` means all runs found in the SRA RunInfo file.

Keep the pilot disk small. A stopped Spot VM still bills for its persistent
disk until you delete the VM.

## Outputs

The VM syncs results to:

```text
gs://<bucket>/results/PRJNA321534/
```

Important files:

```text
manifest/sra_runinfo.csv
manifest/runs.txt
results/all_samples.fasta
results/uniques.fasta
results/otus.fasta
results/otu_table.tsv
results/alpha/breakaway_alpha.tsv
results/alpha/inext_alpha.tsv
results/run_summary.json
logs/
```

## Notes

- Spot VMs can stop at any time. The startup and shutdown scripts try to sync
  completed results and logs back to Cloud Storage.
- Always delete stopped Spot VMs when done. The persistent disk continues to
  bill while the VM exists, even if the VM was stopped by preemption.
- The workflow is intentionally simple and reproducible: SRA Toolkit downloads
  reads, VSEARCH clusters OTUs at `OTU_ID`, and R computes Breakaway/iNEXT.
- For final publication-grade analysis, inspect the BioProject metadata and
  decide whether primer trimming, sample filtering, or a stricter ASV workflow
  is needed before relying on the OTU table.
