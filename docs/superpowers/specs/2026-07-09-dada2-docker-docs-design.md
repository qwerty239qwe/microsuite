# DADA2 Docker documentation (P3) — Design

- **Date:** 2026-07-09
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** Codex's ERP120510 friction list, complaints #5 (Docker build path
  slow/network-fragile via BiocManager) and #6 (container boundaries confusing;
  "run microsuite with Docker" is ambiguous). **P3** of the four-part DADA2
  roadmap (see [[dada2-improvement-roadmap]]); P1 (docker runtime) and P2 (output
  validation) are merged.

## Scope

P3 is documentation only: add a "Running DADA2 in a container" section to
`docs/dada2.md` documenting the three container paths (the `--runtime docker`
option from P1, the whole-CLI `microsuite-dada2` image, and Nextflow), how to
pull the prebuilt GHCR image instead of building, and the container-boundary
model — plus two cross-links and a presence test.

### Out of scope for P3
- Any code, CLI, container, or workflow change.
- The naming-contract end-to-end test (#10) — **P4**.
- Rewriting `docs/containers.md` (it already covers build/registry details; P3
  cross-links to it).

## Verified context

- `docs/dada2.md` exists with sections: `# DADA2 denoising`, `## Backend
  choice`, `## Primers, ITS, and multi-run studies`, `## Diagnostics`. It has
  **no** Docker section.
- `docs/containers.md` already documents the container-boundary model (r-dada2 =
  R+dada2 runtime-only; microsuite-dada2 = CLI + packaged `dada2_denoise.R` in
  one image, "runs end-to-end in one container"), the pinned Bioconductor base,
  and local `docker build` commands — but nothing user-facing about *running*
  DADA2 via containers.
- P1 shipped `microsuite denoise --backend dada2-r --runtime local|docker` with
  `--image` and env `MICROSUITE_R_DADA2_IMAGE` (default
  `ghcr.io/qwerty239qwe/microsuite/r-dada2:latest`). Not yet in user docs.
- The DADA2 images publish to GHCR via the Docker CI workflow; `sha-<commit>`
  tags exist for reproducible pinning.

## Design

### New section in `docs/dada2.md`: "Running DADA2 in a container"

Placed after `## Backend choice` (before `## Diagnostics`), with these parts:

1. **No local R required.** One sentence: the `dada2-r` backend needs R +
   `dada2`; the container options below mean a pip-installed microsuite user
   never installs R locally.

2. **Three container paths:**
   - **`--runtime docker` (recommended).** microsuite runs locally and delegates
     the R step to the `r-dada2` container:
     ```bash
     microsuite denoise --backend dada2-r --runtime docker \
       --demux reads_dir --mode paired \
       --output-table asv_table.tsv --output-rep-seqs rep_seqs.fasta \
       --output-stats stats.tsv --output-plot-dir plots
     ```
     Override the image with `--image` or the `MICROSUITE_R_DADA2_IMAGE` env var.
   - **Whole CLI in one container.** `docker run` the `microsuite-dada2` image
     (bundles the CLI + packaged R script), mounting the data dir:
     ```bash
     docker run --rm -v "$PWD:/work" -w /work \
       ghcr.io/qwerty239qwe/microsuite/microsuite-dada2:latest \
       denoise --backend dada2-r --demux reads_dir --mode paired \
       --output-table asv_table.tsv --output-rep-seqs rep_seqs.fasta \
       --output-stats stats.tsv
     ```
     This is what "run microsuite with Docker" means for the fully-in-container path.
   - **Nextflow.** The `amplicon_qiime2` workflow already runs its DADA2 in the
     qiime2 container; see [multisample runs and concurrency](multisample.md).

3. **Pull the prebuilt image; don't build it (fixes #5).**
   ```bash
   docker pull ghcr.io/qwerty239qwe/microsuite/r-dada2:latest
   ```
   Building the DADA2 image locally installs `dada2` through BiocManager, which
   is slow and network-fragile (low CPU/memory can look like a frozen build).
   Pull the published GHCR image instead; pin a `sha-<commit>` tag for a
   reproducible run.

4. **Container boundaries (fixes #6).** Short paragraph: `r-dada2` = R + `dada2`
   only (what `--runtime docker` invokes); `microsuite-dada2` = the CLI + R
   backend in one image (the whole-CLI path). Cross-link
   [containers.md](containers.md) for build/registry/pinning details.

### Cross-links (one line each)

- `docs/installation.md` — near the Docker users row: a pointer to the new
  DADA2 container section.
- `docs/containers.md` — near the r-dada2 / microsuite-dada2 image entries: a
  pointer to `dada2.md`'s "Running DADA2 in a container" for the how-to.

## Testing

New `tests/test_dada2_docs.py` (string/file-presence, matching the repo's
doc-test convention):
- `docs/dada2.md` contains the Docker section (assert distinctive substrings:
  `--runtime docker`, `docker pull`, `r-dada2`, `MICROSUITE_R_DADA2_IMAGE`).
- `docs/installation.md` links to `dada2.md` (the relative path appears).

## Success criteria

1. `docs/dada2.md` has a "Running DADA2 in a container" section covering the
   three container paths, the prebuilt-image pull (with the BiocManager caveat),
   and the container-boundary model, matching the as-built P1 CLI
   (`--runtime docker`, `--image`, `MICROSUITE_R_DADA2_IMAGE`).
2. `docs/installation.md` and `docs/containers.md` link to it.
3. `tests/test_dada2_docs.py` asserts the section exists and the installation
   cross-link is present; the full offline suite stays green.

## Open questions / follow-ups (not blocking P3)

- If a `docker pull`-and-verify smoke check is ever wanted, it would be opt-in
  (network); P3 keeps to presence tests.
