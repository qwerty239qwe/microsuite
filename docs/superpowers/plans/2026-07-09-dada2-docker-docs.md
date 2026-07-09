# DADA2 Docker Documentation (P3) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a "Running DADA2 in a container" section to `docs/dada2.md` documenting the three container paths (P1's `--runtime docker`, the whole-CLI `microsuite-dada2` image, Nextflow), pulling the prebuilt GHCR image instead of building, and the container-boundary model — with two cross-links and a presence test.

**Architecture:** One new prose section in an existing doc, two additive one-line cross-links, and a small string/file-presence pytest. Documentation only.

**Tech Stack:** Markdown; pytest + `pathlib` presence test.

## Global Constraints

- Documentation only: no changes to `src/`, `scripts/`, `workflows/`, or `containers/`.
- Content must match the as-built P1 CLI: `microsuite denoise --backend dada2-r --runtime local|docker`, `--image`, env `MICROSUITE_R_DADA2_IMAGE`, default image `ghcr.io/qwerty239qwe/microsuite/r-dada2:latest`.
- Cross-links are additive one-liners; do not remove or restructure existing content.
- New python test file starts with `from __future__ import annotations`.

---

### Task 1: DADA2 container section + cross-links + presence test

**Files:**
- Modify: `docs/dada2.md` (add the section after `## Backend choice`)
- Modify: `docs/installation.md` (add one link line near the Docker users row)
- Modify: `docs/containers.md` (add one link line near the r-dada2 / microsuite-dada2 entries)
- Test: `tests/test_dada2_docs.py`

**Interfaces:**
- Consumes: nothing (documentation).
- Produces: a "Running DADA2 in a container" section in `docs/dada2.md`; an `installation.md` link containing `dada2.md`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_dada2_docs.py
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_dada2_doc_has_container_section() -> None:
    text = (ROOT / "docs" / "dada2.md").read_text(encoding="utf-8")
    assert "Running DADA2 in a container" in text
    assert "--runtime docker" in text
    assert "docker pull" in text
    assert "r-dada2" in text
    assert "MICROSUITE_R_DADA2_IMAGE" in text


def test_installation_links_to_dada2_doc() -> None:
    text = (ROOT / "docs" / "installation.md").read_text(encoding="utf-8")
    assert "dada2.md" in text
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_dada2_docs.py -v`
Expected: FAIL (no container section; no installation link).

- [ ] **Step 3: Add the section to `docs/dada2.md`**

Insert this section immediately after the `## Backend choice` section (before the next `##` heading):

````markdown
## Running DADA2 in a container

The `dada2-r` backend needs R and the `dada2` package. You do **not** have to
install R locally — run the R step in a container instead. There are three ways:

**1. `--runtime docker` (recommended).** microsuite runs locally (pip-installed)
and delegates only the R step to the `r-dada2` container:

```bash
microsuite denoise --backend dada2-r --runtime docker \
  --demux reads_dir --mode paired \
  --output-table asv_table.tsv --output-rep-seqs rep_seqs.fasta \
  --output-stats stats.tsv --output-plot-dir plots
```

Override the image with `--image IMAGE` or the `MICROSUITE_R_DADA2_IMAGE`
environment variable (default: `ghcr.io/qwerty239qwe/microsuite/r-dada2:latest`).

**2. The whole CLI in one container.** The `microsuite-dada2` image bundles the
CLI and the packaged `dada2_denoise.R`, so the entire command runs inside one
container — mount your working directory:

```bash
docker run --rm -v "$PWD:/work" -w /work \
  ghcr.io/qwerty239qwe/microsuite/microsuite-dada2:latest \
  denoise --backend dada2-r --demux reads_dir --mode paired \
  --output-table asv_table.tsv --output-rep-seqs rep_seqs.fasta \
  --output-stats stats.tsv
```

**3. Nextflow.** The `amplicon_qiime2` workflow runs its DADA2 in the QIIME 2
container automatically; see [multisample runs and concurrency](multisample.md).

### Pull the prebuilt image; don't build it

```bash
docker pull ghcr.io/qwerty239qwe/microsuite/r-dada2:latest
```

Building the DADA2 image locally installs `dada2` through BiocManager, which is
slow and network-fragile — low CPU/memory during the Bioconductor downloads can
look like a frozen build. Pull the published GHCR image instead. For a
reproducible run, pin a `sha-<commit>` tag rather than `latest`.

### Container boundaries

Two images back the DADA2 path, and they are not the same thing:

- **`r-dada2`** — R + `dada2` only (no microsuite). This is what
  `--runtime docker` invokes for the R step.
- **`microsuite-dada2`** — the microsuite CLI + R backend in one image. This is
  the "run the whole command in Docker" path (option 2 above).

See [containers.md](containers.md) for how these images are built, published to
GHCR, and pinned.
````

- [ ] **Step 4: Add the two cross-links**

- In `docs/installation.md`, near the "Docker users" row / Docker section, add an additive line:
  ```markdown
  Running DADA2 without installing R? See [Running DADA2 in a container](dada2.md#running-dada2-in-a-container).
  ```
- In `docs/containers.md`, near the `r-dada2` / `microsuite-dada2` image entries, add an additive line:
  ```markdown
  For how to run DADA2 with these images, see [Running DADA2 in a container](dada2.md#running-dada2-in-a-container).
  ```

- [ ] **Step 5: Run to verify pass**

Run: `uv run pytest tests/test_dada2_docs.py -v`
Expected: PASS (both tests).

- [ ] **Step 6: Commit**

```bash
git add docs/dada2.md docs/installation.md docs/containers.md tests/test_dada2_docs.py
git commit -m "docs(dada2): document running DADA2 in a container (--runtime docker, prebuilt image, boundaries)"
```

---

## Self-Review

**Spec coverage:**
- "Running DADA2 in a container" section with the three paths + prebuilt-image pull (BiocManager caveat) + container-boundary model → Task 1 Step 3. ✓
- Matches as-built P1 CLI (`--runtime docker`, `--image`, `MICROSUITE_R_DADA2_IMAGE`, GHCR default) → verbatim in the section. ✓
- Cross-links from installation.md + containers.md → Task 1 Step 4. ✓
- Presence test (section substrings + installation link) → Task 1 Steps 1/5. ✓

**Placeholder scan:** none — full section content, exact cross-link lines, and full test provided.

**Consistency:** the test's asserted substrings (`Running DADA2 in a container`, `--runtime docker`, `docker pull`, `r-dada2`, `MICROSUITE_R_DADA2_IMAGE`) all appear verbatim in the Step 3 section; the installation cross-link line contains `dada2.md`, satisfying the `dada2.md` substring assertion.
