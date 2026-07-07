# Multisample & Concurrency Guide Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `docs/multisample.md` documenting the shared sample manifest, the within-sample-threads vs cross-sample-parallelism model, and a decision guide (CLI vs script vs Nextflow), with cross-links from the main docs and a presence test.

**Architecture:** One new prose guide plus four one-line cross-links into existing docs plus a small string/file-presence pytest. Documentation only — no source, CLI, script, or workflow behavior changes.

**Tech Stack:** Markdown; pytest + `pathlib` presence test (matches the repo's doc/skeleton test style).

## Global Constraints

- Documentation only: no changes to `src/`, `scripts/`, or `workflows/` behavior.
- Facts must match the as-built A (`scripts/run_fastp_multiqc.sh`) and B (Nextflow `amplicon_qiime2`): manifest columns `sample_id/layout/read1/read2` (SE rows have empty `read2`); the workflow reads `sample_id/read1/read2` and ignores the extra `layout` column; cross-sample parallelism exists only in the script's `--jobs` and Nextflow's scheduler; `--threads auto` = CPUs − 1; rule of thumb total cores ≈ `jobs × threads`.
- Cross-links are additive one-liners; do not remove or restructure existing doc content.
- New python test file starts with `from __future__ import annotations`.

---

### Task 1: `docs/multisample.md` + cross-links + presence test

**Files:**
- Create: `docs/multisample.md`
- Modify: `README.md` (add one link line near "Three ways to use it")
- Modify: `docs/methods.md` (add one link line near the trim/QC entries)
- Modify: `docs/installation.md` (add one link line near the Nextflow/Docker rows)
- Modify: `docs/api-nextflow.md` (add one link line where the manifest/`amplicon_qiime2` workflow is described)
- Test: `tests/test_multisample_docs.py`

**Interfaces:**
- Consumes: nothing (documentation).
- Produces: `docs/multisample.md`; a README link containing the relative path `multisample.md`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_multisample_docs.py
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_multisample_doc_exists_with_sections() -> None:
    doc = ROOT / "docs" / "multisample.md"
    assert doc.exists(), doc
    text = doc.read_text(encoding="utf-8")
    # manifest section
    assert "sample_id" in text
    # concurrency guidance
    assert ("jobs × threads" in text) or ("jobs x threads" in text)
    # decision guide references all three paths
    assert "run_fastp_multiqc.sh" in text
    assert "amplicon_qiime2" in text


def test_readme_links_to_multisample_doc() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "multisample.md" in readme
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_multisample_docs.py -v`
Expected: FAIL (`docs/multisample.md` missing; README has no link).

- [ ] **Step 3: Create `docs/multisample.md`**

````markdown
# Multisample runs and concurrency

microsuite has three ways to process many samples. They share one **sample
manifest** format and one **concurrency model**; this guide explains both and
helps you pick a path.

## The sample manifest

Batched runs are driven by a tab-separated manifest with these columns:

| Column | Meaning |
| --- | --- |
| `sample_id` | Unique sample name (no whitespace). Names the outputs. |
| `layout` | `PE` (paired-end) or `SE` (single-end). |
| `read1` | Forward (or single) FASTQ path. |
| `read2` | Reverse FASTQ path; empty for single-end. |

Example:

```tsv
sample_id	layout	read1	read2
sampleA	PE	/data/acc/sampleA_R1.fastq.gz	/data/acc/sampleA_R2.fastq.gz
sampleB	SE	/data/acc/sampleB.fastq.gz	
```

`scripts/run_fastp_multiqc.sh` writes this file — run it with `--manifest-only`
to generate the manifest and stop. The Nextflow `amplicon_qiime2` workflow reads
the same `sample_id`/`read1`/`read2` columns (it ignores the extra `layout`
column and treats an empty `read2` as single-end), so a manifest built for the
script also works with the workflow.

## Concurrency: two independent axes

**Within a sample — tool threads.** `microsuite trim`/`qc` pass `--threads` to
the underlying tool (`fastp --thread`, `cutadapt -j`, ...). `--threads auto`
uses the detected CPU count minus one reserved core. fastp scales sublinearly
past roughly 4–8 threads, so piling every core onto a single sample wastes most
of them.

**Across samples — parallel jobs.** A bare `microsuite trim`/`qc` call processes
**one** sample; the CLI itself has no cross-sample parallelism. Two paths add
it:

- `scripts/run_fastp_multiqc.sh --jobs N` runs up to `N` samples at once via
  `xargs -P`.
- The Nextflow workflow runs one process per sample and the executor schedules
  them concurrently (local = host CPUs; HPC/cloud = the cluster).

**Rule of thumb:** total cores used ≈ `jobs × threads`. For many small samples,
prefer more `jobs` with fewer `threads` each; for a few large samples, fewer
jobs with more threads.

## Which path to use

| Situation | Use | What it gives you |
| --- | --- | --- |
| One sample, or a few local files | `microsuite trim` / `qc` directly | Simplest; one task at a time |
| Batched trimming on a single machine | `scripts/run_fastp_multiqc.sh --jobs N` | Cross-sample parallelism on one box |
| Reproducible multisample pipeline (HPC/cloud) | Nextflow `amplicon_qiime2` (`--trim`) | Scheduler parallelism, pinned containers, provenance |

## Examples

Single sample, one CLI call:

```bash
microsuite trim --backend fastp \
  --read1 reads_R1.fastq.gz --read2 reads_R2.fastq.gz \
  --output1 trimmed_R1.fastq.gz --output2 trimmed_R2.fastq.gz \
  --html sample.fastp.html --json-report sample.fastp.json --threads auto
```

Batched trimming on one machine (4 samples at a time, 4 threads each):

```bash
bash scripts/run_fastp_multiqc.sh ACCESSION --input-dir /data/ACCESSION --jobs 4 --threads 4
```

Reproducible pipeline with fastp trimming via Nextflow:

```bash
nextflow run workflows/nextflow/main.nf -profile docker \
  --workflow amplicon_qiime2 \
  --manifest manifest.tsv --metadata metadata.tsv --classifier classifier.qza \
  --trim true
```
````

Note: the example manifest's `sampleB` SE row ends with a trailing tab (empty
`read2`); keep it — it documents the empty-`read2` convention.

- [ ] **Step 4: Add the four cross-links**

Add each as an additive line; do not remove existing content.

- `README.md` — after the "Three ways to use it" table (the `Nextflow workflows / CLI commands / Python SDK` table), add:
  ```markdown
  Processing many samples at once? See the [multisample & concurrency guide](docs/multisample.md).
  ```
- `docs/methods.md` — near the `trim` / `qc` entries, add:
  ```markdown
  > For batching many samples and choosing threads vs parallel jobs, see [multisample runs and concurrency](multisample.md).
  ```
- `docs/installation.md` — near the "Nextflow and HPC users" / "Docker users" rows, add:
  ```markdown
  > Running many samples? See [multisample runs and concurrency](multisample.md).
  ```
- `docs/api-nextflow.md` — where the manifest / `amplicon_qiime2` workflow is described, add:
  ```markdown
  > The manifest format is shared with `scripts/run_fastp_multiqc.sh`; see [multisample runs and concurrency](multisample.md).
  ```

- [ ] **Step 5: Run to verify pass**

Run: `uv run pytest tests/test_multisample_docs.py -v`
Expected: PASS (both tests).

- [ ] **Step 6: Commit**

```bash
git add docs/multisample.md README.md docs/methods.md docs/installation.md docs/api-nextflow.md tests/test_multisample_docs.py
git commit -m "docs(multisample): add multisample & concurrency guide with cross-links"
```

---

## Self-Review

**Spec coverage:**
- `docs/multisample.md` with four sections (manifest, concurrency, decision guide, examples) + example manifest + decision table → Task 1 Step 3. ✓
- Manifest schema + `jobs × threads` guidance match A and B → content copied from the spec's verified context. ✓
- Cross-links from README + methods + installation + api-nextflow → Task 1 Step 4. ✓
- Presence test (doc exists with section substrings; README links to it) → Task 1 Steps 1/5. ✓

**Placeholder scan:** none — full doc content, exact cross-link lines, and full test provided.

**Consistency:** the test's asserted substrings (`sample_id`, `jobs × threads`, `run_fastp_multiqc.sh`, `amplicon_qiime2`) all appear verbatim in the Step 3 doc content; the README link line contains `docs/multisample.md` which satisfies the `multisample.md` substring assertion.
