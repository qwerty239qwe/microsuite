# Nextflow FASTP Trim Module Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add per-sample fastp trimming to the `amplicon_qiime2` Nextflow workflow so trimming runs where the scheduler already parallelizes across samples, feeding trimmed reads into FASTQC/DADA2 and the fastp JSON into MultiQC.

**Architecture:** A new per-sample `FASTP` process emits trimmed reads + a fastp JSON. `main.nf` gates on `params.trim` (default true): when on, it synthesizes a trimmed manifest via `collectFile` and routes trimmed reads + that manifest to FASTQC and DADA2 (DADA2 unchanged — it matches reads by staged basename), and mixes the fastp JSON into MultiQC. `--trim false` reproduces today's raw-read behavior exactly.

**Tech Stack:** Nextflow DSL2 (Groovy), the existing `containers/fastp` GHCR image, pytest string/file-presence tests (the repo does not run Nextflow in CI), plus one opt-in `-stub` smoke test gated on `MICROSUITE_RUN_EXTERNAL_INTEGRATION=1` + `nextflow` on PATH.

## Global Constraints

- Only the `amplicon_qiime2` workflow is touched; `amplicon_microsuite` is untouched.
- Non-breaking: `--trim false` must reproduce today's raw-read behavior (FASTQC + DADA2 on raw reads).
- No new container image — reuse the existing `containers/fastp` (`${params.container_registry}/fastp:${params.container_tag}`).
- Every new Nextflow process ships both a `script:` block and a `stub:` block.
- fastp threads sized modestly (`params.fastp_cpus`, default 4) so the scheduler runs many samples concurrently, not one sample at max threads.
- Nextflow tests are string/file-presence checks in `tests/test_nextflow_skeleton.py` (constants `ROOT`, `NEXTFLOW` already defined there). New python test files start with `from __future__ import annotations`.
- The existing guard assertion `"MULTIQC.out.report_dir," not in main` in `test_nextflow_skeleton.py` must keep passing.

---

### Task 1: FASTP module + params + container wiring

**Files:**
- Create: `workflows/nextflow/modules/fastp.nf`
- Modify: `workflows/nextflow/nextflow.config` (params block)
- Modify: `workflows/nextflow/profiles/docker.config` (add `withLabel: fastp`)
- Modify: `workflows/nextflow/profiles/singularity.config` (add `withLabel: fastp`)
- Test: `tests/test_nextflow_skeleton.py` (append two tests)

**Interfaces:**
- Consumes: nothing.
- Produces: process `FASTP` with input `tuple val(sample_id), path(reads)` and outputs `emit: trimmed` (`tuple val(sample_id), path("${sample_id}*.trim.fastq.gz")`) and `emit: report` (`path "${sample_id}.fastp.json"`). Params `trim` (bool, default true), `fastp_cpus` (int, default 4), `fastp_args` (string, default '').

- [ ] **Step 1: Write the failing tests** (append to `tests/test_nextflow_skeleton.py`)

```python
def test_nextflow_fastp_module_exists_and_declares_process() -> None:
    fastp = NEXTFLOW / "modules" / "fastp.nf"
    assert fastp.exists(), fastp
    text = fastp.read_text(encoding="utf-8")
    assert "process FASTP" in text
    assert "label 'fastp'" in text
    assert "params.fastp_cpus" in text
    assert "emit: trimmed" in text
    assert "emit: report" in text
    # PE and SE fastp invocations
    assert "--in1" in text and "--out1" in text
    assert "--in2" in text and "--out2" in text
    assert "params.fastp_args" in text
    assert "stub:" in text


def test_nextflow_fastp_params_and_container_labels() -> None:
    config = (NEXTFLOW / "nextflow.config").read_text(encoding="utf-8")
    for key in ("trim", "fastp_cpus", "fastp_args"):
        assert key in config, key
    docker = (NEXTFLOW / "profiles" / "docker.config").read_text(encoding="utf-8")
    singularity = (NEXTFLOW / "profiles" / "singularity.config").read_text(encoding="utf-8")
    assert "withLabel: fastp" in docker
    assert "withLabel: fastp" in singularity
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_nextflow_skeleton.py::test_nextflow_fastp_module_exists_and_declares_process tests/test_nextflow_skeleton.py::test_nextflow_fastp_params_and_container_labels -v`
Expected: FAIL (module file missing / assertions unmet).

- [ ] **Step 3: Create `workflows/nextflow/modules/fastp.nf`**

```groovy
process FASTP {
    tag "${sample_id}"
    label 'fastp'
    cpus { params.fastp_cpus as int }
    publishDir "${params.outdir}/trim/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}*.trim.fastq.gz"), emit: trimmed
    path "${sample_id}.fastp.json", emit: report
    path "${sample_id}.fastp.html"

    script:
    def paired = reads instanceof List && reads.size() > 1
    if (paired)
        """
        fastp --in1 ${reads[0]} --in2 ${reads[1]} \
              --out1 ${sample_id}_1.trim.fastq.gz --out2 ${sample_id}_2.trim.fastq.gz \
              --json ${sample_id}.fastp.json --html ${sample_id}.fastp.html \
              --thread ${task.cpus} ${params.fastp_args}
        """
    else
        """
        fastp --in1 ${reads instanceof List ? reads[0] : reads} \
              --out1 ${sample_id}.trim.fastq.gz \
              --json ${sample_id}.fastp.json --html ${sample_id}.fastp.html \
              --thread ${task.cpus} ${params.fastp_args}
        """

    stub:
    def paired = reads instanceof List && reads.size() > 1
    if (paired)
        """
        printf '' | gzip > ${sample_id}_1.trim.fastq.gz
        printf '' | gzip > ${sample_id}_2.trim.fastq.gz
        printf '{"summary":{}}\\n' > ${sample_id}.fastp.json
        printf '<html><body>stub fastp</body></html>\\n' > ${sample_id}.fastp.html
        """
    else
        """
        printf '' | gzip > ${sample_id}.trim.fastq.gz
        printf '{"summary":{}}\\n' > ${sample_id}.fastp.json
        printf '<html><body>stub fastp</body></html>\\n' > ${sample_id}.fastp.html
        """
}
```

- [ ] **Step 4: Add params to `workflows/nextflow/nextflow.config`**

In the `params { ... }` block, after the existing `sampling_depth = 1000` line, add:
```groovy
    trim = true
    fastp_cpus = 4
    fastp_args = ''
```

- [ ] **Step 5: Add the container label to both profiles**

In `workflows/nextflow/profiles/docker.config`, inside `process { ... }`, add:
```groovy
    withLabel: fastp {
        container = "${params.container_registry}/fastp:${params.container_tag}"
    }
```
In `workflows/nextflow/profiles/singularity.config`, inside `process { ... }`, add:
```groovy
    withLabel: fastp {
        container = "docker://${params.container_registry}/fastp:${params.container_tag}"
    }
```

- [ ] **Step 6: Run to verify pass**

Run: `uv run pytest tests/test_nextflow_skeleton.py -v`
Expected: PASS (all skeleton tests, including the two new ones).

- [ ] **Step 7: Commit**

```bash
git add workflows/nextflow/modules/fastp.nf workflows/nextflow/nextflow.config workflows/nextflow/profiles/docker.config workflows/nextflow/profiles/singularity.config tests/test_nextflow_skeleton.py
git commit -m "feat(nextflow): add FASTP trim module + params + container labels"
```

---

### Task 2: Wire FASTP into the `amplicon_qiime2` workflow

**Files:**
- Modify: `workflows/nextflow/main.nf` (add include; replace the FASTQC/MULTIQC/DADA2 block at lines ~112–116)
- Test: `tests/test_nextflow_skeleton.py` (append one test)

**Interfaces:**
- Consumes: process `FASTP` (Task 1) with `FASTP.out.trimmed` (`tuple(sample_id, [reads])`) and `FASTP.out.report` (fastp JSON path).
- Produces: no new symbols; changes the workflow DAG so FASTQC/DADA2 consume trimmed reads + a `trimmed_manifest.tsv` when `params.trim` is true.

- [ ] **Step 1: Write the failing test** (append to `tests/test_nextflow_skeleton.py`)

```python
def test_nextflow_main_wires_fastp_trim() -> None:
    main = (NEXTFLOW / "main.nf").read_text(encoding="utf-8")
    assert "include { FASTP } from './modules/fastp'" in main
    assert "FASTP(samples_ch)" in main
    assert "params.trim" in main
    assert "collectFile" in main
    assert "trimmed_manifest.tsv" in main
    assert ".mix(extra_qc)" in main
    # raw (non-breaking) path preserved in the else branch
    assert "dada2_manifest = manifest_ch" in main
    assert "dada2_reads  = reads_ch" in main or "dada2_reads = reads_ch" in main
    # FASTQC and DADA2 each still invoked exactly once, via the toggled inputs
    assert main.count("FASTQC(") == 1
    assert "QIIME2_DADA2(dada2_manifest, metadata_ch, dada2_reads)" in main
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_nextflow_skeleton.py::test_nextflow_main_wires_fastp_trim -v`
Expected: FAIL (wiring absent).

- [ ] **Step 3: Add the include near the other module includes** (top of `main.nf`, alongside `include { FASTQC } ...`)

```groovy
include { FASTP } from './modules/fastp'
```

- [ ] **Step 4: Replace the FASTQC/MULTIQC/DADA2 block**

Replace these current lines in the `amplicon_qiime2` workflow body:
```groovy
    reads_ch = samples_ch.map { sample_id, reads -> reads }.flatten().collect()

    FASTQC(samples_ch)
    MULTIQC(FASTQC.out.qc_dir.collect())
    QIIME2_DADA2(manifest_ch, metadata_ch, reads_ch)
```
with:
```groovy
    reads_ch = samples_ch.map { sample_id, reads -> reads }.flatten().collect()

    if (params.trim) {
        FASTP(samples_ch)
        trimmed_ch = FASTP.out.trimmed
        dada2_manifest = trimmed_ch
            .map { sid, reads ->
                def r2 = reads.size() > 1 ? reads[1].name : ''
                "${sid}\t${reads[0].name}\t${r2}\n"
            }
            .collectFile(name: 'trimmed_manifest.tsv', seed: 'sample_id\tread1\tread2\n', sort: true)
        fastqc_input = trimmed_ch
        dada2_reads  = trimmed_ch.map { sid, reads -> reads }.flatten().collect()
        extra_qc     = FASTP.out.report
    } else {
        fastqc_input = samples_ch
        dada2_manifest = manifest_ch
        dada2_reads  = reads_ch
        extra_qc     = Channel.empty()
    }

    FASTQC(fastqc_input)
    MULTIQC(FASTQC.out.qc_dir.mix(extra_qc).collect())
    QIIME2_DADA2(dada2_manifest, metadata_ch, dada2_reads)
```

(The `QIIME2_TAXONOMY`, `QIIME2_PHYLOGENY`, `QIIME2_DIVERSITY`, and `REPORT`
calls that follow are unchanged.)

- [ ] **Step 5: Run to verify pass + no skeleton regression**

Run: `uv run pytest tests/test_nextflow_skeleton.py -v`
Expected: PASS (including the existing `test_nextflow_amplicon_modules_are_declared`, which still finds every previously-declared module and the `"MULTIQC.out.report_dir," not in main` guard).

- [ ] **Step 6: Commit**

```bash
git add workflows/nextflow/main.nf tests/test_nextflow_skeleton.py
git commit -m "feat(nextflow): run FASTP per sample and feed trimmed reads to FASTQC/DADA2/MultiQC"
```

---

### Task 3: Opt-in `-stub` smoke test

**Files:**
- Create: `tests/integration/test_nextflow_fastp_stub.py`

**Interfaces:**
- Consumes: the wired `main.nf` (Tasks 1–2), a `nextflow` binary, no network.
- Produces: an opt-in end-to-end DAG check proving the stubbed workflow completes with `--trim true` and publishes fastp outputs.

- [ ] **Step 1: Write the test**

```python
# tests/integration/test_nextflow_fastp_stub.py
from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
NEXTFLOW = ROOT / "workflows" / "nextflow"

pytestmark = pytest.mark.skipif(
    os.environ.get("MICROSUITE_RUN_EXTERNAL_INTEGRATION") != "1",
    reason="set MICROSUITE_RUN_EXTERNAL_INTEGRATION=1 to run external-tool integration tests",
)


@pytest.mark.parametrize("layout", ["single-end", "paired-end"])
def test_nextflow_stub_runs_amplicon_qiime2_with_trim(layout: str, tmp_path: Path) -> None:
    if shutil.which("nextflow") is None:
        pytest.skip("nextflow is not installed on PATH")

    r1 = tmp_path / "s1_R1.fastq.gz"
    r1.write_bytes(gzip.compress(b"@r\nACGT\n+\nIIII\n"))
    if layout == "paired-end":
        r2 = tmp_path / "s1_R2.fastq.gz"
        r2.write_bytes(gzip.compress(b"@r\nACGT\n+\nIIII\n"))
        row = f"s1\t{r1}\t{r2}\n"
    else:
        # single-end: empty read2 -> FASTP emits ONE trimmed file -> FASTP.out.trimmed
        # is a scalar Path (not a list). This is the exact case that regressed.
        row = f"s1\t{r1}\t\n"

    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(f"sample_id\tread1\tread2\n{row}", encoding="utf-8")
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("sample-id\tgroup\ns1\ta\n", encoding="utf-8")
    classifier = tmp_path / "classifier.qza"
    classifier.write_text("stub", encoding="utf-8")
    outdir = tmp_path / "results"

    cmd = [
        "nextflow", "run", str(NEXTFLOW / "main.nf"),
        "-stub-run", "-profile", "local",
        "--workflow", "amplicon_qiime2",
        "--manifest", str(manifest),
        "--metadata", str(metadata),
        "--classifier", str(classifier),
        "--outdir", str(outdir),
        "--trim", "true",
    ]
    result = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True, timeout=900)
    assert result.returncode == 0, result.stderr
    assert (outdir / "trim" / "fastp").exists(), "fastp outputs were not published"
```

Both layouts run under `--trim true`; the single-end case specifically exercises
the scalar-`Path` output of `FASTP.out.trimmed` that the trimmed-manifest closure
must normalize.

- [ ] **Step 2: Verify it skips by default**

Run: `uv run pytest tests/integration/test_nextflow_fastp_stub.py -q`
Expected: 2 skipped (env var unset; parametrized SE + PE).

- [ ] **Step 3: (If `nextflow` is available) verify it passes once**

Run: `MICROSUITE_RUN_EXTERNAL_INTEGRATION=1 uv run pytest tests/integration/test_nextflow_fastp_stub.py -q`
Expected: 2 passed if `nextflow` is on PATH; otherwise 2 SKIP with the "nextflow is not installed" reason. Record which occurred. The single-end case is the regression guard.

- [ ] **Step 4: Commit**

```bash
git add tests/integration/test_nextflow_fastp_stub.py
git commit -m "test(nextflow): opt-in -stub smoke test for FASTP trim path"
```

---

## Self-Review

**Spec coverage:**
- `modules/fastp.nf` with script+stub, PE+SE, `cpus=params.fastp_cpus`, `fastp_args` passthrough → Task 1. ✓
- `main.nf` runs FASTP when `params.trim` (default true), derives trimmed manifest via `collectFile`, feeds trimmed reads+manifest to FASTQC/DADA2, mixes fastp JSON into MultiQC → Task 2. ✓
- `--trim false` reproduces raw behavior → Task 2 else-branch + test asserting `dada2_manifest = manifest_ch`. ✓
- `withLabel: fastp` in docker + singularity; no new image → Task 1. ✓
- params `trim`/`fastp_cpus`/`fastp_args` → Task 1. ✓
- skeleton tests assert all of it; existing guard preserved → Tasks 1–2. ✓
- opt-in `-stub` smoke test, skips cleanly without nextflow → Task 3. ✓

**Placeholder scan:** none — every step has concrete code/commands.

**Type/interface consistency:** `FASTP.out.trimmed` (tuple sid,[reads]) and `FASTP.out.report` (json) named in Task 1 are exactly what Task 2 consumes; `params.trim`/`fastp_cpus`/`fastp_args` names match across Tasks 1–2; `dada2_manifest`/`dada2_reads`/`fastqc_input`/`extra_qc` are internal to the Task 2 block and used consistently.
