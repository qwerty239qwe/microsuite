# Nextflow `FASTP` Trim Module — Design

- **Date:** 2026-07-06
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** A codex-authored `run_fastp_multiqc.sh` performed batched fastp
  trimming by looping the single-sample `microsuite trim` CLI and summarizing
  with MultiQC. Analysis showed the repo's reproducible Nextflow workflow
  (`amplicon_qiime2`) already reads the identical manifest TSV and does
  per-sample FASTQC → MultiQC in parallel, but has **no fastp/trim step** — so
  the script re-implemented, sequentially and non-portably, work the scheduler
  already does, for a step the pipeline lacks. This sub-project (B of three)
  adds fastp trimming where the parallelism already lives.

## Background: the concurrency model (why this design)

- **fastp / cutadapt / trimmomatic / trim-galore** are multithreaded *within a
  single sample* and process **one sample per invocation**. fastp's thread
  scaling is sublinear on amplicon-sized files (little gain past ~4–8 threads).
- **FastQC** `--threads N` is file-parallel within one call; **MultiQC** is a
  single-process aggregation.
- **`microsuite trim`/`qc`** have **no cross-sample parallelism** — one CLI call
  = one sample, tool threads within. There is no Python-level multiprocessing in
  `src/`.
- **Nextflow** is the only place with true cross-sample parallelism: `samples_ch`
  fans out one process per sample, scheduled concurrently by the executor.

Design consequence: put trimming in a per-sample Nextflow process and let the
scheduler parallelize it; size fastp threads modestly so *many samples* run
concurrently rather than one sample hogging all cores.

## Scope

Sub-project **B** of three (A = portable standalone script; C = multisample
API/CLI/workflow docs + ergonomics). This spec covers **B only**.

### Out of scope for B
- The standalone `run_fastp_multiqc.sh` (sub-project A).
- Documentation of the shared manifest schema / concurrency guidance (C).
- Any change to the `amplicon_microsuite` workflow variant (this touches only
  `amplicon_qiime2`).
- New container images (the `containers/fastp` image already builds and
  publishes to GHCR).

## Current integration facts (verified)

- `main.nf` `amplicon_qiime2` builds `samples_ch` from the manifest TSV
  (`sample_id/read1/read2`) as `tuple(sample_id, [read1, read2?])`, then calls
  `FASTQC(samples_ch)`, `MULTIQC(FASTQC.out.qc_dir.collect())`, and
  `QIIME2_DADA2(manifest_ch, metadata_ch, reads_ch)` where
  `reads_ch = samples_ch.map{ _,reads -> reads }.flatten().collect()`.
- `QIIME2_DADA2` rebuilds its QIIME manifest from the microsuite manifest using
  read **basenames staged in the process work dir** (`gsub` strips the path,
  points at `ENVIRON["PWD"]/basename`). So it consumes whatever read files are
  staged by name — it does not care whether they are raw or trimmed, as long as
  the manifest basenames match the staged files.
- Modules get containers via `withLabel: <name>` in `profiles/docker.config`
  and `profiles/singularity.config`, mapping to
  `${params.container_registry}/<image>:${params.container_tag}`.
- Nextflow is not run in CI; `tests/test_nextflow_skeleton.py` asserts file
  existence + string presence in `main.nf`/config. Every module ships a `stub:`
  block for offline `nextflow run -stub`.

## Design

### Architecture: derive a trimmed manifest, leave DADA2 unchanged

Insert a per-sample `FASTP` process. To keep all change isolated to the new
module + one conditional block, FASTP emits `tuple(sample_id, [trimmed_reads])`
plus a `fastp.json`; the workflow synthesizes a **trimmed manifest** (whose
read1/read2 are the trimmed basenames) via Nextflow's `collectFile`, and routes
the **trimmed reads + trimmed manifest** to FASTQC and DADA2. DADA2 works
unchanged because it matches reads by staged basename. MultiQC receives FASTQC
dirs **mixed** with the fastp JSON (MultiQC parses fastp natively).

Gated by `params.trim` (default `true`). With `--trim false` the workflow is
byte-for-byte today's behavior (raw reads → FASTQC/DADA2), so the change is
non-breaking; DADA2's own `trim_left`/`trunc_len` truncation still applies on
top of (or instead of) fastp.

**Rejected alternative:** refactor `QIIME2_DADA2` to take a
`tuple(sample_id, reads)` channel instead of manifest+reads. Cleaner long-term
but modifies a working module and its assumptions; deferred.

### Component 1 — `modules/fastp.nf`

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

`cpus = params.fastp_cpus` (default **4**, not `auto`) so the scheduler runs
many samples concurrently. `params.fastp_args` is a passthrough for
adapter/quality options.

### Component 2 — `main.nf` wiring (amplicon_qiime2)

Add `include { FASTP } from './modules/fastp'`. Replace the unconditional
FASTQC/MULTIQC/DADA2 calls with:

```groovy
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
QIIME2_TAXONOMY(QIIME2_DADA2.out.rep_seqs, classifier_ch)
// ... rest unchanged
```

`FASTQC` and `QIIME2_DADA2` are each still called exactly once; only their
inputs are swapped by the toggle.

### Component 3 — params + container config

`nextflow.config` params block gains:
```groovy
trim = true
fastp_cpus = 4
fastp_args = ''
```
`profiles/docker.config` and `profiles/singularity.config` each gain:
```groovy
withLabel: fastp {
    container = "${params.container_registry}/fastp:${params.container_tag}"
}
```

## Testing

Extend `tests/test_nextflow_skeleton.py` (string/file-presence, matching the
existing convention):
- `modules/fastp.nf` exists;
- `main.nf` contains `include { FASTP }`, a `FASTP(` call, `params.trim`,
  `collectFile`, and the `.mix(` of the fastp report into MultiQC;
- `nextflow.config` contains `trim`, `fastp_cpus`, `fastp_args`;
- `withLabel: fastp` present in BOTH `docker.config` and `singularity.config`;
- the existing assertion `MULTIQC.out.report_dir," not in main` (guarding a
  known past bug) still holds.

Add an **opt-in** offline smoke test (gated like the repo's external-integration
tests, e.g. `MICROSUITE_RUN_EXTERNAL_INTEGRATION=1` + `nextflow` on PATH) that
runs `nextflow run workflows/nextflow/main.nf -stub -profile local --trim true`
against a tiny fixture manifest and asserts the DAG completes and a
`trimmed_manifest.tsv`/fastp outputs are produced. Skips cleanly when `nextflow`
is absent.

## Success criteria

1. `modules/fastp.nf` exists with `script:` and `stub:` blocks, PE and SE paths.
2. `main.nf` runs FASTP per sample when `params.trim` (default true), feeds
   trimmed reads + a derived trimmed manifest to FASTQC and DADA2, and mixes the
   fastp JSON into MultiQC.
3. `--trim false` reproduces today's raw-read behavior exactly (non-breaking).
4. `withLabel: fastp` wired in docker + singularity profiles; no new image.
5. `tests/test_nextflow_skeleton.py` asserts all of the above; the full offline
   suite stays green.
6. `nextflow run -stub` executes the whole DAG offline (verified via the opt-in
   smoke test where `nextflow` is available).

## Open questions / follow-ups (not blocking B)

- Whether DADA2 should eventually consume a `tuple(sample_id, reads)` channel
  directly (removes the trimmed-manifest synthesis) — a later refactor.
- `fastp_cpus`/`fastp_args` defaults may want tuning once run on real data
  (sub-project C's guidance).
