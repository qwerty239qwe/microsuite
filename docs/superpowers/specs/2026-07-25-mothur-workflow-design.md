# mothur Workflow — Design

- **Date:** 2026-07-25
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** microsuite covers DADA2/Deblur (ASV) and VSEARCH/USEARCH (OTU)
  feature-table generation, but has no mothur path. mothur remains a widely used
  16S pipeline with its own alignment-based OTU workflow (the MiSeq SOP) and its
  own naive-Bayes classifier, neither of which is reachable today.

## Background & problem

`microsuite cluster` currently offers `vsearch`, `usearch`, and `qiime2-vsearch`
(`docs/methods.md:87-89`). All three are identity-threshold clustering over
unaligned sequences. mothur's OptiClust path differs materially: it aligns
sequences against a curated reference alignment, screens on alignment
coordinates, denoises with `pre.cluster`, removes chimeras, and only then
clusters on a distance matrix. Users migrating from mothur, or reproducing a
published mothur analysis, cannot do so in microsuite.

`microsuite tax_classify` similarly offers `qiime2`, `kraken2`, `bracken`,
`metaphlan`, and `emu` (`docs/methods.md:108-113`) but not mothur's
`classify.seqs`/`classify.otu` consensus assignment, which is what mothur-based
publications report.

## Scope

Two backends plus one end-to-end workflow composed from them.

| # | Deliverable | Summary |
|---|---|---|
| 1 | `methods/mothur.py` | The `_run_mothur` primitive: invoke one mothur command, parse its `Output File Names:` block, return outputs keyed by extension. |
| 2 | `cluster --backend mothur` | FASTA in, OTU table + representative FASTA out, via the clustering half of the MiSeq SOP. |
| 3 | `tax_classify --backend mothur` | `classify.seqs` → `classify.otu` consensus taxonomy. |
| 4 | `microsuite workflow mothur` | FASTQ directory → `make.contigs` → deliverable 2 → deliverable 3. |
| 5 | `containers/mothur/Dockerfile` | bioconda mothur image. Enables a real run for capturing parser fixtures and for users without a local mothur. No CI job — see [Deferred](#deferred). |

### Out of scope

- A mothur reference-database provider in `refdb`. Reference data is
  user-supplied (see [Reference data](#reference-data)).
- `qc_filter --backend mothur`. mothur's `screen.seqs` operates on an aligned
  FASTA + `.count_table` pair; exposing it as a standalone verb would leak SOP
  intermediates into the CLI for no standalone use case.
- mothur's phylogeny, diversity, and `.biom` export commands. microsuite already
  covers those verbs natively; mothur's versions add no capability.
- `--runtime docker` orchestration beyond what `containers/mothur` provides.

## Decisions

Each of these was an explicit fork, recorded so the plan does not relitigate them.

| Decision | Choice | Rejected alternative and why |
|---|---|---|
| Shape | Backends first, workflow composed on top | A workflow-only implementation would bury mothur's steps where nothing else can reuse them. |
| Verb mapping | `cluster` + `tax_classify` | `qc_filter` also — rejected, see Out of scope. |
| Execution | One mothur process per SOP step | A single batch script is faster and shorter, but collapses ~10 steps into one opaque log entry, so a mid-SOP failure names no step. Per-step invocation is what buys per-step provenance. |
| Filenames | Parse mothur's `Output File Names:` stdout block | Deriving names from mothur's tag rules (`.trim.contigs.good.unique...`) means re-encoding a convention that changes between mothur releases. Globbing a per-step dir is ambiguous when one command emits two files sharing an extension. |
| Reference data | User-supplied paths | A `refdb` provider is better UX but is a second subsystem to build and pin against mothur.org URL churn. Matches the existing kraken2/metaphlan/tax4fun2 convention. |
| Runtime | Container image, external binary also supported. CI smoke deferred. | Shipping the image without a CI job is the cheap half: it still makes a real run reproducible, which is how the parser fixtures get captured. Automated drift detection is the part deliberately given up. |
| `cluster` input | FASTA, like `vsearch`/`usearch` | A FASTQ directory is more faithful to the SOP but gives `cluster` a different input type than its two siblings. `make.contigs` lives in the workflow instead. |
| Code structure | Thin primitive + straight-line SOP | A declarative step list turns per-step quirks into dataclass fields, so understanding one step means reading two places. One function per step yields ~10 permanently single-caller wrappers. |

## Architecture

```
methods/mothur.py          _run_mothur() + stdout parser   [sole owner of "how to talk to mothur"]
        │
        ├── methods/cluster.py        cluster_mothur()
        ├── methods/tax_classify.py   tax_classify_mothur()
        │
workflows/mothur_sop.py    make.contigs, then calls the two backends
containers/mothur/         bioconda image + CI smoke fixture
```

The parser is the load-bearing component and has exactly one home. When a mothur
release changes its stdout format, one test file fails, not ten.

### Component 1 — `methods/mothur.py`

```python
def _run_mothur(
    command: str,
    params: dict[str, str],
    *,
    work_dir: Path,
    run_dir: Path | None,
    timeout: float | None,
) -> dict[str, Path]:
```

Builds `mothur "#set.dir(output=<work_dir>); <command>(<params>)"`, executes it
through the existing `runtime.runner` machinery so each step lands in the run-dir
log, and returns the command's outputs keyed by file extension — for example
`{"fasta": ..., "count_table": ...}`.

`set.dir(output=...)` is required, not cosmetic: without it mothur writes its
`mothur.<timestamp>.logfile` and all intermediates into the process working
directory.

Backend resolution follows the existing convention: `shutil.which("mothur")`,
raising `MicrobiomeSuiteError` with an install hint when absent, mirroring
`_require_qiime` in `methods/denoise.py:1061`.

### Component 2 — `cluster --backend mothur`

Straight-line, one `_run_mothur` call per arrow:

```
seqs.fasta
  → unique.seqs        → .fasta + .count_table
  → align.seqs         → .align                  (user reference alignment)
  → screen.seqs        → .fasta + .count_table   (start/end/maxambig/maxhomop)
  → filter.seqs        → .fasta
  → unique.seqs        → .fasta + .count_table
  → pre.cluster        → .fasta + .count_table   (diffs)
  → chimera.vsearch    → .accnos                 (dereplicate=t)
  → remove.seqs        → .fasta + .count_table
  → dist.seqs          → .dist                   (cutoff)
  → cluster            → .list                   (method=opti)
  → make.shared        → .shared
  → get.oturep         → .fasta                  (representative sequences)
```

`.shared` is sample-major (rows = samples); microsuite's table contract is
feature-major (`feature-id<TAB>sample…`, see
`methods/cluster.py:264 write_otu_table_from_uc`). The final step transposes it.

**Outputs** match the existing `cluster` contract: table TSV at `--output-table`,
representative FASTA at `--output-rep-seqs`, and the raw `.shared` written beside
the output table as `<output-table stem>.shared` — the same sidecar role `.uc`
plays for the vsearch and usearch backends.

`--reference-alignment` is a required option for this backend.

The existing `--identity` option carries over and is converted once, at the top
of `cluster_mothur`, to the distance cutoff both `dist.seqs` and `cluster`
take: `cutoff = round(1.0 - identity, 4)`, so the microsuite default of `0.97`
becomes mothur's conventional `0.03`. The same cutoff value is passed to both
commands; there is no separate distance-cutoff option.

mothur-specific screening parameters (`--maxambig`, `--maxhomop`,
`--pre-cluster-diffs`) are exposed with mothur's own defaults and rejected for
the other `cluster` backends via the `_reject_options` pattern.

### Component 3 — `tax_classify --backend mothur`

`classify.seqs(fasta, count, reference=<user fasta>, taxonomy=<user tax>)` →
`.taxonomy`, then `classify.otu(list, taxonomy, count)` → `.cons.taxonomy` when a
`--list` from the cluster step is supplied. Without a list file, the per-sequence
`.taxonomy` is the output.

`--classifier` is not reused here. mothur's classifier is not a prebuilt
artifact but a *pair* of files — a reference FASTA and a matching tax file — so
`tax_classify` gains two new options:

| Option | Meaning | Backends that accept it |
|---|---|---|
| `--taxonomy-reference` | Reference sequence FASTA | `mothur` only |
| `--taxonomy-map` | Reference taxonomy file (`id<TAB>lineage`) | `mothur` only |
| `--classifier` | Prebuilt classifier artifact or DB directory | every backend **except** `mothur` |

Both directions of mismatch **raise** `MicrobiomeSuiteError` via the
`_reject_options` pattern in `methods/trim.py:381` — passing `--classifier` to
`mothur`, or either new option to a non-mothur backend.

Raising rather than warning-and-ignoring is deliberate. A silently ignored
reference path does not fail; it classifies every sequence against whatever
database the backend fell back to and returns a full, well-formed, *wrong*
taxonomy table. That is indistinguishable from a correct result downstream, so
it must be caught at the argument boundary.

### Component 4 — `microsuite workflow mothur`

A `WorkflowSpec` entry in `workflows/catalog.py` plus `workflows/mothur_sop.py`:

1. Write mothur's stability file from the FASTQ directory (sample name + R1 + R2
   per line).
2. `make.contigs` → contigs FASTA.
3. Call `cluster(backend="mothur", ...)`.
4. Call `tax_classify(backend="mothur", ...)`.
5. Emit `run.json` through the existing metadata/stage machinery.

The workflow is an orchestrator; it contains no mothur command construction of
its own beyond `make.contigs`.

## Error handling

**Corrected 2026-07-26 against real mothur 1.48.5.** This section originally
asserted that mothur exits 0 on failure. That is **false** for 1.48.5: a failed
`align.seqs` (missing reference) returns **exit code 1**, verified twice. The
design premise was wrong and the consequences are recorded here rather than
quietly patched.

What follows from the measurement:

- `run_command` already raises on non-zero exit, so the primary failure path is
  the exit code — not a stdout scan.
- `check_mothur_errors` is therefore **never reached** on an exit-1 failure. It
  is retained as defence in depth, because exactly one failure mode was
  sampled and mothur has historically been reported to continue past some
  non-fatal command errors. It guards a hypothetical, and the spec should not
  pretend otherwise.
- The error scan must anchor on `"[ERROR]: "` (colon and space), **not** the
  bare substring `[ERROR]`. mothur closes a failed run with a summary banner —
  `Detected 1 [ERROR] messages, please review.` — which a bare-substring scan
  double-counts, turning one failure into two reported errors.
- Because `run_command` raises with the whole captured stream, an exit-1 mothur
  failure would otherwise surface as ~60 lines of citation banner and startup
  noise. `run_mothur` catches that failure and re-raises with only the anchored
  `[ERROR]: ` lines, so the user sees `did not complete align.seqs.` rather
  than mothur's front matter.

The remaining four failure modes:

| Condition | Behaviour |
|---|---|
| `Output File Names:` block empty or absent | Raise, naming the step that produced nothing. |
| Requested extension missing from the block | Raise, listing the extensions actually produced. |
| Reference alignment path missing or unreadable | Validated with `ensure_input` **before** step 1, so a bad path fails in a second rather than six steps in. |
| A step removes every sequence | mothur writes an empty FASTA and continues. Detect and raise, naming the step — the same principle as the `filtered_sample_mask` guard in `methods/r/dada2_denoise.R`. |

Per-step stdout is captured into the run directory by the existing
`_run`/`CommandLog` path, which is the provenance benefit the per-step execution
model was chosen for.

## Testing

| Level | Coverage |
|---|---|
| Parser unit tests | Recorded real-mothur stdout: a normal block, exit-0-with-`[ERROR]`, an empty block, a missing extension. |
| Command construction | Per step, with monkeypatched `subprocess.run`, mirroring `test_denoise_qiime2_deblur_builds_command` (`tests/test_denoise_cluster_methods.py:366`). |
| Option rejection | `--classifier` against `mothur`, and each new taxonomy option against a non-mothur backend, both raise. |
| Transpose unit test | `.shared` → feature-major TSV, including a sample with zero counts. |

The parser fixtures are **captured, not written**: run the SOP once inside
`containers/mothur` and save the real stdout under `tests/fixtures/mothur/`. A
hand-invented sample would test the parser against our assumption of mothur's
format rather than the format itself, which is precisely the bug class the
parser exists to survive.

This suite runs entirely without mothur installed, so it belongs in the default
`pytest` run. Validation level for `docs/methods.md` is therefore
**unit-tested wrapper + user environment**, matching the usearch and picrust2
rows — not the CI-smoke-tested rows.

## Deferred

**CI smoke test.** A real end-to-end SOP run in CI is the only check that
catches mothur changing its stdout format between releases. It is deferred
because it requires a committed fixture that cannot be synthesised: a few
hundred reads plus a trimmed SILVA slice that is a *genuine* alignment, since
`align.seqs` rejects unaligned input.

Consequence accepted: format drift surfaces in user environments rather than in
CI. Add it when either a mothur upgrade breaks the parser in the wild, or the
mothur backends move from user-supplied references toward a `refdb` provider
(which would supply the fixture problem's solution as a side effect).

## Documentation

- `docs/methods.md`: rows in *Feature Table Generation*, *Denoising And
  Clustering Backends*, and *Taxonomy And Phylogeny*; a validation-status row for
  the mothur family.
- `docs/api-cli.md`: the new options.
- A `docs/mothur.md` covering reference-data acquisition from mothur.org, since
  user-supplied references are the one place this backend is less turnkey than
  its siblings.
