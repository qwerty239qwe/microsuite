# mothur reference data

The `mothur` backends (`microsuite cluster --backend mothur`,
`microsuite tax_classify --backend mothur`, and `microsuite workflow mothur`) do
not ship or download reference data. Every other microsuite backend that needs
a database — Kraken2, MetaPhlAn, EMU — has the same requirement, but mothur is
the one place it is most visible, because both clustering and classification
need a reference before the first command will run. This is a deliberate
design choice, not an oversight: see the "Reference data" decision in
[`docs/superpowers/specs/2026-07-25-mothur-workflow-design.md`](superpowers/specs/2026-07-25-mothur-workflow-design.md).

This page covers where to get that data and how to trim it down to a usable
size, plus a troubleshooting section for the failure modes that actually
occurred while building these backends.

## 1. The reference alignment (`--reference-alignment`)

`microsuite cluster --backend mothur` aligns your sequences against a curated
reference alignment before clustering (mothur's `align.seqs` step). mothur's
own MiSeq SOP page documents the standard source: the SILVA reference
alignment distributed from mothur.org.

```bash
curl -LO https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v138_1.tgz
tar -xzf silva.seed_v138_1.tgz
```

That archive is the full-length SILVA seed alignment, which is large and
covers the entire 16S gene. Aligning against the full-length reference is slow
and imprecise for an amplicon that only spans one variable region (V4, in a
typical MiSeq SOP run). Trim it down first with mothur's `pcr.seqs`, using
primer coordinates on the reference (515F/806R, i.e. columns ~11894-25319 in
`silva.seed_v138_1.align`):

```bash
mothur "#pcr.seqs(fasta=silva.seed_v138_1.align, start=11894, end=25319, keepdots=F, processors=8)"
```

This produces `silva.seed_v138_1.pcr.align` — pass that trimmed file as
`--reference-alignment`. If your primers target a different region, look up
the corresponding start/end coordinates for that region against the same
reference (the mothur wiki's MiSeq SOP page documents the common 16S primer
sets and their SILVA coordinates).

## 2. The taxonomy trainset (`--taxonomy-reference` / `--taxonomy-map`)

`microsuite tax_classify --backend mothur` runs mothur's naive Bayes
classifier (`classify.seqs`), which needs a matched pair of files rather than
a single prebuilt classifier artifact: a reference FASTA and a taxonomy file
mapping each reference ID to a lineage string. mothur.org distributes
ready-made trainsets for both RDP and SILVA taxonomies.

RDP trainset (v19):

```bash
curl -LO https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset19_072023.rdp.tgz
tar -xzf trainset19_072023.rdp.tgz
```

This unpacks `trainset19_072023.pds/trainset19_072023.pds.fasta` and
`trainset19_072023.pds/trainset19_072023.pds.tax` — pass those as
`--taxonomy-reference` and `--taxonomy-map` respectively.

SILVA taxonomy outline (v138) is distributed similarly under
`https://mothur.s3.us-east-2.amazonaws.com/wiki/` — see the mothur wiki's
"Taxonomy outline" page for the current archive name, since these files are
periodically re-published under version-specific names.

## 3. Worked example: `microsuite workflow mothur`

Given a directory of paired FASTQ files, the trimmed reference alignment from
step 1, and the trainset from step 2, run the full MiSeq SOP — `make.contigs`,
clustering, and classification — in one command. All five paths below are
required:

```bash
microsuite workflow mothur \
  --reads-dir reads/ \
  --reference-alignment silva.seed_v138_1.pcr.align \
  --taxonomy-reference trainset19_072023.pds/trainset19_072023.pds.fasta \
  --taxonomy-map trainset19_072023.pds/trainset19_072023.pds.tax \
  --out runs/mothur-sop
```

This writes `runs/mothur-sop/table.tsv` (feature-major OTU count table),
`runs/mothur-sop/table.shared` (mothur's native sample-major sidecar),
`runs/mothur-sop/rep-seqs.fasta` (one representative sequence per OTU),
`runs/mothur-sop/taxonomy.tsv` (per-OTU consensus taxonomy — see step 4), and
per-step logs under `runs/mothur-sop/logs/`.

`--identity` defaults to `0.97` and is converted internally to mothur's
distance-cutoff convention (`cutoff = 1 - identity`, so `0.97` becomes
`0.03`). Add `--force` to overwrite a previous run's outputs.

## 4. Worked example: standalone `cluster` + `tax_classify`

The workflow above composes two standalone backends. Running them directly
gives more control — for example, reusing an already-assembled contigs FASTA,
or classifying against a different trainset than the one used for the last
clustering run.

Cluster first. `--output-otu-list` and `--output-count-table` are optional,
but supplying them is what makes per-OTU consensus taxonomy possible in the
next step:

```bash
microsuite cluster --backend mothur \
  --rep-seqs contigs.fasta \
  --reference-alignment silva.seed_v138_1.pcr.align \
  --output-table table.tsv \
  --output-rep-seqs rep-seqs.fasta \
  --output-otu-list otu.list \
  --output-count-table table.count_table
```

Then classify. Feeding `cluster`'s `--output-otu-list` back in as
`--otu-list` (and its `--output-count-table` as `--count-table`) makes mothur
collapse per-sequence assignments into one consensus call per OTU via
`classify.otu`; without `--otu-list`, `tax_classify` returns per-sequence
taxonomy instead:

```bash
microsuite tax_classify --backend mothur \
  --rep-seqs rep-seqs.fasta \
  --taxonomy-reference trainset19_072023.pds/trainset19_072023.pds.fasta \
  --taxonomy-map trainset19_072023.pds/trainset19_072023.pds.tax \
  --otu-list otu.list \
  --count-table table.count_table \
  --output taxonomy.tsv
```

`--classifier` is rejected outright for `--backend mothur` — mothur's
classifier is the `--taxonomy-reference`/`--taxonomy-map` pair, not a single
prebuilt artifact, so passing `--classifier` here raises rather than silently
being ignored.

## Troubleshooting

**"contains '(),' which mothur's command syntax ... cannot represent"**

mothur's inline script syntax is `#command(key=value, key=value); command2(...)`.
Any path containing `(`, `)`, `,`, or `;` breaks that parse, and mothur has no
escaping mechanism for them. microsuite checks every path and parameter before
invoking mothur and raises early with the offending value named, rather than
handing mothur a command it will silently misparse.

This is most likely to bite on Windows, where the default 32-bit install
location is `C:\Program Files (x86)\...`. Put reference data, working
directories, and outputs somewhere without parentheses — for example
`C:\mothur-data\` instead of a path under `Program Files (x86)`.

**"Multiple files map to sample '...' mate 1/2"**

`microsuite workflow mothur` builds mothur's stability file (sample name, R1
path, R2 path) from your `--reads-dir`. mothur's `make.contigs` needs exactly
one R1 and one R2 per sample. If your sequencing run split reads across
multiple lanes (common straight off an Illumina bcl2fastq run, e.g.
`sampleA_S1_L001_R1_001.fastq.gz` and `sampleA_S1_L002_R1_001.fastq.gz`), both
lane files map to the same sample and mate, and this raises rather than
silently picking one lane and dropping reads from the other. Concatenate each
sample's lanes into a single R1 file and a single R2 file before running the
workflow:

```bash
cat sampleA_S1_L00*_R1_001.fastq.gz > sampleA_R1.fastq.gz
cat sampleA_S1_L00*_R2_001.fastq.gz > sampleA_R2.fastq.gz
```

**"Cannot determine sample/mate for FASTQ file '...'"**

If a file in `--reads-dir` has a `.fastq`/`.fq` (optionally `.gz`) extension
but its name doesn't match a recognized mate-pair pattern (`sampleA_R1.fastq.gz`,
`sampleA_S1_L001_R1_001.fastq.gz`, `sampleA.R1.fastq.gz`, etc.), this raises
naming the exact file rather than silently skipping it. A silently skipped
FASTQ would mean a sample goes missing from the run with no error — this
fails loudly instead so the file gets renamed or removed deliberately.

## Why this is user-supplied

Unlike Kraken2/Bracken/MetaPhlAn/EMU, mothur has no `refdb:<name>@<version>`
cached-reference shortcut in microsuite (see `--classifier` in
`microsuite tax_classify --help` for backends that do). Building one for
mothur was explicitly out of scope for this feature — see the "Reference
data" row in the Decisions table of
[`docs/superpowers/specs/2026-07-25-mothur-workflow-design.md`](superpowers/specs/2026-07-25-mothur-workflow-design.md) —
because it is a second subsystem to build and pin against mothur.org's URL
layout, which changes between releases independent of microsuite's release
cycle. If the archive names above have moved, check the current MiSeq SOP and
"Taxonomy outline" pages on the mothur wiki for the up-to-date links.
