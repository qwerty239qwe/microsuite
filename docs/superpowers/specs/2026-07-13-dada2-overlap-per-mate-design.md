# Per-mate DADA2 overlap precheck (Round-4 L) — Design

- **Date:** 2026-07-13
- **Status:** Approved (design), pending implementation
- **Origin:** Round-4 complaint **#7** (MEDIUM): the paired overlap precheck reads
  the length of one FASTQ and passes it as both `read_len_f` and `read_len_r`.
  Paired libraries can have unequal mate lengths (esp. after asymmetric primer
  trimming), so the overlap warning can be wrong. See [[microsuite-round4-roadmap]].

## Verified problem

`denoise_dada2_r` (denoise.py ~582-600) picks the first FASTQ in the input dir,
reads its length, and passes it as **both** `read_len_f` and `read_len_r` to
`dada2_qc.check_overlap`. `check_overlap` also ignores `trim_left`. So the
predicted overlap is wrong when R1/R2 differ in length or when `trimLeft` is set.

## Design

### 1. `dada2_qc.check_overlap` → structured, per-mate, trim-aware

Return an `OverlapReport` dataclass instead of `str | None`, and add `trim_left_f`
/ `trim_left_r`:

```python
@dataclass(frozen=True)
class OverlapReport:
    read_len_f: int
    read_len_r: int
    retained_f: int
    retained_r: int
    amplicon_length: int
    min_overlap: int
    predicted_overlap: int
    sufficient: bool
    warning: str | None

def check_overlap(*, trim_left_f, trim_left_r, trunc_len_f, trunc_len_r,
                  read_len_f, read_len_r, amplicon_length, min_overlap) -> OverlapReport:
    ...
```

- `retained_x = max((trunc_len_x if trunc_len_x > 0 else read_len_x) - trim_left_x, 0)`
  (DADA2 applies truncLen then trimLeft; the merge-usable length is what remains).
- `predicted_overlap = retained_f + retained_r - amplicon_length`.
- `sufficient = predicted_overlap >= min_overlap`; `warning` is the existing
  message (now citing both retained lengths) when not sufficient, else `None`.

### 2. Find one R1 and one matching R2 independently (`denoise.py`)

Add `_read_direction(stem) -> int | None` reusing `_READ_PATTERNS` (the `read`
group is `1`/`2`; the `forward|reverse` pattern maps forward→1, reverse→2), and
`_first_paired_reads(input_dir) -> tuple[Path | None, Path | None]` returning the
first R1 and first R2 FASTQ. The overlap wiring reads `read_len_f` from the R1 and
`read_len_r` from the R2 (each `first_read_length`), falling back to `0` when a
mate is missing. It passes `trim_left_f/r` + `trunc_len_f/r` from `tuning`.

The warn/raise behavior is unchanged: `report.warning` drives
`warnings.warn(...)` (or `MicrobiomeSuiteError` under `strict_qc`).

### 3. Record the overlap report in the provenance manifest

Stash the `OverlapReport` from the pre-run check and include it in the manifest's
`run` facts as an `overlap_check` block (read_len_f/r, retained_f/r,
amplicon_length, predicted_overlap, min_overlap, sufficient) — so
`dada2_denoise_manifest.json` records the forward length, reverse length, expected
amplicon length, and predicted overlap (#7's reporting requirement). Only added
when the check ran (paired + `amplicon_length` set).

## Testing (offline)

- `check_overlap`: equal mates (regression); **unequal** R1/R2 lengths →
  per-mate `retained_f`/`retained_r`, correct `predicted_overlap`; `trim_left`
  reduces retained length; `truncLen` overrides read length; sufficient vs
  insufficient → `warning` None vs message; `retained` clamped at 0.
- `_read_direction`/`_first_paired_reads`: `_R1/_R2`, `_1/_2`, `forward/reverse`
  conventions classify correctly; returns the first R1 and R2.
- `denoise_dada2_r` wiring (monkeypatched subprocess + a fixture dir with R1 and
  R2 FASTQs of **different** lengths): the overlap warning uses both mate lengths;
  `strict_qc` raises; and the written `dada2_denoise_manifest.json` has an
  `run.overlap_check` block with the per-mate lengths and predicted overlap.

## Success criteria

1. The overlap precheck reads R1 and R2 lengths independently and accounts for
   per-mate `trimLeft`/`truncLen`; unequal mates no longer give a wrong warning.
2. The QC/provenance manifest records forward length, reverse length, expected
   amplicon length, and predicted overlap.
3. Full offline suite green; `ty check`, `ruff check`, `ruff format --check` pass.

## Out of scope
- The other round-4 items (H-K, M). This is the isolated #7 fix.
