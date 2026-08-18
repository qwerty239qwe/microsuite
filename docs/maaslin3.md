# MaAsLin 3

MaAsLin 3 is available as an additional differential-abundance backend in
microsuite 0.4.0. It does not replace MaAsLin2; keep using MaAsLin2 when
reproducing an analysis that was fitted with that version.

## Basic use

```bash
microsuite diff_abundance \
  --backend maaslin3 \
  --table table.h5ad \
  --formula '~ batch + group / time_point + (1 | subject)' \
  --output maaslin3-results \
  --runtime docker
```

The Python API exposes the same call through `microsuite.api.diff_abundance`:

```python
from microsuite.api import diff_abundance

diff_abundance(
    backend="maaslin3",
    table=table_path,
    formula="~ batch + group / time_point + (1 | subject)",
    output=results_dir,
    runtime="docker",
)
```

`--formula` is passed to MaAsLin 3 as a complete lme4 formula and can contain
fixed effects, interactions, and random effects. As a convenience, the same
model can be split across `--fix-formula 'batch + group / time_point'` and
`--rand-formula '(1 | subject)'`. `--group` remains shorthand for a one-term
fixed-effects formula. A complete `--formula` cannot be combined with those
shorthands, so no term is silently ignored.

## Preprocessing controls

The backend exposes:

- `--normalization TSS|CLR|NONE` (default `TSS`)
- `--transform LOG|PLOG|NONE` (default `LOG`)
- `--min-prevalence` between 0 and 1 (default `0`)
- `--min-abundance` at least 0 (default `0`)

TSS normalization remains the default because the input is normally a raw
count table. Selecting `NONE` is appropriate only when the supplied values are
already on the intended analysis scale. MaAsLin 3 models prevalence as well as
non-zero abundance; include sequencing depth as a model covariate when it can
affect feature detection.

## Output contract

`--output` is a directory for this backend. It contains MaAsLin 3's native
`all_results.tsv` and `significant_results.tsv`, plus two stable microsuite
tables:

- `abundance_results.tsv`, containing only `model=abundance` rows
- `prevalence_results.tsv`, containing only `model=prevalence` rows

The directory also records `microsuite_params.json`; Docker runs add
`maaslin3_container.json`. Keeping the two model families separate prevents a
prevalence log-odds coefficient from being mistaken for an abundance fold
change. Use `--force` to replace an existing result directory; replacement is
performed only after the new run has produced both required tables.
