# Table I/O ergonomics (Round-3 D) — Design

- **Date:** 2026-07-11
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** Codex's ERP120510 friction list, round-3 complaints **#12** (CLR only
  ergonomic for `.h5ad`, not matrix TSVs), **#13** (no direct h5ad→TSV export),
  **#14** (already-collapsed profile matrices are awkward inputs — transforms),
  **#15** (importing collapsed rank tables hits AnnData naming conflicts). First
  sub-project (**D**) of round-3; **E** (assignment QC #16/#17) and **F**
  (metadata-aware viz #18) follow. See [[dada2-improvement-roadmap]].

## Scope

D makes the four table conversions ergonomic and consistent:

| # | Conversion | Command | Status |
|---|------------|---------|--------|
| 1 | tsv + metadata → h5ad | `import tsv` | exists; gains #15 sanitize |
| 2 | h5ad → h5ad (normalized matrix added as a **layer**, X preserved) | `normalize` | behavior change |
| 3 | tsv → tsv (transform a bare matrix) | `table normalize` | new |
| 4 | h5ad → tsv (+ optional metadata) | `table export` | new |

It also fixes the #15 first-column naming clash in the shared TSV reader, which
benefits `import tsv` and the new tsv-native `table` group alike.

### Out of scope for D
- Plotting on profile matrices (the "plots" half of #14) — overlaps E/F viz work
  (#16 assignment plots, #18 metadata-aware viz); deferred there.
- New transform methods beyond `normalize_native`'s
  (`relative`, `total-sum`, `clr`, `prevalence-filter`).
- Non-native normalize backends (only `native` supports these methods).

## Verified context

- `normalize()` (`methods/normalize.py`) is h5ad-in/h5ad-out and today **replaces
  X** with the transformed matrix. Core `normalize_native(adata, *, method,
  target_sum, pseudocount, min_prevalence) -> AnnData` supports `relative`,
  `total-sum`, `clr` (all **shape-preserving**) and `prevalence-filter` (which
  **removes** low-prevalence features → different shape).
- `abundance` writes a TSV from an AnnData via `dense_counts(adata)` +
  `DataFrame.to_csv(sep="\t")`; `obs_names` = samples, `var_names` = features.
- `io/tsv.py:read_tsv(table, metadata, taxonomy)` requires metadata, sets the
  first column as the feature index (its NAME leaks into `var.index.name`), then
  `add_taxonomy_levels(var)` adds `var` columns `kingdom…species`
  (`io/taxonomy.py:LEVELS`). A first column named a rank (e.g. `genus`) collides
  with the `genus` `var` column on `write_h5ad` — the #15 break.
- `read_h5ad`/`write_h5ad` in `io/h5ad.py`; CLI groups registered in `cli/app.py`
  via `app.add_typer(<mod>.app, name=...)`. Fatal → `MicrobiomeSuiteError`
  (`microsuite._errors`); non-fatal → `warnings.warn` (as in sub-projects A/B).

## Design

### Component 1 — shared matrix reader + #15 sanitize (`io/tsv.py`)

```python
_RESERVED_FEATURE_NAMES = set(LEVELS) | {"taxonomy", "taxon"}

def read_count_matrix(path: Path) -> pd.DataFrame:
    """Read a features×samples count matrix TSV (first column = feature IDs):
    numeric-coerce sample columns, reject empty/duplicate IDs, normalize the
    feature index name to 'feature_id'. If the original first-column name is a
    reserved rank/taxonomy name, warn that it was renamed (IDs preserved)."""

def read_matrix_tsv(path: Path) -> ad.AnnData:
    """Build an AnnData from a bare matrix — no metadata. obs is an empty frame
    indexed by the sample columns; var carries add_taxonomy_levels."""
```

`read_count_matrix`: (1) `pd.read_csv(sep="\t")`; reject empty or `<2` cols; (2)
capture `original = str(table.columns[0])`, set it as index, cast to str, reject
duplicate IDs; (3) numeric-coerce the rest (`to_numeric(errors="raise")`), cast
column labels to str; (4) `counts.index.name = "feature_id"`; if
`original.strip().lower() in _RESERVED_FEATURE_NAMES`, `warnings.warn("Renamed
feature-ID column '<original>' to 'feature_id' to avoid a taxonomy-rank naming
conflict; feature IDs are unchanged.")`.

`read_matrix_tsv`: `counts = read_count_matrix(path)`;
`obs = pd.DataFrame(index=counts.columns)`;
`var = add_taxonomy_levels(pd.DataFrame(index=counts.index))`;
`AnnData(X=counts.T.to_numpy(float64), obs=obs, var=var)`, `uns["microsuite"]`
stamp with `importer="matrix-tsv"`.

Refactor `read_tsv` to call `read_count_matrix` for the parse/sanitize, keeping
its metadata join and taxonomy handling — so `import tsv` inherits the #15 fix.

### Component 2 — layer-mode normalize (`methods/normalize.py`)

`normalize()` (h5ad→h5ad) changes so shape-preserving transforms **preserve raw
counts in X and write the result into `layers[<method>]`**:

- `SHAPE_PRESERVING = {"relative", "total-sum", "clr"}`.
- For a shape-preserving method: read the input adata, compute
  `result = normalize_native(adata, method=...)`, then set
  `adata.layers[<layer_name>] = result.X` (X untouched) and write `adata`. The
  layer name is the method (`clr`, `relative`, `total-sum`).
- For `prevalence-filter` (structural — removes features): it cannot be a layer of
  the original X; keep today's behavior (write the filtered adata, X = filtered
  counts). This exception is documented in the command help.
- Re-normalizing an already-normalized h5ad reads the raw X, so multiple layers
  (e.g. `clr` and `relative`) can coexist.

`normalize_native` itself is unchanged (still returns a transformed AnnData); only
the `normalize()` wrapper's write step changes.

### Component 3 — `table` methods (`methods/table_io.py`)

```python
def export_table(*, table: Path, output: Path, layer: str | None = None,
                 metadata: Path | None = None, force: bool = False) -> None:
    """Read an .h5ad; write its X (or layers[layer]) as a features×samples TSV
    (index 'feature_id' = var_names, columns = obs_names). If metadata is given,
    also write obs to that TSV path."""

def normalize_table(*, method: str, input_path: Path, output: Path,
                    target_sum: float = 1_000_000.0, pseudocount: float = 1.0,
                    min_prevalence: float = 0.1, force: bool = False) -> None:
    """Read a matrix TSV, apply normalize_native(method=...), write the
    transformed matrix as a features×samples TSV."""
```

- `export_table`: `adata = read_h5ad(ensure_input(table))`; pick the matrix
  (`adata.layers[layer]` if `layer` else `dense_counts(adata)`), raising
  `MicrobiomeSuiteError` if a requested `layer` is absent (message lists available
  layers); build `pd.DataFrame(matrix.T, index=var_names, columns=obs_names)` with
  `index.name = "feature_id"`; `to_csv(prepare_output(output, force), sep="\t")`.
  If `metadata` is set, `adata.obs` → `to_csv(prepare_output(metadata, force),
  sep="\t")` (index named `sample`).
- `normalize_table`: `adata = read_matrix_tsv(ensure_input(input_path))`;
  `result = normalize_native(adata, method=method, ...)`; write `result` (X) via
  the same features×samples `to_csv` shape (shared inner helper with
  `export_table`). For `prevalence-filter` the written matrix has the surviving
  features (tsv → tsv, single matrix — no layers involved).

### Component 4 — CLI (`cli/table_cmd.py`)

`app = typer.Typer(help="Transform and export feature/profile tables as TSV.",
no_args_is_help=True)`, registered in `cli/app.py` as
`app.add_typer(table_cmd.app, name="table")`.

- `table export`: `--table PATH` (h5ad, required), `--output/-o PATH` (required),
  `--layer TEXT` (optional; default exports X), `--metadata PATH` (optional; also
  dump obs), `--force`. Calls `export_table`.
- `table normalize`: `--method TEXT` (required; `relative`/`total-sum`/`clr`/
  `prevalence-filter`), `--input PATH` (matrix TSV, required), `--output/-o PATH`
  (required), `--target-sum`, `--pseudocount`, `--min-prevalence`, `--force`.
  Calls `normalize_table`.

The existing `normalize` command (`method_tables_cmd.py`) keeps its options; only
its output semantics change per Component 2 (help text updated to say the result
is stored in a layer, X preserved; prevalence-filter filters features).

### Data flow

- (1) `counts.tsv + metadata.tsv → import tsv → x.h5ad`.
- (2) `x.h5ad → normalize --method clr → x_clr.h5ad` (X = counts, `layers["clr"]`
  = CLR).
- (3) `counts.tsv → table normalize --method clr → clr.tsv` (no metadata).
- (4) `x.h5ad → table export [--layer clr] [--metadata meta.tsv] → matrix.tsv`.

## Testing (offline)

- `read_count_matrix`: rank-named first column (`genus`) → `pytest.warns` and
  `index.name == "feature_id"`; normal first column → no warn; empty / duplicate
  IDs → `MicrobiomeSuiteError`.
- `read_matrix_tsv`: builds an AnnData with the right shape and empty obs, no
  metadata arg.
- `import tsv` regression: first column named `genus` now imports without error
  and warns (the exact #15 case).
- `normalize` layer mode: `--method clr` on an h5ad → output `X` equals the
  original counts AND `layers["clr"]` equals `normalize_native(clr).X`;
  `--method relative` likewise adds `layers["relative"]`; `prevalence-filter`
  writes a filtered adata (fewer vars, no layer).
- `export_table`: import a fixture → `table export` → TSV matches original counts
  (features×samples round-trip); `--layer clr` after a layer-normalize exports the
  CLR matrix; `--metadata` writes obs; a missing `--layer` name →
  `MicrobiomeSuiteError` naming available layers.
- `normalize_table`: `--method clr` on a `counts.tsv` equals
  `normalize_native(clr)`; `--method relative` sums to 1 per sample; `--force`.
- CLI smoke (typer `CliRunner`): `table export` / `table normalize` wire options
  through.

## Success criteria

1. `import tsv` handles the four→h5ad path and no longer crashes on a
   rank-named first column (warns instead); feature IDs preserved.
2. `normalize --method clr` on an h5ad writes h5ad with raw counts in `X` and the
   CLR matrix in `layers["clr"]` (X preserved); other shape-preserving methods add
   their own layer; `prevalence-filter` filters features.
3. `table normalize --method clr --input counts.tsv -o clr.tsv` does tsv→tsv CLR
   with no metadata and no h5ad glue; all native methods work.
4. `table export --table x.h5ad -o x.tsv` writes a features×samples TSV, with
   optional `--layer` selection and optional `--metadata` obs dump.
5. Full offline suite green and both CI gates pass
   (`ruff check .`, `ruff format --check .`).

## Open questions / follow-ups (not blocking D)

- `table export` long-format / multi-layer output could be added later.
- A configurable layer name for `normalize` (beyond the method name) is deferred.
- Plotting on profile matrices (#14 second half) is deferred to E/F.
