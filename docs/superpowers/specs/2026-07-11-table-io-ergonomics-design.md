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

D adds a tsv-native `table` command group so downstream workflows can transform
and move profile matrices without the import→h5ad→export glue codex kept writing:
`table export` (h5ad→TSV) and `table normalize` (TSV→TSV, all native methods incl
CLR). It also fixes the #15 first-column naming clash in the shared TSV reader,
which benefits the existing `import tsv` too.

### Out of scope for D
- Plotting on profile matrices (the "plots" half of #14) — overlaps E/F viz work
  (#16 assignment plots, #18 metadata-aware viz); deferred there.
- New transform methods beyond what `normalize_native` already supports
  (relative, total-sum, clr, prevalence-filter).
- Non-native normalize backends (only `native` supports these methods).
- Layer/obs-column selection in `table export` — it exports the primary `X`
  matrix as-is.

## Verified context

- `normalize()` (`methods/normalize.py`) is h5ad-in/h5ad-out; the reusable core
  `normalize_native(adata, *, method, target_sum, pseudocount, min_prevalence)
  -> AnnData` supports `relative`, `total-sum`, `clr`, `prevalence-filter`.
- `abundance` already writes a TSV from an AnnData via `dense_counts(adata)` +
  `DataFrame.to_csv(sep="\t")`; `obs_names` = samples, `var_names` = features.
- `io/tsv.py:read_tsv(table, metadata, taxonomy)` **requires** metadata, sets the
  first column as the feature index (so its NAME leaks into `var.index.name`),
  then `add_taxonomy_levels(var)` unconditionally adds `var` columns
  `kingdom…species` (`io/taxonomy.py:LEVELS`). When the first column is named a
  rank (e.g. `genus`), `var.index.name == "genus"` collides with the `genus`
  `var` column on `write_h5ad` — the #15 break.
- `read_h5ad`/`write_h5ad` exist in `io/h5ad.py`. CLI groups are registered in
  `cli/app.py` via `app.add_typer(<mod>.app, name=...)`.
- Methods/readers raise `MicrobiomeSuiteError` (`microsuite._errors`) for fatal
  cases; the codebase uses `warnings.warn` for non-fatal notices (established by
  sub-projects B and A).

## Design

### Component 1 — shared matrix reader (`io/tsv.py`)

Extract the count-matrix parsing into a reusable helper and add a metadata-free
variant:

```python
_RESERVED_FEATURE_NAMES = set(LEVELS) | {"taxonomy", "taxon"}

def read_count_matrix(path: Path) -> pd.DataFrame:
    """Read a features×samples count matrix TSV (first column = feature IDs).
    Numeric-coerce sample columns, reject empties/duplicates. Normalize the
    feature index name to 'feature_id'; if the original first-column name is a
    reserved rank/taxonomy name, warn that it was renamed (IDs are preserved)."""

def read_matrix_tsv(path: Path) -> ad.AnnData:
    """Build an AnnData from a bare count matrix — no metadata required. obs is an
    empty frame indexed by the sample columns; var carries add_taxonomy_levels."""
```

`read_count_matrix`:
1. `table = pd.read_csv(path, sep="\t")`; reject if empty or `< 2` columns
   (`MicrobiomeSuiteError`, same message as today).
2. Capture `original = str(table.columns[0])`; set that column as index; cast
   index to str; reject duplicate feature IDs.
3. Numeric-coerce the remaining columns (`apply(pd.to_numeric, errors="raise")`),
   cast column labels to str.
4. Set `counts.index.name = "feature_id"`. If
   `original.strip().lower() in _RESERVED_FEATURE_NAMES`,
   `warnings.warn(f"Renamed feature-ID column '{original}' to 'feature_id' to
   avoid a taxonomy-rank naming conflict; feature IDs are unchanged.")`.
5. Return the features×samples frame.

`read_matrix_tsv`: `counts = read_count_matrix(path)`;
`obs = pd.DataFrame(index=counts.columns)`;
`var = add_taxonomy_levels(pd.DataFrame(index=counts.index))`;
`AnnData(X=counts.T.to_numpy(float64), obs=obs, var=var)` with a `uns["microsuite"]`
provenance stamp mirroring `read_tsv` (`importer="matrix-tsv"`).

Refactor `read_tsv` to call `read_count_matrix` for steps 1–4 (so `import tsv`
inherits the #15 sanitize), keeping its metadata join and taxonomy handling.

### Component 2 — `table` methods (`methods/table_io.py`)

```python
def export_table(*, table: Path, output: Path, force: bool = False) -> None:
    """Read an .h5ad and write its X matrix as a features×samples TSV
    (index 'feature_id' = var_names, columns = obs_names)."""

def normalize_table(*, method: str, input_path: Path, output: Path,
                    target_sum: float = 1_000_000.0, pseudocount: float = 1.0,
                    min_prevalence: float = 0.1, force: bool = False) -> None:
    """Read a matrix TSV, apply normalize_native(method=...), write a matrix TSV."""
```

`export_table`: `adata = read_h5ad(ensure_input(table))`;
`frame = pd.DataFrame(dense_counts(adata).T, index=adata.var_names.astype(str),
columns=adata.obs_names.astype(str))`; `frame.index.name = "feature_id"`;
`frame.to_csv(prepare_output(output, force=force), sep="\t")`.

`normalize_table`: `adata = read_matrix_tsv(ensure_input(input_path))`;
`result = normalize_native(adata, method=method, ...)`; write via the same
features×samples `to_csv` shape as `export_table` (shared inner helper).

### Component 3 — CLI (`cli/table_cmd.py`)

`app = typer.Typer(help="Transform and export feature/profile tables as TSV.",
no_args_is_help=True)`, registered in `cli/app.py` as
`app.add_typer(table_cmd.app, name="table")`.

- `table export`: `--table PATH` (h5ad, required), `--output/-o PATH` (required),
  `--force`. Calls `export_table`.
- `table normalize`: `--method TEXT` (required; `relative`/`total-sum`/`clr`/
  `prevalence-filter`), `--input PATH` (matrix TSV, required), `--output/-o PATH`
  (required), `--target-sum`, `--pseudocount`, `--min-prevalence`, `--force`.
  Calls `normalize_table`.

### Data flow

`counts.tsv → read_matrix_tsv → AnnData → normalize_native(clr) → to_csv → clr.tsv`
(one `table normalize` call). `x.h5ad → read_h5ad → to_csv → x.tsv` (one
`table export` call). No metadata anywhere in the `table` group.

## Testing (offline)

- `read_count_matrix`: a rank-named first column (`genus`) → warns
  (`pytest.warns`) and returns the matrix with `index.name == "feature_id"`; a
  normal first column (`feature_id`/`#OTU ID`) → no warn; empty and duplicate-ID
  inputs → `MicrobiomeSuiteError`.
- `read_matrix_tsv`: builds an AnnData with the right shape and empty obs; no
  metadata argument needed.
- `import tsv` regression: a table whose first column is named `genus` now
  imports to `.h5ad` without error and emits the rename warning (the exact #15
  case that previously broke).
- `export_table`: import a fixture (`import tsv`) → `table export` → the TSV
  matches the original counts (round-trip on values, features×samples).
- `normalize_table`: `--method clr` on a small `counts.tsv` yields values equal to
  `normalize_native(clr)`; `--method relative` sums to 1 per sample; `--force`
  overwrite behavior.
- CLI smoke (typer `CliRunner`): `table export` and `table normalize` wire the
  options through to the methods.

## Success criteria

1. `microsuite table export --table x.h5ad -o x.tsv` writes a features×samples
   TSV of the h5ad's `X` matrix.
2. `microsuite table normalize --method clr --input counts.tsv -o clr.tsv`
   produces CLR (and the other native methods) with no import/export glue and no
   metadata.
3. A matrix whose first column is named a taxonomic rank reads and imports with a
   warning and no crash; feature IDs are preserved.
4. The `table` group requires no metadata; existing `import tsv` behavior is
   unchanged except it now also warns-and-sanitizes the rank-named first column.
5. The full offline suite stays green and both CI gates pass
   (`ruff check .`, `ruff format --check .`).

## Open questions / follow-ups (not blocking D)

- `table export` could later grow `--layer`/`--obs` selection or long-format
  output; deferred until a workflow needs it.
- A `table import`/round-trip that also carries metadata/taxonomy is already
  served by `import tsv`; the `table` group intentionally stays metadata-free.
- Plotting on profile matrices (#14 second half) is deferred to E/F.
