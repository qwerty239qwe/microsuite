# Table I/O Ergonomics (Round-3 D) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the four table conversions ergonomic: `import tsv` (tsv+metadata→h5ad, now #15-safe), `normalize` (h5ad→h5ad with the transform stored in a layer, X preserved), `table normalize` (tsv→tsv), `table export` (h5ad→tsv + optional metadata).

**Architecture:** A shared `read_count_matrix` in `io/tsv.py` parses+sanitizes matrices for both `read_tsv` (metadata path) and a new metadata-free `read_matrix_tsv`. `normalize()` writes shape-preserving results into `layers[<method>]`. A new tsv-native `table` CLI group (`export`, `normalize`) reuses `read_h5ad`/`read_matrix_tsv`/`normalize_native`.

**Tech Stack:** Python 3.12, pandas, anndata, numpy, Typer, pytest (`CliRunner`, `pytest.warns`).

## Global Constraints

- Fatal → `MicrobiomeSuiteError` (`microsuite._errors`); non-fatal → `warnings.warn`.
- Feature matrices are **features×samples** on disk (first column = feature IDs, header row = sample IDs). In AnnData, `obs_names` = samples, `var_names` = features; `dense_counts(adata)` is samples×features.
- `#15 sanitize`: the feature index name is always normalized to `feature_id`; warn only when the original first-column name (case-insensitive, stripped) is in `set(LEVELS) | {"taxonomy", "taxon"}`.
- `normalize` layer mode applies to `SHAPE_PRESERVING = {"relative", "total-sum", "clr"}` (preserve X, write `layers[method]`); `prevalence-filter` stays a structural feature filter (writes the filtered adata).
- Native backend only for normalize methods. `NORMALIZE_METHODS = ("relative", "total-sum", "clr", "prevalence-filter")`.
- Both CI gates must pass: `uv run ruff check .` and `uv run ruff format --check .`.
- `from __future__ import annotations` at the top of new modules.

---

### Task 1: shared matrix reader + #15 sanitize (`io/tsv.py`)

**Files:**
- Modify: `src/microsuite/io/tsv.py`
- Test: `tests/test_io_tsv_matrix.py` (new)

**Interfaces:**
- Produces:
  - `read_count_matrix(path: Path) -> pd.DataFrame` (features×samples, index name `feature_id`)
  - `read_matrix_tsv(path: Path) -> ad.AnnData`
  - `read_tsv` unchanged signature, now reusing `read_count_matrix`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_io_tsv_matrix.py
from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.io.tsv import read_count_matrix, read_matrix_tsv


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_read_count_matrix_basic(tmp_path: Path) -> None:
    p = _write(tmp_path / "c.tsv", "feature_id\ts1\ts2\nASV1\t5\t1\nASV2\t0\t3\n")
    m = read_count_matrix(p)
    assert list(m.index) == ["ASV1", "ASV2"]
    assert list(m.columns) == ["s1", "s2"]
    assert m.index.name == "feature_id"
    assert m.loc["ASV1", "s2"] == 1


def test_read_count_matrix_sanitizes_rank_first_column(tmp_path: Path) -> None:
    p = _write(tmp_path / "g.tsv", "genus\ts1\ts2\nBacteroides\t5\t1\nPrevotella\t0\t3\n")
    with pytest.warns(UserWarning, match="feature_id"):
        m = read_count_matrix(p)
    assert m.index.name == "feature_id"
    assert list(m.index) == ["Bacteroides", "Prevotella"]  # IDs preserved


def test_read_count_matrix_no_warn_normal_header(tmp_path: Path, recwarn) -> None:
    p = _write(tmp_path / "c.tsv", "#OTU ID\ts1\nASV1\t5\n")
    read_count_matrix(p)
    assert not any(issubclass(w.category, UserWarning) for w in recwarn.list)


def test_read_count_matrix_rejects_empty_and_dups(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        read_count_matrix(_write(tmp_path / "e.tsv", "feature_id\n"))
    with pytest.raises(MicrobiomeSuiteError):
        read_count_matrix(_write(tmp_path / "d.tsv", "feature_id\ts1\nASV1\t5\nASV1\t3\n"))


def test_read_matrix_tsv_no_metadata(tmp_path: Path) -> None:
    p = _write(tmp_path / "c.tsv", "feature_id\ts1\ts2\nASV1\t5\t1\nASV2\t0\t3\n")
    adata = read_matrix_tsv(p)
    assert adata.shape == (2, 2)  # 2 samples (obs) x 2 features (var)
    assert list(adata.obs_names) == ["s1", "s2"]
    assert list(adata.var_names) == ["ASV1", "ASV2"]
    assert adata.obs.shape[1] == 0  # empty metadata
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_io_tsv_matrix.py -v`
Expected: FAIL (`ImportError: cannot import name 'read_count_matrix'`).

- [ ] **Step 3: Refactor `io/tsv.py`**

The current `read_tsv` inlines the parse. Extract it into `read_count_matrix`, add `read_matrix_tsv`, and make `read_tsv` reuse the helper. Add `import warnings` and `from microsuite.io.taxonomy import LEVELS` to the imports (keep the existing `add_taxonomy_levels`, `normalize_taxonomy_columns`, `read_indexed_tsv`, `__version__`, `datetime` imports).

```python
_RESERVED_FEATURE_NAMES = set(LEVELS) | {"taxonomy", "taxon"}


def read_count_matrix(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    if table.empty or table.shape[1] < 2:
        raise MicrobiomeSuiteError(
            "Feature table must have a feature ID column and sample columns."
        )
    original = str(table.columns[0])
    counts = table.set_index(table.columns[0])
    counts.index = counts.index.astype(str)
    if counts.index.has_duplicates:
        raise MicrobiomeSuiteError("Feature table contains duplicate feature IDs.")
    counts = counts.apply(pd.to_numeric, errors="raise")
    counts.columns = counts.columns.astype(str)
    counts.index.name = "feature_id"
    if original.strip().lower() in _RESERVED_FEATURE_NAMES:
        warnings.warn(
            f"Renamed feature-ID column '{original}' to 'feature_id' to avoid a "
            "taxonomy-rank naming conflict; feature IDs are unchanged.",
            stacklevel=2,
        )
    return counts


def read_matrix_tsv(path: Path) -> ad.AnnData:
    counts = read_count_matrix(path)
    obs = pd.DataFrame(index=counts.columns)
    var = add_taxonomy_levels(pd.DataFrame(index=counts.index))
    adata = ad.AnnData(X=counts.T.to_numpy(dtype=np.float64), obs=obs, var=var)
    adata.uns["microsuite"] = {
        "version": __version__,
        "importer": "matrix-tsv",
        "table": str(path),
        "created_at": datetime.now(UTC).isoformat(),
    }
    return adata
```

Then rewrite `read_tsv` to reuse it (replace its inline parse block):

```python
def read_tsv(
    table_path: Path, metadata_path: Path, taxonomy_path: Path | None = None
) -> ad.AnnData:
    counts = read_count_matrix(table_path)

    metadata = read_indexed_tsv(metadata_path, index_name="sample")
    missing = [sample for sample in counts.columns if sample not in metadata.index]
    if missing:
        raise MicrobiomeSuiteError(f"Metadata is missing samples from table: {missing[:5]}")
    metadata = metadata.loc[counts.columns].copy()

    var = pd.DataFrame(index=counts.index)
    if taxonomy_path is not None:
        taxonomy = read_indexed_tsv(taxonomy_path, index_name="feature")
        taxonomy = taxonomy.rename(columns=normalize_taxonomy_columns(taxonomy.columns))
        taxonomy = taxonomy.reindex(var.index)
        var = var.join(taxonomy)
    var = add_taxonomy_levels(var)

    adata = ad.AnnData(
        X=counts.T.to_numpy(dtype=np.float64),
        obs=metadata,
        var=var,
    )
    adata.uns["microsuite"] = {
        "version": __version__,
        "importer": "tsv",
        "table": str(table_path),
        "metadata": str(metadata_path),
        "taxonomy": str(taxonomy_path) if taxonomy_path else None,
        "created_at": datetime.now(UTC).isoformat(),
    }
    return adata
```

(Note: `var.index` now carries name `feature_id`, so `add_taxonomy_levels` no longer clashes with a rank-named first column — the #15 fix.)

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_io_tsv_matrix.py -v`
Expected: PASS (5).

- [ ] **Step 5: Add the `import tsv` #15 regression test**

Append to `tests/test_io_tsv_matrix.py`:

```python
def test_read_tsv_with_rank_named_first_column(tmp_path: Path) -> None:
    from microsuite.io.h5ad import write_h5ad
    from microsuite.io.tsv import read_tsv

    _write(tmp_path / "g.tsv", "genus\ts1\ts2\nBacteroides\t5\t1\nPrevotella\t0\t3\n")
    _write(tmp_path / "m.tsv", "sample\tgroup\ns1\ta\ns2\tb\n")
    with pytest.warns(UserWarning, match="feature_id"):
        adata = read_tsv(tmp_path / "g.tsv", tmp_path / "m.tsv")
    # the previously-breaking write must now succeed
    write_h5ad(adata, tmp_path / "out.h5ad")
    assert (tmp_path / "out.h5ad").exists()
    assert list(adata.var_names) == ["Bacteroides", "Prevotella"]
```

Run: `uv run pytest tests/test_io_tsv_matrix.py -v` → PASS (6). Then confirm no regression in the existing importer/CLI tests: `uv run pytest tests/test_cli.py -q`.

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/io/tsv.py tests/test_io_tsv_matrix.py
git commit -m "feat(io): shared read_count_matrix + metadata-free read_matrix_tsv + #15 sanitize"
```

---

### Task 2: layer-mode `normalize` (`methods/normalize.py`)

**Files:**
- Modify: `src/microsuite/methods/normalize.py`
- Test: `tests/test_table_ecology_methods.py` (append; update any existing normalize assertions that read X)

**Interfaces:**
- Consumes: `normalize_native` (unchanged).
- Produces: `normalize()` writes shape-preserving results to `layers[method]`, X preserved; `SHAPE_PRESERVING` set exported from the module.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_table_ecology_methods.py`)

```python
def test_normalize_clr_writes_layer_preserves_x(tmp_path) -> None:
    import numpy as np

    from microsuite.io.h5ad import read_h5ad, write_h5ad
    from microsuite.io.tsv import read_matrix_tsv
    from microsuite.methods.normalize import normalize, normalize_native

    counts = tmp_path / "counts.tsv"
    counts.write_text("feature_id\ts1\ts2\nA\t5\t1\nB\t3\t9\n", encoding="utf-8")

    adata_in = read_matrix_tsv(counts)
    src = tmp_path / "in.h5ad"
    write_h5ad(adata_in, src)

    out = tmp_path / "out.h5ad"
    normalize(backend="native", method="clr", table=src, output=out)
    result = read_h5ad(out)
    # X still raw counts
    assert np.allclose(result.X, adata_in.X)
    # clr stored as a layer, matching normalize_native
    assert "clr" in result.layers
    expected = normalize_native(adata_in, method="clr").X
    assert np.allclose(result.layers["clr"], expected)


def test_normalize_prevalence_filter_filters_features(tmp_path) -> None:
    from microsuite.io.h5ad import read_h5ad, write_h5ad
    from microsuite.io.tsv import read_matrix_tsv
    from microsuite.methods.normalize import normalize

    counts = tmp_path / "counts.tsv"
    # feature B present in only 1 of 2 samples -> prevalence 0.5
    counts.write_text("feature_id\ts1\ts2\nA\t5\t1\nB\t0\t9\n", encoding="utf-8")
    src = tmp_path / "in.h5ad"
    write_h5ad(read_matrix_tsv(counts), src)
    out = tmp_path / "out.h5ad"
    normalize(backend="native", method="prevalence-filter", table=src, output=out, min_prevalence=0.75)
    result = read_h5ad(out)
    assert list(result.var_names) == ["A"]  # B filtered
    assert "prevalence-filter" not in result.layers  # structural, not a layer
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_table_ecology_methods.py -k "normalize_clr_writes_layer or prevalence_filter_filters" -v`
Expected: FAIL (today `normalize` replaces X; no `clr` layer).

- [ ] **Step 3: Change the `normalize()` wrapper**

Add `SHAPE_PRESERVING = ("relative", "total-sum", "clr")` near the constants, and rewrite the write step of `normalize()`:

```python
def normalize(
    *,
    backend: str,
    method: str,
    table: Path,
    output: Path,
    target_sum: float = 1_000_000.0,
    pseudocount: float = 1.0,
    min_prevalence: float = 0.1,
    force: bool = False,
) -> None:
    backend = backend.lower()
    if backend != "native":
        raise MicrobiomeSuiteError(
            f"Unsupported normalize backend '{backend}'. "
            f"Choose one of: {', '.join(SUPPORTED_BACKENDS)}"
        )
    source = read_h5ad(ensure_input(table))
    result = normalize_native(
        source,
        method=method,
        target_sum=target_sum,
        pseudocount=pseudocount,
        min_prevalence=min_prevalence,
    )
    normalized_method = method.lower()
    if normalized_method in SHAPE_PRESERVING:
        source.layers[normalized_method] = result.X
        source.uns["microsuite_normalize"] = result.uns["microsuite_normalize"]
        out_adata = source
    else:
        out_adata = result
    write_h5ad(out_adata, prepare_output(output, force=force))
```

`normalize_native` is unchanged.

- [ ] **Step 4: Run to verify pass + update any broken existing tests**

Run: `uv run pytest tests/test_table_ecology_methods.py -v`
Expected: the two new tests PASS. If any pre-existing `normalize` test asserted the transformed values on `result.X` for a shape-preserving method, update it to read `result.layers[method]` (X now holds raw counts). Do NOT change assertions for `prevalence-filter` (still on X/var). Run the full method test module to confirm: `uv run pytest tests/test_table_ecology_methods.py tests/test_methods.py -q`.

- [ ] **Step 5: Update the `normalize` CLI help text**

In `cli/method_tables_cmd.py` `normalize_cmd`, update the `--method` help to note the result is stored in a layer named after the method with raw counts preserved in X (and that `prevalence-filter` removes features). No behavior change in the CLI wiring itself.

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/methods/normalize.py src/microsuite/cli/method_tables_cmd.py tests/test_table_ecology_methods.py
git commit -m "feat(normalize): store shape-preserving transforms in layers, preserve raw X"
```

---

### Task 3: `table` group — export + tsv normalize (`methods/table_io.py`, `cli/table_cmd.py`)

**Files:**
- Create: `src/microsuite/methods/table_io.py`
- Create: `src/microsuite/cli/table_cmd.py`
- Modify: `src/microsuite/cli/app.py` (register the group)
- Test: `tests/test_table_io.py` (new)

**Interfaces:**
- Consumes: `read_h5ad`, `read_matrix_tsv` (Task 1), `normalize_native`, `dense_counts`, `ensure_input`, `prepare_output`.
- Produces: `export_table(*, table, output, layer=None, metadata=None, force=False)`, `normalize_table(*, method, input_path, output, target_sum=..., pseudocount=..., min_prevalence=..., force=False)`; CLI `table export`, `table normalize`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_table_io.py
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.io.h5ad import read_h5ad, write_h5ad
from microsuite.io.tsv import read_count_matrix, read_matrix_tsv
from microsuite.methods.normalize import normalize_native
from microsuite.methods.table_io import export_table, normalize_table


def _counts(tmp_path: Path) -> Path:
    p = tmp_path / "counts.tsv"
    p.write_text("feature_id\ts1\ts2\nA\t5\t1\nB\t3\t9\n", encoding="utf-8")
    return p


def test_export_table_roundtrip(tmp_path: Path) -> None:
    src = tmp_path / "in.h5ad"
    write_h5ad(read_matrix_tsv(_counts(tmp_path)), src)
    out = tmp_path / "x.tsv"
    export_table(table=src, output=out)
    m = read_count_matrix(out)
    assert list(m.index) == ["A", "B"]
    assert list(m.columns) == ["s1", "s2"]
    assert m.loc["B", "s2"] == 9


def test_export_table_layer_and_metadata(tmp_path: Path) -> None:
    adata = read_matrix_tsv(_counts(tmp_path))
    adata.layers["clr"] = normalize_native(adata, method="clr").X
    adata.obs["group"] = ["a", "b"]
    src = tmp_path / "in.h5ad"
    write_h5ad(adata, src)
    out = tmp_path / "clr.tsv"
    meta = tmp_path / "meta.tsv"
    export_table(table=src, output=out, layer="clr", metadata=meta)
    exported = read_count_matrix(out)
    assert np.allclose(exported.to_numpy(), adata.layers["clr"].T)
    assert meta.exists() and "group" in meta.read_text()


def test_export_table_missing_layer_errors(tmp_path: Path) -> None:
    src = tmp_path / "in.h5ad"
    write_h5ad(read_matrix_tsv(_counts(tmp_path)), src)
    with pytest.raises(MicrobiomeSuiteError, match="layer"):
        export_table(table=src, output=tmp_path / "x.tsv", layer="nope")


def test_normalize_table_clr_tsv_to_tsv(tmp_path: Path) -> None:
    out = tmp_path / "clr.tsv"
    normalize_table(method="clr", input_path=_counts(tmp_path), output=out)
    exported = read_count_matrix(out)
    expected = normalize_native(read_matrix_tsv(_counts(tmp_path)), method="clr").X
    assert np.allclose(exported.to_numpy(), expected.T)


def test_table_cli_smoke(tmp_path: Path) -> None:
    runner = CliRunner()
    out = tmp_path / "clr.tsv"
    result = runner.invoke(
        app, ["table", "normalize", "--method", "clr", "--input", str(_counts(tmp_path)), "-o", str(out)]
    )
    assert result.exit_code == 0, result.stdout
    assert out.exists()
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_table_io.py -v`
Expected: FAIL (`ModuleNotFoundError: microsuite.methods.table_io`).

- [ ] **Step 3: Create `methods/table_io.py`**

```python
from __future__ import annotations

from pathlib import Path

import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.diversity._matrix import dense_counts
from microsuite.io.h5ad import read_h5ad
from microsuite.io.tsv import read_matrix_tsv
from microsuite.methods.normalize import normalize_native


def _matrix_frame(matrix, var_names, obs_names) -> pd.DataFrame:
    frame = pd.DataFrame(
        matrix.T,
        index=pd.Index([str(v) for v in var_names], name="feature_id"),
        columns=[str(s) for s in obs_names],
    )
    return frame


def export_table(
    *,
    table: Path,
    output: Path,
    layer: str | None = None,
    metadata: Path | None = None,
    force: bool = False,
) -> None:
    adata = read_h5ad(ensure_input(table))
    if layer is not None:
        if layer not in adata.layers:
            available = ", ".join(adata.layers.keys()) or "(none)"
            raise MicrobiomeSuiteError(
                f"Layer '{layer}' not found; available layers: {available}."
            )
        matrix = adata.layers[layer]
    else:
        matrix = dense_counts(adata)
    frame = _matrix_frame(matrix, adata.var_names, adata.obs_names)
    frame.to_csv(prepare_output(output, force=force), sep="\t")
    if metadata is not None:
        obs = adata.obs.copy()
        obs.index = obs.index.astype(str)
        obs.index.name = "sample"
        obs.to_csv(prepare_output(metadata, force=force), sep="\t")


def normalize_table(
    *,
    method: str,
    input_path: Path,
    output: Path,
    target_sum: float = 1_000_000.0,
    pseudocount: float = 1.0,
    min_prevalence: float = 0.1,
    force: bool = False,
) -> None:
    adata = read_matrix_tsv(ensure_input(input_path))
    result = normalize_native(
        adata,
        method=method,
        target_sum=target_sum,
        pseudocount=pseudocount,
        min_prevalence=min_prevalence,
    )
    frame = _matrix_frame(dense_counts(result), result.var_names, result.obs_names)
    frame.to_csv(prepare_output(output, force=force), sep="\t")
```

- [ ] **Step 4: Create `cli/table_cmd.py`**

```python
from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.methods.table_io import export_table, normalize_table

app = typer.Typer(
    help="Transform and export feature/profile tables as TSV.", no_args_is_help=True
)


@app.command("export")
def export_cmd(
    table: Annotated[Path, typer.Option("--table", help="Input .h5ad table.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output matrix TSV.")],
    layer: Annotated[str | None, typer.Option("--layer", help="Export this layer instead of X.")] = None,
    metadata: Annotated[
        Path | None, typer.Option("--metadata", help="Also write obs metadata to this TSV.")
    ] = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    export_table(table=table, output=output, layer=layer, metadata=metadata, force=force)


@app.command("normalize")
def normalize_cmd(
    method: Annotated[
        str, typer.Option("--method", help="relative, total-sum, clr, or prevalence-filter.")
    ],
    input_path: Annotated[Path, typer.Option("--input", help="Input matrix TSV.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output matrix TSV.")],
    target_sum: Annotated[float, typer.Option("--target-sum")] = 1_000_000.0,
    pseudocount: Annotated[float, typer.Option("--pseudocount")] = 1.0,
    min_prevalence: Annotated[float, typer.Option("--min-prevalence")] = 0.1,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    normalize_table(
        method=method,
        input_path=input_path,
        output=output,
        target_sum=target_sum,
        pseudocount=pseudocount,
        min_prevalence=min_prevalence,
        force=force,
    )
```

- [ ] **Step 5: Register the group in `cli/app.py`**

Add `table_cmd` to the `from microsuite.cli import (...)` block and, in `_install_groups()`, add `app.add_typer(table_cmd.app, name="table")`.

- [ ] **Step 6: Run to verify pass**

Run: `uv run pytest tests/test_table_io.py -v`
Expected: PASS (5).

- [ ] **Step 7: Full suite + lint gates**

Run: `uv run pytest -q` (all green), then `uv run ruff check .` and `uv run ruff format --check .` (both clean; run `uv run ruff format .` if needed and re-check).

- [ ] **Step 8: Commit**

```bash
git add src/microsuite/methods/table_io.py src/microsuite/cli/table_cmd.py src/microsuite/cli/app.py tests/test_table_io.py
git commit -m "feat(table): tsv-native table export + normalize command group"
```

---

## Self-Review

**Spec coverage:**
- #15 sanitize in shared reader, benefits `import tsv` → Task 1 (Steps 3, 5). ✓
- Conversion 1 (tsv+metadata→h5ad) `import tsv` still works, now #15-safe → Task 1 Step 5. ✓
- Conversion 2 (h5ad→h5ad, layer, X preserved; prevalence-filter structural) → Task 2. ✓
- Conversion 3 (tsv→tsv) `table normalize` → Task 3. ✓
- Conversion 4 (h5ad→tsv + optional metadata + layer) `table export` → Task 3. ✓
- Both CI gates → Task 3 Step 7. ✓

**Placeholder scan:** none — full code for readers, normalize change, both `table_io` methods, the CLI module, and all tests.

**Consistency:** matrices are features×samples on disk everywhere; `_matrix_frame` (var_names index `feature_id`, obs_names columns) is the single write shape used by both `export_table` and `normalize_table`, and `read_count_matrix` is the single read shape; `SHAPE_PRESERVING`/`layer[method]` naming aligns between Task 2's `normalize()` and Task 3's `export --layer`. `dense_counts` is samples×features so every disk write transposes (`.T`), matching every disk read's `.T` into AnnData.
