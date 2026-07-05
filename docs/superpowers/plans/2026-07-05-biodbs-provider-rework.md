# biodbs Provider Rework Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `BiodbsProvider`'s placeholder interface (`biodbs.amplicon.fetch_reference`, which does not exist) with real per-DB dispatch against biodbs **v0.4.0**, covering HOMD, SILVA, GTDB, GreenGenes, UNITE, and PR2, so `microsuite refdb fetch <db> --provider biodbs` actually acquires reference sequences + taxonomy.

**Architecture:** biodbs v0.4.0 exposes heterogeneous per-DB fetch functions (each returns a downloaded-file `Path` and/or a `*TableData` table object), not a uniform `(seqs, tax)` call. `BiodbsProvider.fetch` therefore dispatches over a `_DB_ADAPTERS` table of small per-DB adapter callables, each normalizing that DB's outputs into a `RawRefDb(sequences, taxonomy, qza=None)`. A lazy `_load_biodbs()` guard keeps biodbs an optional dependency and makes the real import monkeypatchable in unit tests. Real network fetches are exercised only in opt-in integration tests; unit tests use a fake biodbs namespace.

**Tech Stack:** Python 3.12, biodbs v0.4.0 (optional extra), polars (biodbs's dep), pytest. The qza short-circuit added to `RefDbProvider.build()` (prior fix) is reused for any DB path that yields a pre-built `.qza`.

## Global Constraints

- `from __future__ import annotations` at the top of every new/edited module.
- All failure paths raise `MicrobiomeSuiteError` from `microsuite._errors` with an actionable message.
- biodbs stays an OPTIONAL dependency: core `import microsuite.refdb` must not import biodbs eagerly; the real import happens only inside `_load_biodbs()` at fetch time. A missing biodbs raises `MicrobiomeSuiteError` telling the user to install the `refdb` extra.
- The `biodbs` provider remains the DEFAULT and self-registers on import (unchanged).
- Unit tests run fully offline by monkeypatching the module-level `_load_biodbs`. Any test that hits the real network is marked with the existing external-integration marker and is opt-in.
- `RawRefDb.taxonomy` is a TSV whose first column is the sequence/record id (matching the `merge_raw`/build expectations of the existing subsystem).

---

## Sequencing note

This plan MUST start only after the qza-short-circuit fix (commit on branch
`refdb-subsystem-design`) is committed and reviewed. Task R1 installs biodbs and
PROBES the real return shapes; later tasks' exact per-DB calls are finalized from
the probe artifact it writes, so R1 is a hard prerequisite for R2–R4.

---

### Task R1: Add biodbs optional dependency and capture a probe artifact

**Files:**
- Modify: `pyproject.toml` (add optional extra `refdb = ["biodbs>=0.4.0"]`)
- Create: `docs/superpowers/specs/biodbs-v040-probe.md` (captured real return shapes)
- Test: `tests/test_refdb_biodbs_available.py`

**Why a probe:** biodbs's per-DB return types are heterogeneous and some base
serializers (`TableData.to_csv`/`as_dataframe`) raise `NotImplementedError` unless
the subclass overrides them. R2–R4 need the EXACT shape (does `homd_get_hmt_lineage()`
support `.as_dataframe()`? what columns? what does `silva_list_current_files` return?)
to be written without guesswork.

- [ ] **Step 1: Add the optional extra**

In `pyproject.toml`, under `[project.optional-dependencies]`, add:
```toml
refdb = ["biodbs>=0.4.0"]
```
If biodbs is not yet on PyPI at 0.4.0, install from the local tag instead and
record the exact command used in the probe doc:
`uv pip install "git+https://github.com/qwerty239qwe/biodbs.git@biodbs_v0.4.0"`

- [ ] **Step 2: Install and verify import**

Run: `uv pip install "git+https://github.com/qwerty239qwe/biodbs.git@biodbs_v0.4.0"` (or `uv sync --extra refdb` once published)
Then: `uv run python -c "import biodbs; print(biodbs.__version__)"`
Expected: prints `0.4.0`

- [ ] **Step 3: Write the availability test**

```python
# tests/test_refdb_biodbs_available.py
from __future__ import annotations

import pytest


def test_biodbs_v040_importable() -> None:
    biodbs = pytest.importorskip("biodbs")
    assert biodbs.__version__.startswith("0.4")
    # the amplicon fetch functions the rework depends on must exist
    for fn in (
        "homd_download_16s_refseq",
        "homd_get_hmt_lineage",
        "silva_download_file",
        "silva_list_current_files",
        "gtdb_download_taxonomy",
        "greengenes_download_file",
        "unite_download",
        "pr2_download_asset",
    ):
        assert hasattr(biodbs, fn), fn
```

- [ ] **Step 4: Run the probe and capture shapes**

Run a throwaway probe script (do NOT commit the script) that, for each of HOMD /
SILVA / GTDB / GreenGenes / UNITE / PR2, calls the cheapest listing/metadata
function and one taxonomy-table getter, and records: the return type, whether
`.as_dataframe()` works and its columns, whether `.to_csv(path)` works, and the
smallest real file that can be downloaded for an integration test. Write the
findings into `docs/superpowers/specs/biodbs-v040-probe.md` as a table per DB.
This doc is the source of truth for R2–R4.

- [ ] **Step 5: Run the test**

Run: `uv run pytest tests/test_refdb_biodbs_available.py -v`
Expected: PASS (or SKIP if biodbs unavailable in CI — `importorskip` handles that).

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml docs/superpowers/specs/biodbs-v040-probe.md tests/test_refdb_biodbs_available.py
git commit -m "feat(refdb): add biodbs v0.4.0 optional extra + capture API probe"
```

---

### Task R2: Rework BiodbsProvider to per-DB dispatch + HOMD adapter

**Files:**
- Rewrite: `src/microsuite/refdb/providers/biodbs.py`
- Test: `tests/test_refdb_provider_biodbs.py` (replace the stub-era tests)

**Interfaces:**
- Produces (module-level, monkeypatchable):
  - `_load_biodbs()` — returns the imported `biodbs` module; raises `MicrobiomeSuiteError` on `ImportError`.
  - `_DB_ADAPTERS: dict[str, Callable[[object, RefDbSpec, Path], RawRefDb]]` keyed by lowercased db name (first key: `"homd"`).
  - `BiodbsProvider.fetch(spec, out_dir)` — looks up `spec.name.lower()` in `_DB_ADAPTERS`; unknown name raises `MicrobiomeSuiteError` listing supported DBs; otherwise calls the adapter with the loaded biodbs module.

**HOMD adapter behavior — CORRECTED per the R1 probe (`docs/superpowers/specs/biodbs-v040-probe.md`):**
The convenience wrappers `homd_download_16s_refseq()` / `homd_get_hmt_lineage()` are BROKEN against the live HOMD site (see probe TL;DR #2 and #5). The working path is the generic `biodbs.homd_download_file(path_or_url, dest, overwrite=False) -> Path` against the real `current/` subdirectory, downloading BOTH:
- sequences: `ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/current/HOMD_16S_rRNA_RefSeq_V16.03.fasta`
- taxonomy: the sibling `.../current/HOMD_16S_rRNA_RefSeq_V16.03.qiime.taxonomy` — this file is ALREADY a QIIME-style id-first `seqID<TAB>lineage` TSV, so it is used directly as `RawRefDb.taxonomy` with NO `*TableData`/`as_dataframe` serialization.
`homd_download_file`'s `dest` may be a directory (it appends the URL basename) — pass `out_dir` (already created by `fetch()`). Store the current filename stem as a module constant `_HOMD_STEM = "HOMD_16S_rRNA_RefSeq_V16.03"`; `spec.version` is ignored (HOMD serves a single `current` set; the stem tracks HOMD's `current/` symlink and R5's live test will catch drift). The `_write_taxonomy_tsv` helper is NOT needed here — it is introduced in R3 for GTDB, the first DB whose taxonomy arrives as a `*TableData`.

- [ ] **Step 1: Write the failing unit tests (offline, fake biodbs)**

```python
# tests/test_refdb_provider_biodbs.py
from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import get_provider
from microsuite.refdb.providers import biodbs as _biodbs  # noqa: F401  (registration)
from microsuite.refdb.spec import RefDbSpec


class _FakeBiodbs:
    """Mimics biodbs.homd_download_file: dest is a dir, basename taken from the URL path.
    Writes fixture content for the two HOMD files the adapter requests."""

    _CONTENT = {
        "HOMD_16S_rRNA_RefSeq_V16.03.fasta": ">seqA\nACGT\n>seqB\nTGCA\n",
        "HOMD_16S_rRNA_RefSeq_V16.03.qiime.taxonomy": "seqA\tk__B;s__x\nseqB\tk__B;s__y\n",
    }

    def homd_download_file(self, path_or_url, dest, overwrite=False) -> Path:
        name = Path(path_or_url).name
        target = Path(dest) / name
        target.write_text(self._CONTENT[name], encoding="utf-8")
        return target


def test_biodbs_is_default_provider() -> None:
    assert RefDbSpec(name="homd", version="15.22").provider == "biodbs"


def test_homd_adapter_produces_seqs_and_taxonomy(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(_biodbs, "_load_biodbs", lambda: _FakeBiodbs())
    provider = get_provider("biodbs")
    raw = provider.fetch(RefDbSpec(name="homd", version="15.22"), out_dir=tmp_path)
    assert raw.sequences.exists()
    ids = [l[1:].strip() for l in raw.sequences.read_text().splitlines() if l.startswith(">")]
    assert ids == ["seqA", "seqB"]
    tax_first_col = [r.split("\t")[0] for r in raw.taxonomy.read_text().splitlines() if r.strip()]
    assert tax_first_col == ["seqA", "seqB"]


def test_unknown_db_raises(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(_biodbs, "_load_biodbs", lambda: _FakeBiodbs())
    provider = get_provider("biodbs")
    with pytest.raises(MicrobiomeSuiteError):
        provider.fetch(RefDbSpec(name="not-a-db", version="1"), out_dir=tmp_path)


def test_missing_biodbs_raises(tmp_path: Path, monkeypatch) -> None:
    def boom():
        raise MicrobiomeSuiteError("biodbs not installed")
    monkeypatch.setattr(_biodbs, "_load_biodbs", boom)
    provider = get_provider("biodbs")
    with pytest.raises(MicrobiomeSuiteError):
        provider.fetch(RefDbSpec(name="homd", version="15.22"), out_dir=tmp_path)
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_refdb_provider_biodbs.py -v`
Expected: FAIL (old stub `_load_biodbs_fetch` gone / new dispatch not present).

- [ ] **Step 3: Rewrite the provider**

```python
# src/microsuite/refdb/providers/biodbs.py
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import register_provider
from microsuite.refdb.providers._base import RefDbProvider
from microsuite.refdb.spec import RawRefDb, RefDbSpec


def _load_biodbs():
    try:
        import biodbs  # type: ignore[import-not-found]
    except ImportError as exc:
        raise MicrobiomeSuiteError(
            "The default 'biodbs' provider requires biodbs>=0.4.0. Install the "
            "refdb extra (e.g. `uv sync --extra refdb`), pass a raw --classifier "
            "path, or use --provider rescript."
        ) from exc
    return biodbs


# HOMD serves 16S RefSeq under a `current/` symlink dir; the wrapper
# homd_download_16s_refseq() is broken (probe TL;DR #2/#5), so download the
# fasta and its sibling QIIME-style taxonomy via the generic homd_download_file.
_HOMD_BASE = "ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/current"
_HOMD_STEM = "HOMD_16S_rRNA_RefSeq_V16.03"


def _homd_adapter(bd, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
    seqs = Path(bd.homd_download_file(f"{_HOMD_BASE}/{_HOMD_STEM}.fasta", str(out_dir)))
    # The .qiime.taxonomy file is already an id-first `seqID<TAB>lineage` TSV.
    tax = Path(bd.homd_download_file(f"{_HOMD_BASE}/{_HOMD_STEM}.qiime.taxonomy", str(out_dir)))
    return RawRefDb(sequences=seqs, taxonomy=tax)


_DB_ADAPTERS: dict[str, Callable[[object, RefDbSpec, Path], RawRefDb]] = {
    "homd": _homd_adapter,
}


class BiodbsProvider(RefDbProvider):
    name = "biodbs"

    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
        out_dir.mkdir(parents=True, exist_ok=True)
        adapter = _DB_ADAPTERS.get(spec.name.lower())
        if adapter is None:
            supported = ", ".join(sorted(_DB_ADAPTERS))
            raise MicrobiomeSuiteError(
                f"biodbs provider has no adapter for DB '{spec.name}'. Supported: {supported}."
            )
        bd = _load_biodbs()
        return adapter(bd, spec, out_dir)


register_provider(BiodbsProvider())
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_refdb_provider_biodbs.py -v`
Expected: PASS (4 passed)

- [ ] **Step 5: Full refdb suite regression check**

Run: `uv run pytest tests/test_refdb_*.py -q`
Expected: PASS (no regression). Report the count.

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/refdb/providers/biodbs.py tests/test_refdb_provider_biodbs.py
git commit -m "feat(refdb): biodbs per-DB dispatch with HOMD adapter"
```

---

### Task R3: SILVA + GTDB adapters

**Files:**
- Modify: `src/microsuite/refdb/providers/biodbs.py` (add `_silva_adapter`, `_gtdb_adapter`, register in `_DB_ADAPTERS`)
- Test: `tests/test_refdb_provider_biodbs.py` (add SILVA + GTDB cases with fake biodbs)

**Interfaces (finalize exact file paths from the R1 probe doc):**
- SILVA: use `bd.silva_list_current_files(...)` to locate the SSU Ref NR99 FASTA + taxonomy, then `bd.silva_download_file(path, out_dir)` for each; `spec.version` selects release. Return `RawRefDb(sequences, taxonomy)`.
- GTDB: `bd.gtdb_download_file(<seq path>, out_dir)` for sequences and `bd.gtdb_download_taxonomy(domain="bac120", dest=out_dir, release=spec.version)` for taxonomy; return `RawRefDb(sequences, taxonomy)`.

- [ ] **Step 1: Write failing tests (fake biodbs exposing silva_*/gtdb_* used by the adapters)** — assert each adapter writes a sequences FASTA and an id-first taxonomy TSV; unknown-DB path already covered in R2.
- [ ] **Step 2: Run to verify fail.** `uv run pytest tests/test_refdb_provider_biodbs.py -v`
- [ ] **Step 3: Implement `_silva_adapter` and `_gtdb_adapter`; add both to `_DB_ADAPTERS`.** Use the exact biodbs call arguments recorded in `docs/superpowers/specs/biodbs-v040-probe.md`.
- [ ] **Step 4: Run to verify pass.**
- [ ] **Step 5: Full refdb suite regression.** `uv run pytest tests/test_refdb_*.py -q`
- [ ] **Step 6: Commit.** `git commit -m "feat(refdb): SILVA and GTDB biodbs adapters"`

---

### Task R4: GreenGenes + UNITE + PR2 adapters

**Files:**
- Modify: `src/microsuite/refdb/providers/biodbs.py` (add three adapters, register them)
- Test: `tests/test_refdb_provider_biodbs.py` (three fake-biodbs cases)

**Interfaces (finalize exact args from the R1 probe doc):**
- GreenGenes: `bd.greengenes_list_files(...)` → locate seqs + taxonomy, `bd.greengenes_download_file(path, out_dir)`.
- UNITE: `bd.unite_download(version=spec.version, dest=out_dir, taxon_group="fungi")` → a bundled archive; adapter extracts/points sequences + taxonomy.
- PR2: `bd.pr2_download_asset(...)` → seqs + taxonomy assets.

- [ ] **Step 1: Write failing tests** (fake biodbs per DB; assert seqs FASTA + id-first taxonomy TSV).
- [ ] **Step 2: Run to verify fail.**
- [ ] **Step 3: Implement the three adapters; register in `_DB_ADAPTERS`.**
- [ ] **Step 4: Run to verify pass.**
- [ ] **Step 5: Full refdb suite regression.**
- [ ] **Step 6: Commit.** `git commit -m "feat(refdb): GreenGenes, UNITE, PR2 biodbs adapters"`

---

### Task R5: Opt-in live-integration tests

**Files:**
- Create: `tests/integration/test_refdb_biodbs_live.py`

**Interfaces:**
- Consumes: `microsuite.refdb.service.fetch_refdb`, the real biodbs package, network.
- Behavior: one test per DB, marked with the repo's existing external-integration marker (inspect `tests/integration/` and `pyproject.toml`/`conftest.py` for the exact marker name and skip convention). Each test does the smallest real acquisition identified in the R1 probe (e.g. HOMD 16S RefSeq → build vsearch → assert the artifact FASTA is non-empty and the registry records it). Guarded so a normal `uv run pytest` (no marker/opt-in flag) skips them.

- [ ] **Step 1: Confirm the marker convention.** Read `tests/integration/` + `conftest.py` for how existing external-tool tests opt in; reuse that exact mechanism (do not invent a new one).
- [ ] **Step 2: Write one live test (HOMD end-to-end).** `fetch_refdb(RefDbSpec(name="homd", version="15.22"), "vsearch")` → assert artifact exists and FASTA has ≥1 record.
- [ ] **Step 3: Verify it SKIPS by default.** `uv run pytest tests/integration/test_refdb_biodbs_live.py -q` → skipped without the opt-in flag.
- [ ] **Step 4: Add the remaining DB live tests** (SILVA/GTDB/GreenGenes/UNITE/PR2), each smallest-possible.
- [ ] **Step 5: Commit.** `git commit -m "test(refdb): opt-in live biodbs integration tests"`

---

## Self-Review

**Spec coverage (against the user decision "fix + full biodbs rework"):**
- Per-DB dispatch replacing the stub → R2 (architecture + HOMD), R3, R4. ✓
- All six target DBs (HOMD/SILVA/GTDB/GreenGenes/UNITE/PR2) → R2/R3/R4. ✓
- biodbs as real optional dependency → R1. ✓
- qza short-circuit reuse for any `.qza`-yielding path → provided by the prior fix (base `build()`); SILVA classifier variant deferred unless the probe shows it's needed. ✓
- Live-integration tests behind opt-in marker → R5. ✓
- Offline unit tests via monkeypatched `_load_biodbs` → R2/R3/R4. ✓

**Known gaps recorded for the human (not in this rework's scope):** biodbs v0.4.0 has NO MOMD fetcher and NO NCBI 16S Microbial *sequence* download, so the full FOMC combined DB (HOMD+MOMD+NCBI) still cannot be built from biodbs alone — the HOMD arm works; MOMD + NCBI 16S need a separate source. This is a subject for the FOMC-combined follow-up, not this rework.

**Placeholder note:** R3/R4/R5 intentionally defer EXACT per-DB call arguments to the R1 probe artifact rather than guessing file paths that would be wrong. R1 is a hard prerequisite; its probe doc makes R3/R4 concrete before they are implemented.
