# refdb Subsystem Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a reference-database subsystem (`microsuite.refdb`) that fetches, builds, caches, and registers 16S reference DBs behind a provider interface, with biodbs as the default provider.

**Architecture:** A new `src/microsuite/refdb/` package defines immutable dataclasses (`RefDbSpec`, `RawRefDb`, `BuiltArtifact`), a JSON-manifest registry for caching built artifacts, a `build.py` layer that turns FASTA+taxonomy into vsearch/BLAST/QIIME2 artifacts, and a provider interface with two implementations (`biodbs` default, `rescript` optional). A service orchestrator wires fetch→merge→build→register and resolves `refdb:` classifier references. A `microsuite refdb` CLI group and a `tax_classify --classifier` integration expose it.

**Tech Stack:** Python 3.12, Typer (CLI), pytest + `typer.testing.CliRunner`, dataclasses, `subprocess` via the repo's `run_command`. External tools (`makeblastdb`, `qiime`, `vsearch`, biodbs, q2-rescript) are invoked but never run in unit tests — they are argv-asserted with monkeypatched `shutil.which`/`subprocess.run`, matching `tests/test_denoise_cluster_methods.py`.

## Global Constraints

- All raised errors use `MicrobiomeSuiteError` from `microsuite._errors` with an actionable message (missing tool, missing DB, bad version).
- External-tool invocation goes through `run_command` from `microsuite.runtime.runner` with a `CommandLog`; never call `subprocess.run` directly in library code.
- Use `ensure_input` / `prepare_output` from `microsuite._paths` for path checks.
- Backwards compatibility: a raw `--classifier` filesystem path must keep working unchanged.
- Unit tests must run fully offline; real external-tool paths are opt-in and gated by the existing integration marker.
- `from __future__ import annotations` at the top of every new module (repo convention).
- Every new public function callable from Python is exported through `microsuite.api` where the repo already does so for methods.

---

### Task 1: Core types and cache directory

**Files:**
- Create: `src/microsuite/refdb/__init__.py`
- Create: `src/microsuite/refdb/spec.py`
- Create: `src/microsuite/refdb/paths.py`
- Test: `tests/test_refdb_types.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `RefDbSource(name: str, version: str, locator: str | None = None)` — frozen dataclass.
  - `RefDbSpec(name: str, version: str, provider: str = "biodbs", target: str = "16S", build_targets: tuple[str, ...] = ("vsearch",), sources: tuple[RefDbSource, ...] = ())` — frozen dataclass.
  - `RawRefDb(sequences: Path, taxonomy: Path, qza: Path | None = None)` — frozen dataclass.
  - `BuiltArtifact(path: Path, build_target: str, checksum: str)` — frozen dataclass.
  - `refdb_cache_dir() -> Path` — returns `$MICROSUITE_REFDB_DIR` if set, else `~/.cache/microsuite/refdb`; creates it.
  - `VALID_BUILD_TARGETS: tuple[str, ...] = ("vsearch", "blast", "qiime2")`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_refdb_types.py
from __future__ import annotations

from pathlib import Path

from microsuite.refdb.paths import VALID_BUILD_TARGETS, refdb_cache_dir
from microsuite.refdb.spec import BuiltArtifact, RawRefDb, RefDbSource, RefDbSpec


def test_refdbspec_defaults() -> None:
    spec = RefDbSpec(name="homd", version="15.22")
    assert spec.provider == "biodbs"
    assert spec.target == "16S"
    assert spec.build_targets == ("vsearch",)
    assert spec.sources == ()


def test_refdbspec_is_frozen() -> None:
    spec = RefDbSpec(name="homd", version="15.22")
    try:
        spec.name = "other"  # type: ignore[misc]
    except AttributeError:
        return
    raise AssertionError("RefDbSpec should be frozen")


def test_valid_build_targets() -> None:
    assert VALID_BUILD_TARGETS == ("vsearch", "blast", "qiime2")


def test_cache_dir_honors_env(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_REFDB_DIR", str(tmp_path / "cache"))
    root = refdb_cache_dir()
    assert root == tmp_path / "cache"
    assert root.is_dir()


def test_cache_dir_default(monkeypatch) -> None:
    monkeypatch.delenv("MICROSUITE_REFDB_DIR", raising=False)
    root = refdb_cache_dir()
    assert root == Path.home() / ".cache" / "microsuite" / "refdb"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_refdb_types.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'microsuite.refdb'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/microsuite/refdb/__init__.py
from __future__ import annotations
```

```python
# src/microsuite/refdb/spec.py
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class RefDbSource:
    name: str
    version: str
    locator: str | None = None


@dataclass(frozen=True)
class RefDbSpec:
    name: str
    version: str
    provider: str = "biodbs"
    target: str = "16S"
    build_targets: tuple[str, ...] = ("vsearch",)
    sources: tuple[RefDbSource, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class RawRefDb:
    sequences: Path
    taxonomy: Path
    qza: Path | None = None


@dataclass(frozen=True)
class BuiltArtifact:
    path: Path
    build_target: str
    checksum: str
```

```python
# src/microsuite/refdb/paths.py
from __future__ import annotations

import os
from pathlib import Path

VALID_BUILD_TARGETS: tuple[str, ...] = ("vsearch", "blast", "qiime2")


def refdb_cache_dir() -> Path:
    env = os.environ.get("MICROSUITE_REFDB_DIR")
    root = Path(env) if env else Path.home() / ".cache" / "microsuite" / "refdb"
    root.mkdir(parents=True, exist_ok=True)
    return root
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_refdb_types.py -v`
Expected: PASS (5 passed)

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/refdb/__init__.py src/microsuite/refdb/spec.py src/microsuite/refdb/paths.py tests/test_refdb_types.py
git commit -m "feat(refdb): core dataclasses and cache-dir resolution"
```

---

### Task 2: Registry (JSON manifest cache)

**Files:**
- Create: `src/microsuite/refdb/registry.py`
- Test: `tests/test_refdb_registry.py`

**Interfaces:**
- Consumes: `BuiltArtifact` (Task 1).
- Produces:
  - `sha256_file(path: Path) -> str` — hex digest of a file's bytes.
  - `class RefDbRegistry` with:
    - `__init__(self, root: Path)` — manifest lives at `root / "manifest.json"`.
    - `resolve(self, name: str, version: str, build_target: str) -> BuiltArtifact | None` — returns the recorded artifact if present AND the file exists AND its checksum matches; otherwise `None`.
    - `record(self, name: str, version: str, artifact: BuiltArtifact, provider: str) -> None` — writes/updates the manifest entry.
  - Manifest key format: `f"{name}@{version}:{build_target}"`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_refdb_registry.py
from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.registry import RefDbRegistry, sha256_file
from microsuite.refdb.spec import BuiltArtifact


def _artifact(tmp_path: Path, text: str = "ACGT") -> BuiltArtifact:
    path = tmp_path / "db.fasta"
    path.write_text(text, encoding="utf-8")
    return BuiltArtifact(path=path, build_target="vsearch", checksum=sha256_file(path))


def test_record_then_resolve_round_trip(tmp_path: Path) -> None:
    reg = RefDbRegistry(tmp_path / "cache")
    art = _artifact(tmp_path)
    reg.record("homd", "15.22", art, provider="biodbs")

    resolved = reg.resolve("homd", "15.22", "vsearch")
    assert resolved is not None
    assert resolved.path == art.path
    assert resolved.checksum == art.checksum


def test_resolve_missing_returns_none(tmp_path: Path) -> None:
    reg = RefDbRegistry(tmp_path / "cache")
    assert reg.resolve("nope", "1", "vsearch") is None


def test_resolve_checksum_mismatch_returns_none(tmp_path: Path) -> None:
    reg = RefDbRegistry(tmp_path / "cache")
    art = _artifact(tmp_path)
    reg.record("homd", "15.22", art, provider="biodbs")
    art.path.write_text("MUTATED", encoding="utf-8")  # invalidate
    assert reg.resolve("homd", "15.22", "vsearch") is None


def test_corrupt_manifest_raises(tmp_path: Path) -> None:
    root = tmp_path / "cache"
    root.mkdir(parents=True)
    (root / "manifest.json").write_text("{ not json", encoding="utf-8")
    reg = RefDbRegistry(root)
    with pytest.raises(MicrobiomeSuiteError):
        reg.resolve("homd", "15.22", "vsearch")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_refdb_registry.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'microsuite.refdb.registry'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/microsuite/refdb/registry.py
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.spec import BuiltArtifact


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


class RefDbRegistry:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)
        self.manifest_path = self.root / "manifest.json"

    def _key(self, name: str, version: str, build_target: str) -> str:
        return f"{name}@{version}:{build_target}"

    def _load(self) -> dict[str, dict[str, str]]:
        if not self.manifest_path.exists():
            return {}
        try:
            return json.loads(self.manifest_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError) as exc:
            raise MicrobiomeSuiteError(
                f"Reference-DB manifest is corrupt: {self.manifest_path}. "
                "Delete it to rebuild the cache."
            ) from exc

    def resolve(self, name: str, version: str, build_target: str) -> BuiltArtifact | None:
        entry = self._load().get(self._key(name, version, build_target))
        if entry is None:
            return None
        path = Path(entry["path"])
        if not path.exists() or sha256_file(path) != entry["checksum"]:
            return None
        return BuiltArtifact(path=path, build_target=build_target, checksum=entry["checksum"])

    def record(
        self, name: str, version: str, artifact: BuiltArtifact, provider: str
    ) -> None:
        manifest = self._load()
        manifest[self._key(name, version, artifact.build_target)] = {
            "path": str(artifact.path),
            "checksum": artifact.checksum,
            "provider": provider,
        }
        self.manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
        )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_refdb_registry.py -v`
Expected: PASS (4 passed)

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/refdb/registry.py tests/test_refdb_registry.py
git commit -m "feat(refdb): JSON-manifest registry with checksum validation"
```

---

### Task 3: Build layer (merge + vsearch real; blast/qiime2 argv-asserted) and fixtures

**Files:**
- Create: `src/microsuite/refdb/build.py`
- Create: `src/microsuite/data/fixtures/refdb_mock/source_a.fasta`
- Create: `src/microsuite/data/fixtures/refdb_mock/source_a.tax.tsv`
- Create: `src/microsuite/data/fixtures/refdb_mock/source_b.fasta`
- Create: `src/microsuite/data/fixtures/refdb_mock/source_b.tax.tsv`
- Create: `src/microsuite/data/fixtures/refdb_mock/README.md`
- Test: `tests/test_refdb_build.py`

**Interfaces:**
- Consumes: `RawRefDb`, `BuiltArtifact` (Task 1), `sha256_file` (Task 2).
- Produces:
  - `merge_raw(raws: list[RawRefDb], out_dir: Path) -> RawRefDb` — concatenates FASTAs and taxonomy TSVs, de-duplicating by sequence id (first occurrence wins); writes `merged.fasta` + `merged.tax.tsv` under `out_dir`.
  - `build_artifact(raw: RawRefDb, build_target: str, out_dir: Path, *, force: bool = False, run_dir: Path | None = None, timeout: float | None = None) -> BuiltArtifact` — `vsearch` writes/validates a reference FASTA and returns it (no external tool); `blast` invokes `makeblastdb`; `qiime2` invokes `qiime tools import`. Unknown target raises `MicrobiomeSuiteError`.

**Fixture format:** FASTA headers are `>seqID`; taxonomy TSV rows are `seqID<TAB>k__...;p__...;...;s__...`. `source_a` and `source_b` share one overlapping id (`seq1`) so the merge/dedup test has something to collapse.

- [ ] **Step 1: Write the fixtures**

```text
# src/microsuite/data/fixtures/refdb_mock/source_a.fasta
>seq1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>seq2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
```

```text
# src/microsuite/data/fixtures/refdb_mock/source_a.tax.tsv
seq1	k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus_mutans
seq2	k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella;s__Prevotella_melaninogenica
```

```text
# src/microsuite/data/fixtures/refdb_mock/source_b.fasta
>seq1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>seq3
GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCC
```

```text
# src/microsuite/data/fixtures/refdb_mock/source_b.tax.tsv
seq1	k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus_mutans
seq3	k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Actinomyces;s__Actinomyces_naeslundii
```

```markdown
# src/microsuite/data/fixtures/refdb_mock/README.md
Tiny mock reference-DB fixture (HOMD-like) for refdb unit tests. Two sources
share `seq1` so merge/dedup collapses to seq1, seq2, seq3. Not real biology.
```

- [ ] **Step 2: Write the failing test**

```python
# tests/test_refdb_build.py
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.build import build_artifact, merge_raw
from microsuite.refdb.spec import RawRefDb

FIXTURE = Path("src/microsuite/data/fixtures/refdb_mock")


def _raw(tag: str) -> RawRefDb:
    return RawRefDb(
        sequences=FIXTURE / f"{tag}.fasta",
        taxonomy=FIXTURE / f"{tag}.tax.tsv",
    )


def test_merge_dedups_by_seq_id(tmp_path: Path) -> None:
    merged = merge_raw([_raw("source_a"), _raw("source_b")], out_dir=tmp_path)
    ids = [
        line[1:].strip()
        for line in merged.sequences.read_text().splitlines()
        if line.startswith(">")
    ]
    assert ids == ["seq1", "seq2", "seq3"]
    tax_ids = [row.split("\t")[0] for row in merged.taxonomy.read_text().splitlines() if row]
    assert tax_ids == ["seq1", "seq2", "seq3"]


def test_build_vsearch_is_offline_and_checksummed(tmp_path: Path) -> None:
    art = build_artifact(_raw("source_a"), "vsearch", out_dir=tmp_path)
    assert art.build_target == "vsearch"
    assert art.path.exists()
    assert len(art.checksum) == 64


def test_build_blast_invokes_makeblastdb(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(
        "shutil.which", lambda name: "makeblastdb" if name == "makeblastdb" else None
    )
    calls: list[list[str]] = []

    def fake_run(command, **kwargs):
        calls.append(command)
        (tmp_path / "blastdb.nhr").write_text("x", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    build_artifact(_raw("source_a"), "blast", out_dir=tmp_path)
    assert calls[0][0] == "makeblastdb"
    assert "-dbtype" in calls[0] and "nucl" in calls[0]


def test_build_unknown_target_raises(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        build_artifact(_raw("source_a"), "bowtie", out_dir=tmp_path)
```

- [ ] **Step 3: Run test to verify it fails**

Run: `uv run pytest tests/test_refdb_build.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'microsuite.refdb.build'`

- [ ] **Step 4: Write minimal implementation**

```python
# src/microsuite/refdb/build.py
from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input
from microsuite.refdb.registry import sha256_file
from microsuite.refdb.spec import BuiltArtifact, RawRefDb
from microsuite.runtime.runner import CommandLog, run_command


def _iter_fasta(path: Path):
    seq_id: str | None = None
    lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if seq_id is not None:
                yield seq_id, lines
            seq_id = line[1:].strip()
            lines = []
        elif seq_id is not None:
            lines.append(line)
    if seq_id is not None:
        yield seq_id, lines


def merge_raw(raws: list[RawRefDb], out_dir: Path) -> RawRefDb:
    out_dir.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    seq_out = out_dir / "merged.fasta"
    tax_out = out_dir / "merged.tax.tsv"
    with seq_out.open("w", encoding="utf-8") as seq_fh, tax_out.open(
        "w", encoding="utf-8"
    ) as tax_fh:
        for raw in raws:
            ensure_input(raw.sequences)
            ensure_input(raw.taxonomy)
            tax_by_id = {
                row.split("\t", 1)[0]: row
                for row in raw.taxonomy.read_text(encoding="utf-8").splitlines()
                if row.strip()
            }
            for seq_id, body in _iter_fasta(raw.sequences):
                if seq_id in seen:
                    continue
                seen.add(seq_id)
                seq_fh.write(f">{seq_id}\n")
                seq_fh.write("\n".join(body) + "\n")
                if seq_id in tax_by_id:
                    tax_fh.write(tax_by_id[seq_id].rstrip("\n") + "\n")
    return RawRefDb(sequences=seq_out, taxonomy=tax_out)


def build_artifact(
    raw: RawRefDb,
    build_target: str,
    out_dir: Path,
    *,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> BuiltArtifact:
    out_dir.mkdir(parents=True, exist_ok=True)
    ensure_input(raw.sequences)
    if build_target == "vsearch":
        target = out_dir / "reference.fasta"
        shutil.copyfile(raw.sequences, target)
        return BuiltArtifact(target, "vsearch", sha256_file(target))
    if build_target == "blast":
        tool = shutil.which("makeblastdb")
        if tool is None:
            raise MicrobiomeSuiteError(
                "BLAST DB build requires 'makeblastdb'. Install BLAST+ or use the "
                "microsuite BLAST container and rerun."
            )
        db_prefix = out_dir / "blastdb"
        run_command(
            [tool, "-in", str(raw.sequences), "-dbtype", "nucl", "-out", str(db_prefix)],
            "makeblastdb failed.",
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(task="refdb_build", backend="blast"),
        )
        marker = db_prefix.with_suffix(".nhr")
        return BuiltArtifact(marker, "blast", sha256_file(marker))
    if build_target == "qiime2":
        tool = shutil.which("qiime")
        if tool is None:
            raise MicrobiomeSuiteError(
                "QIIME2 artifact build requires the 'qiime' command. Activate a "
                "QIIME 2 environment and rerun."
            )
        artifact = out_dir / "reference-seqs.qza"
        run_command(
            [
                tool,
                "tools",
                "import",
                "--type",
                "FeatureData[Sequence]",
                "--input-path",
                str(raw.sequences),
                "--output-path",
                str(artifact),
            ],
            "QIIME 2 reference import failed.",
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(task="refdb_build", backend="qiime2"),
        )
        return BuiltArtifact(artifact, "qiime2", sha256_file(artifact))
    raise MicrobiomeSuiteError(
        f"Unknown build target '{build_target}'. Choose one of: vsearch, blast, qiime2."
    )
```

Note: the `blast`/`qiime2` branches call `run_command`, which calls `subprocess.run`; the test monkeypatches `subprocess.run` to create the expected output marker file so `sha256_file` succeeds.

- [ ] **Step 5: Run test to verify it passes**

Run: `uv run pytest tests/test_refdb_build.py -v`
Expected: PASS (4 passed)

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/refdb/build.py src/microsuite/data/fixtures/refdb_mock tests/test_refdb_build.py
git commit -m "feat(refdb): merge+dedup and vsearch/blast/qiime2 build layer with fixtures"
```

---

### Task 4: Provider interface, registry, and FakeProvider

**Files:**
- Create: `src/microsuite/refdb/providers/__init__.py`
- Create: `src/microsuite/refdb/providers/_base.py`
- Test: `tests/test_refdb_providers.py`

**Interfaces:**
- Consumes: `RefDbSpec`, `RawRefDb`, `BuiltArtifact` (Task 1); `build_artifact` (Task 3).
- Produces:
  - `class RefDbProvider(ABC)` with attribute `name: str`, abstract `fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb`, and concrete `build(self, raw: RawRefDb, build_target: str, out_dir: Path, **kw) -> BuiltArtifact` delegating to `build_artifact`.
  - `get_provider(name: str) -> RefDbProvider` — looks up `_PROVIDERS`; raises `MicrobiomeSuiteError` for unknown names.
  - `register_provider(provider: RefDbProvider) -> None` — inserts into `_PROVIDERS` (used by tests and by Tasks 5–6).
  - `_PROVIDERS: dict[str, RefDbProvider]` — module-level registry, initially empty (populated by Tasks 5–6 on import).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_refdb_providers.py
from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import get_provider, register_provider
from microsuite.refdb.providers._base import RefDbProvider
from microsuite.refdb.spec import RawRefDb, RefDbSpec

FIXTURE = Path("src/microsuite/data/fixtures/refdb_mock")


class FakeProvider(RefDbProvider):
    name = "fake"

    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
        return RawRefDb(
            sequences=FIXTURE / "source_a.fasta",
            taxonomy=FIXTURE / "source_a.tax.tsv",
        )


def test_register_and_get_provider() -> None:
    register_provider(FakeProvider())
    assert isinstance(get_provider("fake"), FakeProvider)


def test_get_unknown_provider_raises() -> None:
    with pytest.raises(MicrobiomeSuiteError):
        get_provider("does-not-exist")


def test_default_build_delegates_to_build_artifact(tmp_path: Path) -> None:
    provider = FakeProvider()
    raw = provider.fetch(RefDbSpec(name="x", version="1"), out_dir=tmp_path)
    art = provider.build(raw, "vsearch", out_dir=tmp_path)
    assert art.build_target == "vsearch"
    assert art.path.exists()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_refdb_providers.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'microsuite.refdb.providers'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/microsuite/refdb/providers/_base.py
from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

from microsuite.refdb.build import build_artifact
from microsuite.refdb.spec import BuiltArtifact, RawRefDb, RefDbSpec


class RefDbProvider(ABC):
    name: str

    @abstractmethod
    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb: ...

    def build(
        self,
        raw: RawRefDb,
        build_target: str,
        out_dir: Path,
        *,
        force: bool = False,
        run_dir: Path | None = None,
        timeout: float | None = None,
    ) -> BuiltArtifact:
        return build_artifact(
            raw, build_target, out_dir, force=force, run_dir=run_dir, timeout=timeout
        )
```

```python
# src/microsuite/refdb/providers/__init__.py
from __future__ import annotations

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers._base import RefDbProvider

_PROVIDERS: dict[str, RefDbProvider] = {}


def register_provider(provider: RefDbProvider) -> None:
    _PROVIDERS[provider.name] = provider


def get_provider(name: str) -> RefDbProvider:
    provider = _PROVIDERS.get(name)
    if provider is None:
        choices = ", ".join(sorted(_PROVIDERS)) or "(none registered)"
        raise MicrobiomeSuiteError(
            f"Unknown reference-DB provider '{name}'. Available: {choices}"
        )
    return provider
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_refdb_providers.py -v`
Expected: PASS (3 passed)

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/refdb/providers/__init__.py src/microsuite/refdb/providers/_base.py tests/test_refdb_providers.py
git commit -m "feat(refdb): provider ABC, registry, and default build delegation"
```

---

### Task 5: RESCRIPt provider (optional, argv-asserted)

**Files:**
- Create: `src/microsuite/refdb/providers/rescript.py`
- Modify: `src/microsuite/refdb/providers/__init__.py` (auto-register on import)
- Test: `tests/test_refdb_provider_rescript.py`

**Interfaces:**
- Consumes: `RefDbProvider`, `register_provider` (Task 4).
- Produces: `class RescriptProvider(RefDbProvider)` with `name = "rescript"`. `fetch` requires `qiime`; for `spec.name == "silva"` it invokes `qiime rescript get-silva-data`, writing `silva-seqs.qza`/`silva-tax.qza`, and returns a `RawRefDb` whose `qza` field points at the seqs artifact. Missing `qiime` raises `MicrobiomeSuiteError`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_refdb_provider_rescript.py
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import get_provider
from microsuite.refdb.providers import rescript as _rescript  # noqa: F401  (force registration)
from microsuite.refdb.spec import RefDbSpec


def test_rescript_silva_builds_command(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)
    calls: list[list[str]] = []

    def fake_run(command, **kwargs):
        calls.append(command)
        for i, tok in enumerate(command):
            if tok == "--o-silva-sequences":
                Path(command[i + 1]).write_text("x", encoding="utf-8")
            if tok == "--o-silva-taxonomy":
                Path(command[i + 1]).write_text("x", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    provider = get_provider("rescript")
    provider.fetch(RefDbSpec(name="silva", version="138.1", provider="rescript"), out_dir=tmp_path)

    assert calls[0][:3] == ["qiime", "rescript", "get-silva-data"]


def test_rescript_requires_qiime(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr("shutil.which", lambda name: None)
    provider = get_provider("rescript")
    with pytest.raises(MicrobiomeSuiteError):
        provider.fetch(RefDbSpec(name="silva", version="138.1", provider="rescript"), out_dir=tmp_path)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_refdb_provider_rescript.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'microsuite.refdb.providers.rescript'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/microsuite/refdb/providers/rescript.py
from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import register_provider
from microsuite.refdb.providers._base import RefDbProvider
from microsuite.refdb.spec import RawRefDb, RefDbSpec
from microsuite.runtime.runner import CommandLog, run_command


class RescriptProvider(RefDbProvider):
    name = "rescript"

    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
        qiime = shutil.which("qiime")
        if qiime is None:
            raise MicrobiomeSuiteError(
                "The 'rescript' provider requires a QIIME 2 environment with the "
                "RESCRIPt plugin (the 'qiime' command was not found)."
            )
        out_dir.mkdir(parents=True, exist_ok=True)
        if spec.name != "silva":
            raise MicrobiomeSuiteError(
                f"The 'rescript' provider does not yet support DB '{spec.name}'. "
                "Supported: silva."
            )
        seqs = out_dir / "silva-seqs.qza"
        tax = out_dir / "silva-tax.qza"
        run_command(
            [
                qiime,
                "rescript",
                "get-silva-data",
                "--p-version",
                spec.version,
                "--o-silva-sequences",
                str(seqs),
                "--o-silva-taxonomy",
                str(tax),
            ],
            "RESCRIPt get-silva-data failed.",
            log=CommandLog(task="refdb_fetch", backend="rescript"),
        )
        return RawRefDb(sequences=seqs, taxonomy=tax, qza=seqs)


register_provider(RescriptProvider())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_refdb_provider_rescript.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/refdb/providers/rescript.py tests/test_refdb_provider_rescript.py
git commit -m "feat(refdb): optional RESCRIPt provider (SILVA)"
```

---

### Task 6: biodbs provider (default, lazy import with guard)

**Files:**
- Create: `src/microsuite/refdb/providers/biodbs.py`
- Test: `tests/test_refdb_provider_biodbs.py`

**Interfaces:**
- Consumes: `RefDbProvider`, `register_provider` (Task 4).
- Produces: `class BiodbsProvider(RefDbProvider)` with `name = "biodbs"`. `fetch` lazily imports biodbs's amplicon-reference API via a module-level indirection `_load_biodbs_fetch()` (so tests can monkeypatch it) and calls `fetch(name, version, out_dir) -> (seqs_path, tax_path)`. A missing biodbs (ImportError) raises `MicrobiomeSuiteError`. The expected upstream biodbs entry point is `biodbs.amplicon.fetch_reference(name: str, version: str, out_dir: str) -> tuple[str, str]`; this is the interface biodbs must add.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_refdb_provider_biodbs.py
from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import get_provider
from microsuite.refdb.providers import biodbs as _biodbs  # noqa: F401  (force registration)
from microsuite.refdb.spec import RefDbSpec

FIXTURE = Path("src/microsuite/data/fixtures/refdb_mock")


def test_biodbs_is_default_provider() -> None:
    assert RefDbSpec(name="homd", version="15.22").provider == "biodbs"


def test_biodbs_fetch_uses_upstream_api(tmp_path: Path, monkeypatch) -> None:
    def fake_fetch(name, version, out_dir):
        return (str(FIXTURE / "source_a.fasta"), str(FIXTURE / "source_a.tax.tsv"))

    monkeypatch.setattr(_biodbs, "_load_biodbs_fetch", lambda: fake_fetch)
    provider = get_provider("biodbs")
    raw = provider.fetch(RefDbSpec(name="homd", version="15.22"), out_dir=tmp_path)
    assert raw.sequences.name == "source_a.fasta"
    assert raw.taxonomy.name == "source_a.tax.tsv"


def test_biodbs_missing_dependency_raises(tmp_path: Path, monkeypatch) -> None:
    def boom():
        raise ImportError("no biodbs")

    monkeypatch.setattr(_biodbs, "_load_biodbs_fetch", boom)
    provider = get_provider("biodbs")
    with pytest.raises(MicrobiomeSuiteError):
        provider.fetch(RefDbSpec(name="homd", version="15.22"), out_dir=tmp_path)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_refdb_provider_biodbs.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'microsuite.refdb.providers.biodbs'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/microsuite/refdb/providers/biodbs.py
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import register_provider
from microsuite.refdb.providers._base import RefDbProvider
from microsuite.refdb.spec import RawRefDb, RefDbSpec


def _load_biodbs_fetch() -> Callable[[str, str, str], tuple[str, str]]:
    from biodbs.amplicon import fetch_reference  # type: ignore[import-not-found]

    return fetch_reference


class BiodbsProvider(RefDbProvider):
    name = "biodbs"

    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            fetch = _load_biodbs_fetch()
        except ImportError as exc:
            raise MicrobiomeSuiteError(
                "The default 'biodbs' provider requires the biodbs package with "
                "amplicon-reference support. Install/upgrade biodbs, or pass a raw "
                "--classifier path, or use --provider rescript."
            ) from exc
        seqs, tax = fetch(spec.name, spec.version, str(out_dir))
        return RawRefDb(sequences=Path(seqs), taxonomy=Path(tax))


register_provider(BiodbsProvider())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_refdb_provider_biodbs.py -v`
Expected: PASS (3 passed)

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/refdb/providers/biodbs.py tests/test_refdb_provider_biodbs.py
git commit -m "feat(refdb): default biodbs provider with lazy-import guard"
```

---

### Task 7: Service orchestrator (fetch→merge→build→cache, and classifier resolution)

**Files:**
- Create: `src/microsuite/refdb/service.py`
- Modify: `src/microsuite/refdb/__init__.py` (import providers so both auto-register)
- Test: `tests/test_refdb_service.py`

**Interfaces:**
- Consumes: `RefDbSpec`, `RefDbSource` (Task 1); `RefDbRegistry` (Task 2); `merge_raw` (Task 3); `get_provider` (Task 4); `refdb_cache_dir` (Task 1).
- Produces:
  - `fetch_refdb(spec: RefDbSpec, build_target: str, *, force: bool = False, registry: RefDbRegistry | None = None, run_dir: Path | None = None, timeout: float | None = None) -> BuiltArtifact` — returns cached artifact when present and `not force`; otherwise fetches each source (or the spec itself if `sources` is empty), merges when >1 source, builds, records, and returns.
  - `resolve_classifier(value: str, *, registry: RefDbRegistry | None = None) -> Path` — parses `refdb:<name>@<version>[:<build>]` (default build `vsearch`) into a registry lookup; any other string is returned as `Path(value)`. An unresolved `refdb:` ref raises `MicrobiomeSuiteError`.

**`__init__.py` change:** add `from microsuite.refdb.providers import biodbs as _biodbs, rescript as _rescript  # noqa: F401` so importing `microsuite.refdb` registers both providers.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_refdb_service.py
from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import register_provider
from microsuite.refdb.providers._base import RefDbProvider
from microsuite.refdb.registry import RefDbRegistry
from microsuite.refdb.service import fetch_refdb, resolve_classifier
from microsuite.refdb.spec import RawRefDb, RefDbSource, RefDbSpec

FIXTURE = Path("src/microsuite/data/fixtures/refdb_mock")


class CountingProvider(RefDbProvider):
    name = "counting"

    def __init__(self) -> None:
        self.fetch_calls = 0

    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
        self.fetch_calls += 1
        tag = spec.name if (FIXTURE / f"{spec.name}.fasta").exists() else "source_a"
        return RawRefDb(
            sequences=FIXTURE / f"{tag}.fasta",
            taxonomy=FIXTURE / f"{tag}.tax.tsv",
        )


def test_fetch_then_cache_skips_second_fetch(tmp_path: Path) -> None:
    provider = CountingProvider()
    register_provider(provider)
    reg = RefDbRegistry(tmp_path / "cache")
    spec = RefDbSpec(name="source_a", version="1", provider="counting", build_targets=("vsearch",))

    first = fetch_refdb(spec, "vsearch", registry=reg)
    second = fetch_refdb(spec, "vsearch", registry=reg)

    assert first.path == second.path
    assert provider.fetch_calls == 1  # second call served from cache


def test_fetch_merges_multiple_sources(tmp_path: Path) -> None:
    register_provider(CountingProvider())
    reg = RefDbRegistry(tmp_path / "cache")
    spec = RefDbSpec(
        name="fomc-combined",
        version="20221029",
        provider="counting",
        sources=(RefDbSource("source_a", "1"), RefDbSource("source_b", "1")),
    )
    art = fetch_refdb(spec, "vsearch", registry=reg)
    ids = [
        line[1:].strip()
        for line in art.path.read_text().splitlines()
        if line.startswith(">")
    ]
    assert ids == ["seq1", "seq2", "seq3"]


def test_resolve_classifier_raw_path_passthrough(tmp_path: Path) -> None:
    raw = tmp_path / "my.qza"
    raw.write_text("x", encoding="utf-8")
    assert resolve_classifier(str(raw)) == raw


def test_resolve_classifier_registry_ref(tmp_path: Path) -> None:
    register_provider(CountingProvider())
    reg = RefDbRegistry(tmp_path / "cache")
    spec = RefDbSpec(name="source_a", version="1", provider="counting")
    art = fetch_refdb(spec, "vsearch", registry=reg)
    resolved = resolve_classifier("refdb:source_a@1", registry=reg)
    assert resolved == art.path


def test_resolve_classifier_unknown_ref_raises(tmp_path: Path) -> None:
    reg = RefDbRegistry(tmp_path / "cache")
    with pytest.raises(MicrobiomeSuiteError):
        resolve_classifier("refdb:ghost@9", registry=reg)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_refdb_service.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'microsuite.refdb.service'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/microsuite/refdb/service.py
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.build import merge_raw
from microsuite.refdb.paths import refdb_cache_dir
from microsuite.refdb.providers import get_provider
from microsuite.refdb.registry import RefDbRegistry
from microsuite.refdb.spec import BuiltArtifact, RefDbSpec


def _work_dir(registry: RefDbRegistry, spec: RefDbSpec, build_target: str) -> Path:
    work = registry.root / f"{spec.name}@{spec.version}" / build_target
    work.mkdir(parents=True, exist_ok=True)
    return work


def fetch_refdb(
    spec: RefDbSpec,
    build_target: str,
    *,
    force: bool = False,
    registry: RefDbRegistry | None = None,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> BuiltArtifact:
    registry = registry or RefDbRegistry(refdb_cache_dir())
    if not force:
        cached = registry.resolve(spec.name, spec.version, build_target)
        if cached is not None:
            return cached

    provider = get_provider(spec.provider)
    work = _work_dir(registry, spec, build_target)
    if spec.sources:
        raws = []
        for source in spec.sources:
            sub = replace(spec, name=source.name, version=source.version, sources=())
            sub_dir = work / "sources" / source.name
            sub_dir.mkdir(parents=True, exist_ok=True)
            raws.append(provider.fetch(sub, out_dir=sub_dir))
        raw = merge_raw(raws, out_dir=work / "merged")
    else:
        raw = provider.fetch(spec, out_dir=work / "fetch")

    artifact = provider.build(
        raw, build_target, out_dir=work / "build", run_dir=run_dir, timeout=timeout
    )
    registry.record(spec.name, spec.version, artifact, spec.provider)
    return artifact


def resolve_classifier(value: str, *, registry: RefDbRegistry | None = None) -> Path:
    if not value.startswith("refdb:"):
        return Path(value)
    body = value[len("refdb:") :]
    build_target = "vsearch"
    if ":" in body:
        body, build_target = body.rsplit(":", 1)
    if "@" not in body:
        raise MicrobiomeSuiteError(
            f"Malformed refdb reference '{value}'. Expected refdb:<name>@<version>[:<build>]."
        )
    name, version = body.split("@", 1)
    registry = registry or RefDbRegistry(refdb_cache_dir())
    art = registry.resolve(name, version, build_target)
    if art is None:
        raise MicrobiomeSuiteError(
            f"Reference DB '{name}@{version}:{build_target}' is not in the cache. "
            f"Run: microsuite refdb fetch {name} --version {version} --build {build_target}"
        )
    return art.path
```

Then update the package init:

```python
# src/microsuite/refdb/__init__.py
from __future__ import annotations

from microsuite.refdb.providers import biodbs as _biodbs  # noqa: F401
from microsuite.refdb.providers import rescript as _rescript  # noqa: F401
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_refdb_service.py -v`
Expected: PASS (5 passed)

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/refdb/service.py src/microsuite/refdb/__init__.py tests/test_refdb_service.py
git commit -m "feat(refdb): orchestrator with caching, source merge, and classifier resolution"
```

---

### Task 8: `microsuite refdb` CLI group and app wiring

**Files:**
- Create: `src/microsuite/cli/refdb_cmd.py`
- Modify: `src/microsuite/cli/app.py` (import + `add_typer`)
- Test: `tests/test_refdb_cli.py`

**Interfaces:**
- Consumes: `fetch_refdb` (Task 7), `RefDbSpec` (Task 1), `RefDbRegistry` (Task 2), `refdb_cache_dir` (Task 1), `VALID_BUILD_TARGETS` (Task 1).
- Produces: a Typer app with command `fetch` that prints the built artifact path. Signature: `microsuite refdb fetch NAME --version V [--provider biodbs|rescript] [--build vsearch|blast|qiime2] [--force]`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_refdb_cli.py
from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from microsuite.cli.app import app

FIXTURE = Path("src/microsuite/data/fixtures/refdb_mock")
runner = CliRunner()


def test_refdb_fetch_prints_artifact_path(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_REFDB_DIR", str(tmp_path / "cache"))

    def fake_fetch(name, version, out_dir):
        return (str(FIXTURE / "source_a.fasta"), str(FIXTURE / "source_a.tax.tsv"))

    monkeypatch.setattr(
        "microsuite.refdb.providers.biodbs._load_biodbs_fetch", lambda: fake_fetch
    )
    result = runner.invoke(
        app, ["refdb", "fetch", "homd", "--version", "15.22", "--build", "vsearch"]
    )
    assert result.exit_code == 0, result.output
    assert "reference.fasta" in result.output


def test_refdb_fetch_rejects_bad_build(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_REFDB_DIR", str(tmp_path / "cache"))
    result = runner.invoke(
        app, ["refdb", "fetch", "homd", "--version", "15.22", "--build", "bowtie"]
    )
    assert result.exit_code != 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_refdb_cli.py -v`
Expected: FAIL — `refdb` is not a registered command (nonzero exit / "No such command").

- [ ] **Step 3: Write minimal implementation**

```python
# src/microsuite/cli/refdb_cmd.py
from __future__ import annotations

from typing import Annotated

import typer

from microsuite.refdb.paths import VALID_BUILD_TARGETS
from microsuite.refdb.service import fetch_refdb
from microsuite.refdb.spec import RefDbSpec

app = typer.Typer(help="Fetch, build, and cache reference databases.", no_args_is_help=True)


@app.command("fetch")
def fetch(
    name: Annotated[str, typer.Argument(help="Reference DB name, e.g. homd, silva.")],
    version: Annotated[str, typer.Option("--version", help="DB version.")],
    provider: Annotated[
        str, typer.Option("--provider", help="Acquisition provider.")
    ] = "biodbs",
    build: Annotated[
        str, typer.Option("--build", help="Build target: vsearch, blast, or qiime2.")
    ] = "vsearch",
    force: Annotated[bool, typer.Option("--force", help="Rebuild even if cached.")] = False,
) -> None:
    if build not in VALID_BUILD_TARGETS:
        raise typer.BadParameter(f"--build must be one of: {', '.join(VALID_BUILD_TARGETS)}")
    spec = RefDbSpec(name=name, version=version, provider=provider, build_targets=(build,))
    artifact = fetch_refdb(spec, build, force=force)
    typer.echo(str(artifact.path))
```

Wire into the app:

```python
# src/microsuite/cli/app.py  — add to the import block
from microsuite.cli import (
    data_cmd,
    diffab_cmd,
    example_cmd,
    import_cmd,
    method_cmd,
    ml_cmd,
    network_cmd,
    ordination_cmd,
    qiime_cmd,
    refdb_cmd,
    viz_cmd,
    workflow_cmd,
    diversity_cmd,
)
```

```python
# src/microsuite/cli/app.py  — add inside _install_groups()
    app.add_typer(refdb_cmd.app, name="refdb")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_refdb_cli.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/cli/refdb_cmd.py src/microsuite/cli/app.py tests/test_refdb_cli.py
git commit -m "feat(refdb): microsuite refdb fetch CLI command"
```

---

### Task 9: Wire `refdb:` refs into `tax_classify --classifier`

**Files:**
- Modify: `src/microsuite/methods/tax_classify.py` (resolve classifier at entry)
- Modify: `src/microsuite/cli/method_taxonomy_cmd.py:34-42` (accept `str` and resolve)
- Test: `tests/test_refdb_tax_integration.py`

**Interfaces:**
- Consumes: `resolve_classifier` (Task 7).
- Produces: `tax_classify(..., classifier: Path | str | None, ...)` resolves a `refdb:` string to a cached path before dispatch; a raw path string/Path is unchanged.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_refdb_tax_integration.py
from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.tax_classify import tax_classify
from microsuite.refdb.providers import register_provider
from microsuite.refdb.providers._base import RefDbProvider
from microsuite.refdb.registry import RefDbRegistry
from microsuite.refdb.service import fetch_refdb
from microsuite.refdb.spec import RawRefDb, RefDbSpec

FIXTURE = Path("src/microsuite/data/fixtures/refdb_mock")


class FixtureProvider(RefDbProvider):
    name = "fixture"

    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
        return RawRefDb(
            sequences=FIXTURE / "source_a.fasta",
            taxonomy=FIXTURE / "source_a.tax.tsv",
        )


def test_tax_classify_resolves_refdb_ref(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_REFDB_DIR", str(tmp_path / "cache"))
    register_provider(FixtureProvider())
    reg = RefDbRegistry(tmp_path / "cache")
    art = fetch_refdb(
        RefDbSpec(name="source_a", version="1", provider="fixture"), "vsearch", registry=reg
    )

    captured: dict[str, Path | None] = {}

    def fake_qiime(*, rep_seqs, classifier, output, threads, force, run_dir, timeout):
        captured["classifier"] = classifier

    monkeypatch.setattr("microsuite.methods.tax_classify.tax_classify_qiime2", fake_qiime)
    rep = tmp_path / "rep.qza"
    rep.write_text("x", encoding="utf-8")

    tax_classify(
        backend="qiime2",
        rep_seqs=rep,
        classifier="refdb:source_a@1",
        output=tmp_path / "out.qza",
        threads=1,
        force=True,
    )
    assert captured["classifier"] == art.path


def test_tax_classify_unknown_ref_raises(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_REFDB_DIR", str(tmp_path / "cache"))
    rep = tmp_path / "rep.qza"
    rep.write_text("x", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError):
        tax_classify(
            backend="qiime2",
            rep_seqs=rep,
            classifier="refdb:ghost@9",
            output=tmp_path / "out.qza",
            threads=1,
            force=True,
        )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_refdb_tax_integration.py -v`
Expected: FAIL — `tax_classify` passes the string through unchanged, so `captured["classifier"]` is the string `"refdb:source_a@1"`, not `art.path` (assertion error); the second test does not raise.

- [ ] **Step 3: Write minimal implementation**

Read the current signature/head of `tax_classify` (`src/microsuite/methods/tax_classify.py:14-30`). Change the `classifier` parameter annotation to `Path | str | None` and resolve at the top of the function body, before `require_backend`:

```python
# src/microsuite/methods/tax_classify.py  — near the top of tax_classify(), after the docstring/first line
    from microsuite.refdb.service import resolve_classifier

    if classifier is not None:
        classifier = resolve_classifier(str(classifier))
```

Place this block as the first statements inside `tax_classify(...)` so every backend receives an already-resolved `Path`. (Import is local to avoid any import cycle between `methods` and `refdb`.)

Update the CLI option type so a `refdb:` string is accepted (`src/microsuite/cli/method_taxonomy_cmd.py:33-42`):

```python
        classifier: Annotated[
            str | None,
            typer.Option(
                "--classifier",
                help=(
                    "QIIME classifier / Kraken2 / MetaPhlAn / EMU database as a path, "
                    "OR a cached reference as 'refdb:<name>@<version>[:<build>]'."
                ),
            ),
        ] = None,
```

The command body already forwards `classifier=classifier` to `tax_classify`; passing a `str` is now correct because `tax_classify` resolves it. Any backend that needs a `Path` receives the resolved `Path` from `resolve_classifier`.

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_refdb_tax_integration.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Run the full refdb suite + a broad sanity check**

Run: `uv run pytest tests/test_refdb_*.py tests/test_cli.py -q`
Expected: PASS (no regressions in CLI registration).

- [ ] **Step 6: Commit**

```bash
git add src/microsuite/methods/tax_classify.py src/microsuite/cli/method_taxonomy_cmd.py tests/test_refdb_tax_integration.py
git commit -m "feat(refdb): resolve refdb: references in tax_classify --classifier"
```

---

## Self-Review

**Spec coverage:**
- Provider interface + two providers (biodbs default, rescript optional) → Tasks 4, 5, 6. ✓
- JSON registry with checksum/caching/validation → Task 2; caching behavior asserted in Task 7. ✓
- Build layer (vsearch/blast/qiime2) + FOMC merge → Task 3; merge orchestration in Task 7. ✓
- `$MICROSUITE_REFDB_DIR` cache dir → Task 1. ✓
- CLI `refdb fetch` with `--provider` defaulting to biodbs → Task 8. ✓
- `--classifier` backward-compatible + `refdb:` refs → Task 9 (raw-path passthrough asserted in Task 7 and Task 9). ✓
- Offline tests, external tools argv-asserted → every task. ✓
- biodbs upstream interface documented (`biodbs.amplicon.fetch_reference`) → Task 6 interface block; this is the contract biodbs must implement. ✓
- `MicrobiomeSuiteError` on all failure paths → Tasks 2, 3, 5, 6, 7. ✓

**Deferred to biodbs repo (out of this plan's scope, noted in spec):** implementing `biodbs.amplicon.fetch_reference` for SILVA/NCBI/GTDB/UNITE/GG2/HOMD/MOMD. The microsuite side is fully testable now via the monkeypatched `_load_biodbs_fetch` seam.

**Type consistency check:** `RawRefDb(sequences, taxonomy, qza)`, `BuiltArtifact(path, build_target, checksum)`, `fetch(spec, out_dir)`, `build(raw, build_target, out_dir, ...)`, `fetch_refdb(spec, build_target, ...)`, `resolve_classifier(value, ...)` — names/signatures match across Tasks 1, 3, 4, 5, 6, 7, 9. ✓

**Placeholder scan:** no TBD/TODO; every code step shows full code. ✓
