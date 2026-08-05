# Batch Effect Correction — Phases 1–2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship `microsuite batch correct` with five R backends (`mmuphin`, `combat-seq`, `conqur`, `plsda-batch`, `metadict`), a three-valued output-scale contract that stops corrected tables from being silently misused downstream, and a real-execution smoke test that asserts the batch effect shrank and the biological signal survived.

**Architecture:** A new `microsuite.batch` package holds one capability record and one R script per backend. The container/mount/provenance machinery currently trapped in `diffab/_runner.py` moves to `runtime/r_backend.py` and is shared by both packages. Every corrected table carries `uns["microsuite"]["value_type"]`, and five existing call sites refuse tables whose scale they cannot consume.

**Tech Stack:** Python 3.11+, anndata, pandas, numpy, typer, pytest; R 4.3 via micromamba containers (Bioconductor MMUPHin and sva; GitHub ConQuR, PLSDAbatch, MetaDICT).

**Source spec:** `docs/superpowers/specs/2026-08-05-batch-effect-correction-design.md`

**Not in this plan:** `debias-m` (spec phase 3) and `microsuite batch diagnose` (spec phase 4). They get their own plan.

## Global Constraints

- Target release **0.3.0**, which is in `CHANGELOG.md` but not tagged. Fold entries into the existing unreleased 0.3.0 section; do not open a 0.4.0 section.
- `value_type` is exactly one of `counts`, `relative`, `clr`. No other value is ever written or accepted.
- An absent `uns["microsuite"]["value_type"]` **always** means "skip the check". Pre-0.3.0 tables must behave exactly as they do today.
- Every R backend script takes exactly four positional arguments: `counts.tsv metadata.tsv params.json corrected.tsv`. Options travel in the JSON, never positionally.
- Every container Dockerfile carries `org.opencontainers.image.title`, `org.opencontainers.image.description`, and a `# Expected commands:` comment — `tests/test_container_skeletons.py` asserts all three.
- Every container runs a build-time smoke through the real backend script on real toy data and fails the build on an empty result. A package-import check is not sufficient.
- The three GitHub-sourced R packages (ConQuR, PLSDAbatch, MetaDICT) pin to a **commit SHA**, never a branch.
- `ruff format`, `ruff check`, and `ty check` must pass. Line length 100. `from __future__ import annotations` at the top of every new Python module.
- Commit messages carry **no** `Co-Authored-By:` and no `Claude-Session:` trailer.
- Run `uv sync --all-extras --locked` before `ty check`, or you will see a spurious `biodbs` unresolved-import error that CI does not have.
- `pyproject.toml` sets `addopts = "-q"`, so pytest prints no "N passed" line. Use `--junitxml` and read the XML when you need a reliable count.

---

## File Structure

**Created:**

| File | Responsibility |
|---|---|
| `src/microsuite/runtime/r_backend.py` | The one place that knows how to run an R script locally or in a container, mount its paths, and write the image-digest sidecar. Shared by `diffab` and `batch`. |
| `src/microsuite/batch/__init__.py` | Package marker. |
| `src/microsuite/batch/value_type.py` | The output-scale contract: write it, read it, enforce it. No knowledge of any backend. |
| `src/microsuite/batch/backends.py` | One `BatchBackend` record per backend — script name, R package, emitted `value_type`, covariate and target support. The dispatch reads this table instead of branching per backend. |
| `src/microsuite/batch/correct.py` | Marshals an AnnData to TSV, invokes the backend, reads the corrected table back, and rebuilds an aligned AnnData with provenance. |
| `src/microsuite/batch/r/{mmuphin,combat_seq,conqur,plsda_batch,metadict}.R` | One script per backend. |
| `src/microsuite/methods/batch_correct.py` | Public API and `SUPPORTED_BACKENDS`. |
| `src/microsuite/cli/batch_cmd.py` | `microsuite batch correct`. |
| `containers/r-batch-{mmuphin,combatseq,conqur,plsdabatch,metadict}/` | One image per backend, each with a `smoke/` directory. |
| `tests/test_batch_value_type.py`, `tests/test_batch_backends.py`, `tests/test_batch_correct.py`, `tests/test_batch_cli.py` | Unit tests. |
| `tests/integration/test_batch_correct_smoke.py` | Real-execution smoke, gated. |
| `docs/batch_correction.md` | User-facing guidance. |

**Modified:**

| File | Change |
|---|---|
| `src/microsuite/diffab/_runner.py` | Becomes a thin shim over `runtime/r_backend.py`. Public signature unchanged. |
| `src/microsuite/runtime/container.py:215` | Gains `resolve_batch_image` beside `resolve_diffab_image`. |
| `src/microsuite/methods/diff_abundance.py` | `ancombc` and `aldex2` require counts. |
| `src/microsuite/methods/rarefy.py` | `rarefy_native` requires counts. |
| `src/microsuite/methods/normalize.py` | `normalize_native` guards `relative`, `total-sum`, `clr`. |
| `src/microsuite/cli/app.py:51-64` | Registers the `batch` sub-app. |
| `src/microsuite/cli/_method_api.py:66` | Adds `batch_correct` to `METHOD_BACKENDS`. |
| `tests/test_diffab_runner.py` | Monkeypatch target moves with the code. |
| `tests/test_container_skeletons.py:15` | Five new expected images. |
| `.github/workflows/docker.yml:117` | Five new heavy matrix entries. |
| `docs/methods.md`, `CHANGELOG.md` | Documentation. |

**One decision the spec left implicit.** The spec's guard table names `normalize --method relative` but not `--method total-sum`. They are the same operation with a scale factor, so `total-sum` gets the identical guard. `--method prevalence-filter` gets no guard: filtering by prevalence is valid at any scale.

---

### Task 1: Extract the R runner out of `diffab/`

`diffab/_runner.py` holds the bind-mount layout, caller-UID execution, and digest sidecar that any containerized R backend needs, but it is hardcoded to scripts in `microsuite.diffab.r` and to `resolve_diffab_image`. The batch package needs the same 110 lines. Move them once.

Two details that make this more than a file move:

1. Backend names and script names diverge. The backend is `combat-seq`; the script is `combat_seq.R`. The current code derives the script path from the backend name. Script name becomes a separate parameter.
2. `tests/test_diffab_runner.py` monkeypatches `_runner.run_command` and `_runner.shutil`. Once the body moves, those patches bind to a module that no longer executes anything, and the tests would pass while testing nothing. The patch target moves with the code in this same task.

**Files:**
- Create: `src/microsuite/runtime/r_backend.py`
- Modify: `src/microsuite/diffab/_runner.py` (replace all 111 lines)
- Modify: `src/microsuite/runtime/container.py:215-222` (add `resolve_batch_image` after `resolve_diffab_image`)
- Test: `tests/test_diffab_runner.py` (repoint monkeypatches), `tests/test_batch_image_resolution.py` (create)

**Interfaces:**
- Consumes: `microsuite.runtime.container.{PathMapper, build_container_command, host_user_spec, require_engine, resolve_diffab_image, resolve_image_digest}`, `microsuite.runtime.runner.{CommandLog, run_command}` — all existing.
- Produces:
  - `runtime.r_backend.invoke_r_script(*, backend: str, script_package: str, script_name: str, resolve_image: Callable[[str, str | None], str], positional: list[str | Path], runtime: str = "local", image: str | None = None, engine: str = "docker", run_dir: Path | None = None, timeout: float | None = None, log: CommandLog, local_missing_message: str) -> None`
  - `runtime.container.resolve_batch_image(backend: str, override: str | None) -> str`
  - `diffab._runner.invoke_r_backend(...)` — signature unchanged from today.

- [ ] **Step 1: Write the failing test for image resolution**

Create `tests/test_batch_image_resolution.py`:

```python
from __future__ import annotations

from microsuite.runtime.container import resolve_batch_image


def test_default_image_is_per_backend() -> None:
    assert resolve_batch_image("mmuphin", None) == (
        "ghcr.io/qwerty239qwe/microsuite/r-batch-mmuphin:latest"
    )


def test_explicit_override_wins_over_env(monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_R_BATCH_MMUPHIN_IMAGE", "from-env:1")
    assert resolve_batch_image("mmuphin", "explicit:2") == "explicit:2"


def test_env_override_is_uppercased_and_underscored(monkeypatch) -> None:
    # The backend is 'combat-seq'; the env var cannot contain a hyphen.
    monkeypatch.setenv("MICROSUITE_R_BATCH_COMBAT_SEQ_IMAGE", "from-env:1")
    assert resolve_batch_image("combat-seq", None) == "from-env:1"
```

- [ ] **Step 2: Write the failing test for script-name/backend divergence**

Append to `tests/test_diffab_runner.py`:

```python
from microsuite.runtime import r_backend
from microsuite.runtime.r_backend import invoke_r_script


def test_script_name_can_differ_from_backend_name(tmp_path: Path, monkeypatch) -> None:
    captured: dict = {}
    monkeypatch.setattr(r_backend, "run_command", lambda command, **kw: captured.update(cmd=command))
    monkeypatch.setattr(r_backend.shutil, "which", lambda name: "/usr/bin/Rscript")
    invoke_r_script(
        backend="combat-seq",
        script_package="microsuite.diffab.r",
        script_name="ancombc",  # any script that exists; the point is the two names differ
        resolve_image=lambda backend, override: "unused",
        positional=[tmp_path / "counts.tsv", tmp_path / "out.tsv"],
        runtime="local",
        log=CommandLog(task="batch_correct", backend="combat-seq"),
        local_missing_message="need R",
    )
    assert captured["cmd"][1].endswith("ancombc.R")
```

- [ ] **Step 3: Run both to verify they fail**

Run: `uv run pytest tests/test_batch_image_resolution.py tests/test_diffab_runner.py -v`
Expected: FAIL — `ImportError: cannot import name 'resolve_batch_image'` and `ModuleNotFoundError: microsuite.runtime.r_backend`.

- [ ] **Step 4: Add `resolve_batch_image`**

In `src/microsuite/runtime/container.py`, beside the existing prefix constants near line 18:

```python
DEFAULT_BATCH_IMAGE_PREFIX = "ghcr.io/qwerty239qwe/microsuite/r-batch-"
_BATCH_IMAGE_ENV_PREFIX = "MICROSUITE_R_BATCH_"
```

and after `resolve_diffab_image`:

```python
def resolve_batch_image(backend: str, override: str | None) -> str:
    """Resolve the per-backend r-batch image: override, then env, then default."""
    if override:
        return override
    env_name = f"{_BATCH_IMAGE_ENV_PREFIX}{backend.upper().replace('-', '_')}_IMAGE"
    env = os.environ.get(env_name)
    if env:
        return env
    return f"{DEFAULT_BATCH_IMAGE_PREFIX}{backend}:latest"
```

- [ ] **Step 5: Create `runtime/r_backend.py`**

Move the body of `diffab/_runner.py` verbatim, changing only what the two new parameters require:

```python
"""Shared local/container runner for microsuite's R backends.

Callers build an ordered positional argument list for their ``.R`` script and
hand it here. Convention: the last ``Path`` is the output (its parent is
bind-mounted read-write); every earlier ``Path`` is an input (its parent is
bind-mounted read-only); ``str`` items pass through verbatim. Container runs
execute the per-backend image as the caller's UID/GID (writable outputs) and
record the resolved image + digest in ``<output-dir>/<backend>_container.json``.

``script_name`` is separate from ``backend`` because the two diverge: the
``combat-seq`` backend runs ``combat_seq.R``.
"""

from __future__ import annotations

import json
import shutil
from collections.abc import Callable
from importlib.resources import files
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.runtime.container import (
    PathMapper,
    build_container_command,
    host_user_spec,
    require_engine,
    resolve_image_digest,
)
from microsuite.runtime.runner import CommandLog, run_command


def invoke_r_script(
    *,
    backend: str,
    script_package: str,
    script_name: str,
    resolve_image: Callable[[str, str | None], str],
    positional: list[str | Path],
    runtime: str = "local",
    image: str | None = None,
    engine: str = "docker",
    run_dir: Path | None = None,
    timeout: float | None = None,
    log: CommandLog,
    local_missing_message: str,
) -> None:
    if runtime not in ("local", "docker"):
        raise MicrobiomeSuiteError(
            f"Unsupported --runtime '{runtime}' for {backend}; choose 'local' or 'docker'."
        )

    if runtime == "local":
        rscript = shutil.which("Rscript")
        if rscript is None:
            raise MicrobiomeSuiteError(local_missing_message)
        script = files(script_package).joinpath(f"{script_name}.R")
        command = [rscript, str(script), *[str(arg) for arg in positional]]
        run_command(
            command,
            failure_message=f"{backend} failed.",
            run_dir=run_dir,
            log=log,
            timeout=timeout,
        )
        return

    resolved_image = resolve_image(backend, image)
    require_engine(engine)
    paths = [arg for arg in positional if isinstance(arg, Path)]
    if not paths:
        raise MicrobiomeSuiteError(f"{backend} requires at least one file argument.")
    output_path = paths[-1]

    mapper = PathMapper()
    mountpoints: dict[Path, str] = {}

    def _mount(host_dir: Path, mode: str) -> None:
        resolved = host_dir.resolve()
        if resolved not in mountpoints:
            mountpoints[resolved] = f"/mnt/d{len(mountpoints)}"
        mapper.add_dir(host_dir, mode, mountpoints[resolved])

    for path in paths[:-1]:
        _mount(path.parent, "ro")
    _mount(output_path.parent, "rw")

    inner = [f"/opt/microsuite/{script_name}.R"]
    for arg in positional:
        inner.append(mapper.to_container(arg) if isinstance(arg, Path) else arg)
    command = build_container_command(
        inner, resolved_image, mapper.mounts(), engine=engine, user=host_user_spec()
    )
    run_command(
        command,
        failure_message=f"{backend} failed.",
        run_dir=run_dir,
        log=log,
        timeout=timeout,
    )

    sidecar = output_path.parent / f"{backend}_container.json"
    sidecar.write_text(
        json.dumps(
            {
                "runtime": "docker",
                "engine": engine,
                "image": resolved_image,
                "digest": resolve_image_digest(engine, resolved_image),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
```

- [ ] **Step 6: Replace `diffab/_runner.py` with the shim**

```python
"""Compatibility shim: the diffab R runner now lives in ``runtime.r_backend``.

Kept so every existing diffab caller and its tests are unaffected by the move.
New backends should call ``invoke_r_script`` directly.
"""

from __future__ import annotations

from pathlib import Path

from microsuite.runtime.container import resolve_diffab_image
from microsuite.runtime.r_backend import invoke_r_script
from microsuite.runtime.runner import CommandLog


def invoke_r_backend(
    *,
    backend: str,
    positional: list[str | Path],
    runtime: str = "local",
    image: str | None = None,
    engine: str = "docker",
    run_dir: Path | None = None,
    timeout: float | None = None,
    log: CommandLog,
    local_missing_message: str,
) -> None:
    invoke_r_script(
        backend=backend,
        script_package="microsuite.diffab.r",
        script_name=backend,
        resolve_image=resolve_diffab_image,
        positional=positional,
        runtime=runtime,
        image=image,
        engine=engine,
        run_dir=run_dir,
        timeout=timeout,
        log=log,
        local_missing_message=local_missing_message,
    )
```

- [ ] **Step 7: Repoint the existing monkeypatches**

In `tests/test_diffab_runner.py`, every `monkeypatch.setattr(_runner, "run_command", ...)` becomes `monkeypatch.setattr(r_backend, "run_command", ...)`, and every `monkeypatch.setattr(_runner.shutil, ...)` becomes `monkeypatch.setattr(r_backend.shutil, ...)`. Add `from microsuite.runtime import r_backend` at the top. Leave the `invoke_r_backend` call sites alone — the point is that the public signature did not change.

- [ ] **Step 8: Prove the repointed patches actually bite**

Temporarily change `run_command` in `r_backend.py` to `raise AssertionError("unpatched")`. Run `uv run pytest tests/test_diffab_runner.py -v`. If any test still passes, its patch is not reaching the executing code — fix it. Revert the deliberate break.

This step exists because a monkeypatch aimed at a module that no longer runs anything fails silently and leaves a green suite testing nothing.

- [ ] **Step 9: Run the full affected suite**

Run: `uv run pytest tests/test_diffab_runner.py tests/test_batch_image_resolution.py tests/test_diffab_ancombc.py tests/test_diff_abundance_method.py -v`
Expected: all PASS.

- [ ] **Step 10: Lint and commit**

```bash
uv run ruff format src/microsuite/runtime/r_backend.py src/microsuite/diffab/_runner.py src/microsuite/runtime/container.py tests/test_batch_image_resolution.py tests/test_diffab_runner.py
uv run ruff check src tests
git add src/microsuite/runtime/r_backend.py src/microsuite/diffab/_runner.py src/microsuite/runtime/container.py tests/test_batch_image_resolution.py tests/test_diffab_runner.py
git commit -m "refactor(runtime): share the R backend runner between diffab and batch

Parametrizes the script package, script name, and image resolver so a second
backend family can reuse the mount, UID, and digest-sidecar logic instead of
copying it. diffab/_runner.py keeps its public signature."
```

---

### Task 2: The `value_type` contract

The six backends return data on three different scales. Nothing downstream can tell them apart by inspection: a CLR matrix and a relative-abundance matrix are both float arrays of the right shape. This module is where the scale becomes explicit.

**Files:**
- Create: `src/microsuite/batch/__init__.py`, `src/microsuite/batch/value_type.py`
- Test: `tests/test_batch_value_type.py`

**Interfaces:**
- Consumes: `microsuite._errors.MicrobiomeSuiteError`.
- Produces:
  - `VALUE_TYPES: tuple[str, ...] = ("counts", "relative", "clr")`
  - `record_batch_correction(adata: ad.AnnData, *, value_type: str, backend: str, batch: str, covariates: list[str], target: str | None) -> None`
  - `read_value_type(adata: ad.AnnData) -> str | None`
  - `require_value_types(adata: ad.AnnData, allowed: tuple[str, ...], *, operation: str) -> None`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_batch_value_type.py`:

```python
from __future__ import annotations

import anndata as ad
import numpy as np
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.batch.value_type import (
    read_value_type,
    record_batch_correction,
    require_value_types,
)


def _adata() -> ad.AnnData:
    return ad.AnnData(np.array([[1.0, 2.0], [3.0, 4.0]]))


def test_absent_key_skips_the_check() -> None:
    # Any table written before 0.3.0 has no key at all. Behaviour must not change.
    require_value_types(_adata(), ("counts",), operation="rarefy")


def test_microsuite_key_without_value_type_skips_the_check() -> None:
    adata = _adata()
    adata.uns["microsuite"] = {"source": "tsv"}  # io/tsv.py writes this shape
    require_value_types(adata, ("counts",), operation="rarefy")


def test_matching_value_type_passes() -> None:
    adata = _adata()
    record_batch_correction(
        adata, value_type="counts", backend="combat-seq", batch="run_id",
        covariates=[], target=None,
    )
    require_value_types(adata, ("counts",), operation="rarefy")


def test_mismatched_value_type_names_the_producing_backend() -> None:
    adata = _adata()
    record_batch_correction(
        adata, value_type="clr", backend="plsda-batch", batch="run_id",
        covariates=[], target="disease",
    )
    with pytest.raises(MicrobiomeSuiteError) as excinfo:
        require_value_types(adata, ("counts",), operation="diff_abundance --backend ancombc")
    message = str(excinfo.value)
    assert "diff_abundance --backend ancombc" in message
    assert "clr" in message
    assert "plsda-batch" in message


def test_provenance_is_recorded_alongside_the_scale() -> None:
    adata = _adata()
    record_batch_correction(
        adata, value_type="relative", backend="mmuphin", batch="run_id",
        covariates=["sex", "age"], target=None,
    )
    assert read_value_type(adata) == "relative"
    provenance = adata.uns["microsuite"]["batch_correct"]
    assert provenance["backend"] == "mmuphin"
    assert provenance["batch"] == "run_id"
    assert provenance["covariates"] == ["sex", "age"]
    assert provenance["target"] is None


def test_recording_preserves_existing_microsuite_keys() -> None:
    adata = _adata()
    adata.uns["microsuite"] = {"source": "tsv"}
    record_batch_correction(
        adata, value_type="counts", backend="conqur", batch="run_id",
        covariates=[], target=None,
    )
    assert adata.uns["microsuite"]["source"] == "tsv"
    assert adata.uns["microsuite"]["value_type"] == "counts"


def test_unknown_value_type_is_rejected_at_write_time() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="abundance"):
        record_batch_correction(
            _adata(), value_type="abundance", backend="mmuphin", batch="run_id",
            covariates=[], target=None,
        )
```

The last test matters: `abundance` was the 2026-08-02 spec's term for what is now `relative`. Rejecting it at write time stops the old vocabulary leaking back in.

- [ ] **Step 2: Run to verify they fail**

Run: `uv run pytest tests/test_batch_value_type.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'microsuite.batch'`.

- [ ] **Step 3: Create the package and the module**

`src/microsuite/batch/__init__.py`:

```python
"""Batch effect correction backends and the output-scale contract."""
```

`src/microsuite/batch/value_type.py`:

```python
"""The output-scale contract for batch-corrected tables.

Batch correction backends disagree about what they return: ComBat-seq and
ConQuR emit integer counts, MMUPHin and MetaDICT emit relative abundances, and
PLSDA-batch emits CLR log-ratios. Downstream, ANCOM-BC and ALDEx2 require
counts, ``rarefy`` requires counts, and ``normalize`` must not transform data
that is already transformed.

Nothing about a float matrix reveals which of the three it is, so the producing
command records it and the consuming commands assert on it. An absent key means
"unknown" and is always allowed: every table written before 0.3.0 lacks it, and
this contract must never change the behaviour of existing pipelines.
"""

from __future__ import annotations

import anndata as ad

from microsuite._errors import MicrobiomeSuiteError

VALUE_TYPES: tuple[str, ...] = ("counts", "relative", "clr")


def record_batch_correction(
    adata: ad.AnnData,
    *,
    value_type: str,
    backend: str,
    batch: str,
    covariates: list[str],
    target: str | None,
) -> None:
    """Stamp the corrected table with its scale and how it was produced."""
    if value_type not in VALUE_TYPES:
        raise MicrobiomeSuiteError(
            f"Unknown value_type '{value_type}'. Choose one of: {', '.join(VALUE_TYPES)}"
        )
    info = adata.uns.get("microsuite")
    if not isinstance(info, dict):
        info = {}
    info = dict(info)
    info["value_type"] = value_type
    info["batch_correct"] = {
        "backend": backend,
        "batch": batch,
        "covariates": list(covariates),
        "target": target,
    }
    adata.uns["microsuite"] = info


def read_value_type(adata: ad.AnnData) -> str | None:
    """Return the recorded scale, or None when the table does not declare one."""
    info = adata.uns.get("microsuite")
    if not isinstance(info, dict):
        return None
    value_type = info.get("value_type")
    return value_type if isinstance(value_type, str) else None


def require_value_types(
    adata: ad.AnnData, allowed: tuple[str, ...], *, operation: str
) -> None:
    """Raise when the table declares a scale that ``operation`` cannot consume."""
    value_type = read_value_type(adata)
    if value_type is None or value_type in allowed:
        return
    info = adata.uns.get("microsuite")
    provenance = info.get("batch_correct") if isinstance(info, dict) else None
    origin = ""
    if isinstance(provenance, dict) and provenance.get("backend"):
        origin = (
            f" It was produced by 'batch correct --backend {provenance['backend']}', "
            f"which emits '{value_type}'."
        )
    raise MicrobiomeSuiteError(
        f"{operation} requires a table of type {' or '.join(allowed)}, "
        f"but this table is '{value_type}'.{origin} "
        f"Use a backend that accepts '{value_type}', or correct a different way."
    )
```

- [ ] **Step 4: Run to verify they pass**

Run: `uv run pytest tests/test_batch_value_type.py -v`
Expected: all PASS.

- [ ] **Step 5: Lint and commit**

```bash
uv run ruff format src/microsuite/batch tests/test_batch_value_type.py
uv run ruff check src tests
git add src/microsuite/batch tests/test_batch_value_type.py
git commit -m "feat(batch): add the counts/relative/clr output-scale contract

Batch correction backends return data on three different scales that are
indistinguishable by inspection. The producing command records the scale; the
consuming commands assert on it. An absent key always means unknown and is
always allowed, so pre-0.3.0 tables are unaffected."
```

---

### Task 3: Enforce the contract at the five call sites

The contract is worthless until something reads it. Feeding a CLR matrix to ANCOM-BC produces a complete, plausible, wrong table — the same failure shape as the 0.2.1 depth-confounding defect.

**Files:**
- Modify: `src/microsuite/methods/diff_abundance.py:60-85`
- Modify: `src/microsuite/methods/rarefy.py:32-41`
- Modify: `src/microsuite/methods/normalize.py:62-96`
- Test: `tests/test_batch_value_type_guards.py` (create)

**Interfaces:**
- Consumes: `microsuite.batch.value_type.require_value_types` from Task 2.
- Produces: no new public symbols. Behaviour change only.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_batch_value_type_guards.py`:

```python
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.batch.value_type import record_batch_correction
from microsuite.io.h5ad import write_h5ad
from microsuite.methods.normalize import normalize_native
from microsuite.methods.rarefy import rarefy_native


def _table(value_type: str | None, backend: str = "mmuphin") -> ad.AnnData:
    adata = ad.AnnData(np.array([[10.0, 30.0], [20.0, 20.0]]))
    adata.obs_names = ["s1", "s2"]
    adata.var_names = ["f1", "f2"]
    if value_type is not None:
        record_batch_correction(
            adata, value_type=value_type, backend=backend, batch="run_id",
            covariates=[], target=None,
        )
    return adata


@pytest.mark.parametrize("value_type", ["relative", "clr"])
def test_rarefy_rejects_non_counts(value_type: str) -> None:
    with pytest.raises(MicrobiomeSuiteError, match="rarefy"):
        rarefy_native(_table(value_type), depth=10)


def test_rarefy_accepts_counts_and_unmarked_tables() -> None:
    rarefy_native(_table("counts"), depth=10)
    rarefy_native(_table(None), depth=10)


@pytest.mark.parametrize("method", ["relative", "total-sum"])
@pytest.mark.parametrize("value_type", ["relative", "clr"])
def test_normalize_rejects_already_scaled_input(method: str, value_type: str) -> None:
    with pytest.raises(MicrobiomeSuiteError, match="normalize"):
        normalize_native(_table(value_type), method=method)


def test_normalize_clr_rejects_clr_but_accepts_relative() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="normalize"):
        normalize_native(_table("clr"), method="clr")
    normalize_native(_table("relative"), method="clr")


def test_prevalence_filter_is_scale_agnostic() -> None:
    # Filtering by prevalence is valid at any scale, so it carries no guard.
    normalize_native(_table("clr"), method="prevalence-filter")


@pytest.mark.parametrize("backend", ["ancombc", "aldex2"])
def test_diff_abundance_count_backends_reject_clr(
    backend: str, tmp_path: Path, monkeypatch
) -> None:
    from microsuite.methods import diff_abundance as module

    table = tmp_path / "corrected.h5ad"
    write_h5ad(_table("clr", backend="plsda-batch"), table)
    monkeypatch.setattr(
        module, "run_ancombc", lambda *a, **kw: pytest.fail("backend must not be invoked")
    )
    monkeypatch.setattr(
        module,
        "run_r_diffab_backend",
        lambda *a, **kw: pytest.fail("backend must not be invoked"),
    )
    with pytest.raises(MicrobiomeSuiteError, match="plsda-batch"):
        module.diff_abundance(
            backend=backend, table=table, group="g", output=tmp_path / "out.tsv"
        )


@pytest.mark.parametrize("backend", ["maaslin2", "lefse"])
def test_diff_abundance_internally_normalizing_backends_do_not_check(
    backend: str, tmp_path: Path, monkeypatch
) -> None:
    from microsuite.methods import diff_abundance as module

    invoked: dict = {}
    table = tmp_path / "corrected.h5ad"
    write_h5ad(_table("clr", backend="plsda-batch"), table)
    monkeypatch.setattr(
        module, "run_r_diffab_backend", lambda *a, **kw: invoked.update(ran=True)
    )
    module.diff_abundance(
        backend=backend, table=table, group="g", output=tmp_path / "out.tsv"
    )
    assert invoked["ran"] is True
```

- [ ] **Step 2: Run to verify they fail**

Run: `uv run pytest tests/test_batch_value_type_guards.py -v`
Expected: FAIL — no guard exists, so every rejection test reports `DID NOT RAISE`.

- [ ] **Step 3: Guard `rarefy_native`**

In `src/microsuite/methods/rarefy.py`, add the import and put the check first, before the depth check:

```python
from microsuite.batch.value_type import require_value_types
```

```python
def rarefy_native(adata: ad.AnnData, *, depth: int, seed: int = 0) -> ad.AnnData:
    # Subsampling reads from a matrix that no longer holds read counts is not
    # meaningful, and produces a plausible table rather than an error.
    require_value_types(adata, ("counts",), operation="rarefy")
    if depth <= 0:
        raise MicrobiomeSuiteError("Rarefaction depth must be greater than zero.")
```

- [ ] **Step 4: Guard `normalize_native`**

In `src/microsuite/methods/normalize.py`, add the import and insert after the method-validation block, before `counts = dense_counts(adata)`:

```python
from microsuite.batch.value_type import require_value_types
```

```python
    # Re-scaling data that is already relative, or CLR-transforming data that is
    # already CLR, returns a full table computed twice over.
    if method in ("relative", "total-sum"):
        require_value_types(adata, ("counts",), operation=f"normalize --method {method}")
    elif method == "clr":
        require_value_types(adata, ("counts", "relative"), operation="normalize --method clr")
```

`prevalence-filter` deliberately has no branch.

- [ ] **Step 5: Guard the count-requiring diff_abundance backends**

In `src/microsuite/methods/diff_abundance.py`, add:

```python
from microsuite.batch.value_type import require_value_types

COUNT_REQUIRING_BACKENDS = ("ancombc", "aldex2")
```

In the `backend in R_BACKENDS` branch, after `adata = read_h5ad(...)`:

```python
        adata = read_h5ad(ensure_input(table))
        if backend in COUNT_REQUIRING_BACKENDS:
            require_value_types(
                adata, ("counts",), operation=f"diff_abundance --backend {backend}"
            )
```

and in the trailing native-ancombc path, after its `read_h5ad`:

```python
    adata = read_h5ad(ensure_input(table))
    require_value_types(adata, ("counts",), operation="diff_abundance --backend ancombc")
```

`maaslin2` and `lefse` fall through unguarded because both normalize internally.

- [ ] **Step 6: Run to verify they pass**

Run: `uv run pytest tests/test_batch_value_type_guards.py -v`
Expected: all PASS.

- [ ] **Step 7: Confirm each guard actually fires**

For each of the three modified files in turn, comment out its `require_value_types` call, run `uv run pytest tests/test_batch_value_type_guards.py -v`, and confirm the corresponding tests **fail**. Restore the call.

A guard whose test passes with the guard removed is testing nothing. Do all three; do not sample.

- [ ] **Step 8: Confirm nothing regressed**

Run: `uv run pytest tests/test_methods.py tests/test_diff_abundance_method.py tests/test_alpha.py tests/test_beta.py -v`
Expected: all PASS. Every existing table lacks the key, so every existing path skips every check.

- [ ] **Step 9: Lint and commit**

```bash
uv run ruff format src/microsuite/methods/rarefy.py src/microsuite/methods/normalize.py src/microsuite/methods/diff_abundance.py tests/test_batch_value_type_guards.py
uv run ruff check src tests
git add src/microsuite/methods/rarefy.py src/microsuite/methods/normalize.py src/microsuite/methods/diff_abundance.py tests/test_batch_value_type_guards.py
git commit -m "feat(methods): refuse tables whose scale the operation cannot consume

ancombc, aldex2, and rarefy now require counts; normalize refuses input that is
already relative or already CLR. maaslin2 and lefse are unguarded because both
normalize internally. Tables without a declared scale are unaffected."
```

---

### Task 4: The backend capability table

Five backends differ along four axes: whether they accept covariates, whether they need an outcome label, what scale they emit, and which R package to name when `Rscript` is missing. Putting that in a table means adding a sixth backend is a table entry rather than a new branch in three functions.

**Files:**
- Create: `src/microsuite/batch/backends.py`
- Test: `tests/test_batch_backends.py`

**Interfaces:**
- Consumes: `microsuite.methods._dispatch.{require_backend, reject_options}`.
- Produces:
  - `@dataclass(frozen=True) class BatchBackend` with fields `name: str`, `script: str`, `package: str`, `install_hint: str`, `value_type: str`, `supports_covariates: bool`, `requires_target: bool`
  - `BATCH_BACKENDS: dict[str, BatchBackend]`
  - `SUPPORTED_BACKENDS: tuple[str, ...]`
  - `resolve_backend(backend: str, *, covariates: list[str] | None, target: str | None) -> BatchBackend`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_batch_backends.py`:

```python
from __future__ import annotations

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.batch.backends import BATCH_BACKENDS, SUPPORTED_BACKENDS, resolve_backend


def test_every_backend_declares_a_known_scale() -> None:
    from microsuite.batch.value_type import VALUE_TYPES

    for backend in BATCH_BACKENDS.values():
        assert backend.value_type in VALUE_TYPES


def test_supported_backends_matches_the_table() -> None:
    assert set(SUPPORTED_BACKENDS) == set(BATCH_BACKENDS)
    assert "mmuphin" in SUPPORTED_BACKENDS


@pytest.mark.parametrize(
    ("backend", "value_type"),
    [
        ("mmuphin", "relative"),
        ("combat-seq", "counts"),
        ("conqur", "counts"),
        ("plsda-batch", "clr"),
        ("metadict", "relative"),
    ],
)
def test_declared_scales(backend: str, value_type: str) -> None:
    assert BATCH_BACKENDS[backend].value_type == value_type


def test_unknown_backend_lists_the_alternatives() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="mmuphin"):
        resolve_backend("combat", covariates=None, target=None)


def test_covariates_rejected_where_unsupported() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="covariates"):
        resolve_backend("plsda-batch", covariates=["sex"], target="disease")


def test_target_rejected_where_unsupported() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="target"):
        resolve_backend("mmuphin", covariates=None, target="disease")


def test_supervised_backend_requires_a_target() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="supervised"):
        resolve_backend("plsda-batch", covariates=None, target=None)


def test_supported_combination_resolves() -> None:
    backend = resolve_backend("mmuphin", covariates=["sex"], target=None)
    assert backend.name == "mmuphin"
    assert backend.script == "mmuphin"


def test_script_name_differs_from_backend_name_for_combat_seq() -> None:
    assert BATCH_BACKENDS["combat-seq"].script == "combat_seq"
    assert BATCH_BACKENDS["plsda-batch"].script == "plsda_batch"
```

- [ ] **Step 2: Run to verify they fail**

Run: `uv run pytest tests/test_batch_backends.py -v`
Expected: FAIL — `ModuleNotFoundError: microsuite.batch.backends`.

- [ ] **Step 3: Write the table**

`src/microsuite/batch/backends.py`:

```python
"""Capability table for the batch-correction backends.

Adding a backend is a row here, an R script beside it, and a container. The
dispatch reads these records rather than branching on backend names, so a new
backend cannot be half-wired into one function and forgotten in another.
"""

from __future__ import annotations

from dataclasses import dataclass

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods._dispatch import reject_options, require_backend


@dataclass(frozen=True)
class BatchBackend:
    name: str
    script: str
    package: str
    install_hint: str
    value_type: str
    supports_covariates: bool
    requires_target: bool


BATCH_BACKENDS: dict[str, BatchBackend] = {
    "mmuphin": BatchBackend(
        name="mmuphin",
        script="mmuphin",
        package="MMUPHin",
        install_hint="BiocManager::install('MMUPHin')",
        value_type="relative",
        supports_covariates=True,
        requires_target=False,
    ),
    "combat-seq": BatchBackend(
        name="combat-seq",
        script="combat_seq",
        package="sva",
        install_hint="BiocManager::install('sva')",
        value_type="counts",
        supports_covariates=True,
        requires_target=False,
    ),
    "conqur": BatchBackend(
        name="conqur",
        script="conqur",
        package="ConQuR",
        install_hint="remotes::install_github('wdl2459/ConQuR')",
        value_type="counts",
        supports_covariates=True,
        requires_target=False,
    ),
    "plsda-batch": BatchBackend(
        name="plsda-batch",
        script="plsda_batch",
        package="PLSDAbatch",
        install_hint="remotes::install_github('EvaYiwenWang/PLSDAbatch')",
        value_type="clr",
        supports_covariates=False,
        requires_target=True,
    ),
    "metadict": BatchBackend(
        name="metadict",
        script="metadict",
        package="MetaDICT",
        install_hint="remotes::install_github('wangyf1996/MetaDICT')",
        value_type="relative",
        supports_covariates=True,
        requires_target=False,
    ),
}

SUPPORTED_BACKENDS: tuple[str, ...] = tuple(BATCH_BACKENDS)


def resolve_backend(
    backend: str, *, covariates: list[str] | None, target: str | None
) -> BatchBackend:
    """Validate the backend/option combination and return its capability record."""
    name = require_backend(backend, SUPPORTED_BACKENDS, "batch correction")
    record = BATCH_BACKENDS[name]

    unsupported: dict[str, object | None] = {}
    if not record.supports_covariates:
        unsupported["--covariates"] = covariates
    if not record.requires_target:
        unsupported["--target-col"] = target
    if unsupported:
        reject_options(name, unsupported)

    if record.requires_target and not target:
        raise MicrobiomeSuiteError(
            f"--backend {name} is supervised: it fits using the outcome labels, so "
            f"--target-col is required. Note that correcting with the same outcome you "
            f"later test inflates significance; see docs/batch_correction.md."
        )
    return record
```

- [ ] **Step 4: Run to verify they pass**

Run: `uv run pytest tests/test_batch_backends.py -v`
Expected: all PASS.

- [ ] **Step 5: Verify the install hints name real packages**

Run: `uv run python -c "from microsuite.batch.backends import BATCH_BACKENDS; [print(b.name, b.install_hint) for b in BATCH_BACKENDS.values()]"`

Confirm each GitHub path resolves in a browser before Task 7 pins its SHA. If `wangyf1996/MetaDICT` is not the correct repository, correct the hint here and record the real one — do not leave a guess in the error message a user will read.

- [ ] **Step 6: Lint and commit**

```bash
uv run ruff format src/microsuite/batch/backends.py tests/test_batch_backends.py
uv run ruff check src tests
git add src/microsuite/batch/backends.py tests/test_batch_backends.py
git commit -m "feat(batch): add the backend capability table

Covariate support, target requirement, emitted scale, and install hint per
backend, with option rejection driven from the record rather than from
per-backend branches."
```

---

### Task 5: `batch_correct` and the `mmuphin` backend

The first end-to-end path. The riskiest part is not the R call — it is reading the corrected table back. Backends may drop features, reorder rows, or mangle names through `read.delim`, and a misaligned rebuild produces a full table with the wrong numbers against the right labels.

**Files:**
- Create: `src/microsuite/batch/correct.py`, `src/microsuite/batch/r/__init__.py`, `src/microsuite/batch/r/mmuphin.R`, `src/microsuite/methods/batch_correct.py`
- Test: `tests/test_batch_correct.py`

**Interfaces:**
- Consumes: `batch.backends.resolve_backend`, `batch.value_type.record_batch_correction`, `runtime.r_backend.invoke_r_script`, `runtime.container.resolve_batch_image`, `diversity._matrix.dense_counts`, `io.h5ad.{read_h5ad, write_h5ad}`, `_paths.{ensure_input, prepare_output}`.
- Produces:
  - `batch.correct.run_batch_correction(adata: ad.AnnData, *, backend: str, batch: str, covariates: list[str] | None = None, target: str | None = None, extra_params: dict | None = None, run_dir: Path | None = None, timeout: float | None = None, runtime: str = "local", image: str | None = None, engine: str = "docker") -> ad.AnnData`
  - `methods.batch_correct.batch_correct(*, backend: str, table: Path, output: Path, batch: str, covariates: list[str] | None = None, target: str | None = None, force: bool = False, run_dir: Path | None = None, timeout: float | None = None, runtime: str = "local", image: str | None = None, engine: str = "docker") -> None`
  - `methods.batch_correct.SUPPORTED_BACKENDS`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_batch_correct.py`:

```python
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.batch import correct as correct_module
from microsuite.batch.correct import run_batch_correction


def _adata() -> ad.AnnData:
    adata = ad.AnnData(np.array([[10.0, 30.0, 5.0], [20.0, 20.0, 5.0]]))
    adata.obs_names = ["s1", "s2"]
    adata.var_names = ["f1", "f2", "f3"]
    adata.obs = pd.DataFrame(
        {"run_id": ["A", "B"], "sex": ["m", "f"]}, index=adata.obs_names
    )
    return adata


def _fake_backend(write: dict[str, list[float]] | None = None, capture: dict | None = None):
    """Stand in for the R script: record the call, write a corrected table."""

    def _invoke(**kwargs) -> None:
        positional = kwargs["positional"]
        params = json.loads(Path(positional[2]).read_text(encoding="utf-8"))
        if capture is not None:
            capture.update(params=params, kwargs=kwargs)
        payload = write or {"f1": [11.0, 21.0], "f2": [31.0, 21.0], "f3": [6.0, 6.0]}
        frame = pd.DataFrame(payload, index=["s1", "s2"]).T
        frame.index.name = "feature_id"
        frame.to_csv(positional[3], sep="\t")

    return _invoke


def test_params_json_carries_batch_and_covariates(monkeypatch) -> None:
    capture: dict = {}
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend(capture=capture))
    run_batch_correction(_adata(), backend="mmuphin", batch="run_id", covariates=["sex"])
    assert capture["params"]["batch"] == "run_id"
    assert capture["params"]["covariates"] == ["sex"]
    assert capture["kwargs"]["script_name"] == "mmuphin"
    assert capture["kwargs"]["backend"] == "mmuphin"


def test_corrected_values_land_on_the_right_labels(monkeypatch) -> None:
    # The R script returns features as rows in an arbitrary order. A rebuild that
    # trusts position rather than labels puts f3's values under f1.
    payload = {"f3": [6.0, 6.0], "f1": [11.0, 21.0], "f2": [31.0, 21.0]}
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend(write=payload))
    result = run_batch_correction(_adata(), backend="mmuphin", batch="run_id")
    assert list(result.var_names) == ["f1", "f2", "f3"]
    np.testing.assert_allclose(result.X[0], [11.0, 31.0, 6.0])
    np.testing.assert_allclose(result.X[1], [21.0, 21.0, 6.0])


def test_dropped_features_subset_var_rather_than_silently_realigning(monkeypatch) -> None:
    payload = {"f1": [11.0, 21.0], "f2": [31.0, 21.0]}
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend(write=payload))
    result = run_batch_correction(_adata(), backend="mmuphin", batch="run_id")
    assert list(result.var_names) == ["f1", "f2"]
    assert result.shape == (2, 2)


def test_unknown_feature_in_output_raises(monkeypatch) -> None:
    payload = {"f1": [11.0, 21.0], "f9": [1.0, 1.0]}
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend(write=payload))
    with pytest.raises(MicrobiomeSuiteError, match="f9"):
        run_batch_correction(_adata(), backend="mmuphin", batch="run_id")


def test_result_records_its_scale_and_provenance(monkeypatch) -> None:
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend())
    result = run_batch_correction(
        _adata(), backend="mmuphin", batch="run_id", covariates=["sex"]
    )
    assert result.uns["microsuite"]["value_type"] == "relative"
    assert result.uns["microsuite"]["batch_correct"]["backend"] == "mmuphin"
    assert result.uns["microsuite"]["batch_correct"]["covariates"] == ["sex"]


def test_missing_batch_column_lists_available_columns() -> None:
    with pytest.raises(MicrobiomeSuiteError) as excinfo:
        run_batch_correction(_adata(), backend="mmuphin", batch="plate")
    assert "plate" in str(excinfo.value)
    assert "run_id" in str(excinfo.value)


def test_missing_covariate_column_raises() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="age"):
        run_batch_correction(
            _adata(), backend="mmuphin", batch="run_id", covariates=["age"]
        )


def test_single_batch_level_raises() -> None:
    adata = _adata()
    adata.obs["run_id"] = ["A", "A"]
    with pytest.raises(MicrobiomeSuiteError, match="one batch"):
        run_batch_correction(adata, backend="mmuphin", batch="run_id")


def test_covariate_confounded_with_batch_raises() -> None:
    # 'sex' varies exactly with 'run_id' here, so no model can separate them.
    with pytest.raises(MicrobiomeSuiteError, match="confounded"):
        run_batch_correction(
            _adata(), backend="mmuphin", batch="run_id", covariates=["sex"]
        )
```

The last test and `test_params_json_carries_batch_and_covariates` both use `sex`, which is confounded with `run_id` in the two-sample fixture. Give the params test a four-sample fixture where `sex` crosses `run_id`; keep the confounding fixture at two samples. Write `_adata()` with four samples (`run_id = [A, A, B, B]`, `sex = [m, f, m, f]`) and add a separate `_confounded()` helper for the last test.

- [ ] **Step 2: Run to verify they fail**

Run: `uv run pytest tests/test_batch_correct.py -v`
Expected: FAIL — `ModuleNotFoundError: microsuite.batch.correct`.

- [ ] **Step 3: Write `batch/correct.py`**

```python
"""Marshal an AnnData through an R batch-correction backend and back.

The read-back is the delicate half. Backends may drop features and return rows
in their own order, so the corrected matrix is realigned by feature label, never
by position: a positional rebuild yields a complete table with the right labels
on the wrong numbers.
"""

from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, cast

import anndata as ad
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.batch.backends import BatchBackend, resolve_backend
from microsuite.batch.value_type import record_batch_correction
from microsuite.diversity._matrix import dense_counts
from microsuite.runtime.container import resolve_batch_image
from microsuite.runtime.r_backend import invoke_r_script
from microsuite.runtime.runner import CommandLog


def run_batch_correction(
    adata: ad.AnnData,
    *,
    backend: str,
    batch: str,
    covariates: list[str] | None = None,
    target: str | None = None,
    extra_params: dict[str, Any] | None = None,
    run_dir: Path | None = None,
    timeout: float | None = None,
    runtime: str = "local",
    image: str | None = None,
    engine: str = "docker",
) -> ad.AnnData:
    record = resolve_backend(backend, covariates=covariates, target=target)
    covariate_list = list(covariates or [])
    _validate_design(adata, batch=batch, covariates=covariate_list, target=target)

    with TemporaryDirectory() as temp_dir:
        temp = Path(temp_dir)
        counts_path = temp / "counts.tsv"
        metadata_path = temp / "metadata.tsv"
        params_path = temp / "params.json"
        corrected_path = temp / "corrected.tsv"

        pd.DataFrame(
            dense_counts(adata).T, index=adata.var_names, columns=adata.obs_names
        ).to_csv(counts_path, sep="\t")
        pd.DataFrame(adata.obs).to_csv(metadata_path, sep="\t")

        params: dict[str, Any] = {
            "batch": batch,
            "covariates": covariate_list,
            "target": target,
        }
        params.update(extra_params or {})
        params_path.write_text(json.dumps(params, indent=2, sort_keys=True), encoding="utf-8")

        invoke_r_script(
            backend=record.name,
            script_package="microsuite.batch.r",
            script_name=record.script,
            resolve_image=resolve_batch_image,
            positional=[counts_path, metadata_path, params_path, corrected_path],
            runtime=runtime,
            image=image,
            engine=engine,
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(
                task="batch_correct",
                backend=record.name,
                inputs={"batch": batch, "covariates": ",".join(covariate_list)},
                outputs={"output": str(corrected_path)},
            ),
            local_missing_message=(
                f"batch correct --backend {record.name} requires external Rscript and the "
                f"R package '{record.package}'. Install R, then {record.install_hint}, or "
                f"use --runtime docker with the r-batch-{record.name} image."
            ),
        )

        corrected = _read_corrected(corrected_path)

    return _rebuild(adata, corrected, record=record, batch=batch,
                    covariates=covariate_list, target=target)


def _validate_design(
    adata: ad.AnnData, *, batch: str, covariates: list[str], target: str | None
) -> None:
    available = list(adata.obs.columns)
    for label, column in [("--batch-col", batch), ("--target-col", target)]:
        if column is not None and column not in available:
            raise MicrobiomeSuiteError(
                f"{label} '{column}' not found in sample metadata. Available: "
                f"{', '.join(map(str, available))}"
            )
    missing = [name for name in covariates if name not in available]
    if missing:
        raise MicrobiomeSuiteError(
            f"--covariates not found in sample metadata: {', '.join(missing)}. "
            f"Available: {', '.join(map(str, available))}"
        )

    batch_values = adata.obs[batch].astype(str)
    if batch_values.nunique() < 2:
        raise MicrobiomeSuiteError(
            f"'{batch}' has one batch level ({batch_values.iloc[0]}); there is nothing "
            f"to correct."
        )
    for name in covariates:
        values = adata.obs[name].astype(str)
        crosstab = pd.crosstab(batch_values, values)
        # Perfectly confounded: each batch sees exactly one covariate level.
        if (crosstab > 0).sum(axis=1).max() == 1:
            raise MicrobiomeSuiteError(
                f"Covariate '{name}' is perfectly confounded with batch '{batch}': every "
                f"batch contains a single '{name}' level, so no model can separate their "
                f"effects. Drop the covariate, or correct a different grouping."
            )


def _read_corrected(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        raise MicrobiomeSuiteError(
            "The batch-correction backend produced no output table."
        )
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    return frame


def _rebuild(
    adata: ad.AnnData,
    corrected: pd.DataFrame,
    *,
    record: BatchBackend,
    batch: str,
    covariates: list[str],
    target: str | None,
) -> ad.AnnData:
    unknown = [name for name in corrected.index if name not in set(map(str, adata.var_names))]
    if unknown:
        raise MicrobiomeSuiteError(
            f"The backend returned features absent from the input table: "
            f"{', '.join(unknown[:5])}"
        )
    missing_samples = [
        name for name in map(str, adata.obs_names) if name not in set(corrected.columns)
    ]
    if missing_samples:
        raise MicrobiomeSuiteError(
            f"The backend dropped samples from the corrected table: "
            f"{', '.join(missing_samples[:5])}"
        )

    kept = [name for name in map(str, adata.var_names) if name in set(corrected.index)]
    aligned = corrected.loc[kept, [str(name) for name in adata.obs_names]]

    result = cast(Any, adata[:, kept]).copy()
    result.X = aligned.to_numpy(dtype=float).T
    record_batch_correction(
        result,
        value_type=record.value_type,
        backend=record.name,
        batch=batch,
        covariates=covariates,
        target=target,
    )
    return result
```

- [ ] **Step 4: Write `methods/batch_correct.py`**

```python
from __future__ import annotations

from pathlib import Path

from microsuite._paths import ensure_input, prepare_output
from microsuite.batch.backends import SUPPORTED_BACKENDS
from microsuite.batch.correct import run_batch_correction
from microsuite.io.h5ad import read_h5ad, write_h5ad

__all__ = ["SUPPORTED_BACKENDS", "batch_correct"]


def batch_correct(
    *,
    backend: str,
    table: Path,
    output: Path,
    batch: str,
    covariates: list[str] | None = None,
    target: str | None = None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
    runtime: str = "local",
    image: str | None = None,
    engine: str = "docker",
) -> None:
    adata = read_h5ad(ensure_input(table))
    corrected = run_batch_correction(
        adata,
        backend=backend,
        batch=batch,
        covariates=covariates,
        target=target,
        run_dir=run_dir,
        timeout=timeout,
        runtime=runtime,
        image=image,
        engine=engine,
    )
    write_h5ad(corrected, prepare_output(output, force=force))
```

- [ ] **Step 5: Write `batch/r/mmuphin.R`**

Create `src/microsuite/batch/r/__init__.py` (empty file, so `importlib.resources` can find the scripts — mirror `src/microsuite/diffab/r/__init__.py`), then `src/microsuite/batch/r/mmuphin.R`:

```r
#!/usr/bin/env Rscript
# MMUPHin adjust_batch: ComBat extended to zero-inflated microbiome profiles.
# Usage: mmuphin.R counts.tsv metadata.tsv params.json corrected.tsv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: mmuphin.R counts.tsv metadata.tsv params.json corrected.tsv")
}
suppressPackageStartupMessages({
  library(jsonlite)
  library(MMUPHin)
})

counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
meta <- read.delim(args[2], row.names = 1, check.names = FALSE)
params <- fromJSON(args[3], simplifyVector = TRUE)

# The caller writes features as rows and samples as columns; align metadata to
# the count columns rather than trusting the two files to share an order.
meta <- meta[colnames(counts), , drop = FALSE]
meta[[params$batch]] <- factor(meta[[params$batch]])

covariates <- if (length(params$covariates) > 0) as.character(params$covariates) else NULL

fit <- adjust_batch(
  feature_abd = counts,
  batch = params$batch,
  covariates = covariates,
  data = meta
)
adjusted <- fit$feature_abd_adj

out <- data.frame(feature_id = rownames(adjusted), adjusted, check.names = FALSE)
write.table(out, args[4], sep = "\t", quote = FALSE, row.names = FALSE)
```

- [ ] **Step 6: Make the R scripts ship with the package**

Check `pyproject.toml` for how `src/microsuite/diffab/r/*.R` is included in the wheel. Hatchling includes package data under `src/` by default; if `[tool.hatch.build]` names files explicitly, add `src/microsuite/batch/r/*.R` there. Verify with:

```bash
uv build --wheel && uv run python -c "
import zipfile, glob
names = zipfile.ZipFile(sorted(glob.glob('dist/*.whl'))[-1]).namelist()
assert any(n.endswith('batch/r/mmuphin.R') for n in names), 'mmuphin.R missing from wheel'
print('mmuphin.R packaged')
"
```

Expected: `mmuphin.R packaged`. A script that is not in the wheel makes the local runtime fail only after installation, where no test looks.

- [ ] **Step 7: Register in the method API**

In `src/microsuite/cli/_method_api.py`, add beside the other imports:

```python
from microsuite.methods.batch_correct import SUPPORTED_BACKENDS as BATCH_CORRECT_BACKENDS
from microsuite.methods.batch_correct import batch_correct
```

and add to `METHOD_BACKENDS`:

```python
    "batch_correct": BATCH_CORRECT_BACKENDS,
```

- [ ] **Step 8: Run the tests**

Run: `uv run pytest tests/test_batch_correct.py tests/test_methods.py -v`
Expected: all PASS.

- [ ] **Step 9: Break the alignment and confirm the test catches it**

In `_rebuild`, temporarily replace `aligned = corrected.loc[kept, [...]]` with `aligned = corrected.iloc[: len(kept), : adata.n_obs]` — a positional rebuild. Run `uv run pytest tests/test_batch_correct.py -v` and confirm `test_corrected_values_land_on_the_right_labels` **fails**. Revert.

This is the defect class the whole task exists to prevent, so the test must be shown to detect it.

- [ ] **Step 10: Lint and commit**

```bash
uv run ruff format src/microsuite/batch src/microsuite/methods/batch_correct.py src/microsuite/cli/_method_api.py tests/test_batch_correct.py
uv run ruff check src tests
uv run ty check
git add src/microsuite/batch src/microsuite/methods/batch_correct.py src/microsuite/cli/_method_api.py tests/test_batch_correct.py
git commit -m "feat(batch): add batch_correct with the mmuphin backend

Corrected tables are rebuilt by feature label rather than by position, so a
backend that drops or reorders features cannot put the right labels on the
wrong numbers. Covariates perfectly confounded with batch are refused."
```

---

### Task 6: The `r-batch-mmuphin` container

**Files:**
- Create: `containers/r-batch-mmuphin/Dockerfile`, `containers/r-batch-mmuphin/smoke/{counts.tsv,metadata.tsv,params.json}`
- Modify: `tests/test_container_skeletons.py:15-37`, `.github/workflows/docker.yml:117-120`

**Interfaces:**
- Consumes: `src/microsuite/batch/r/mmuphin.R` from Task 5.
- Produces: image `ghcr.io/qwerty239qwe/microsuite/r-batch-mmuphin:latest`, which `resolve_batch_image("mmuphin", None)` already returns.

- [ ] **Step 1: Add the expectation to the skeleton test**

In `tests/test_container_skeletons.py`, add to the `expected` dict:

```python
        "r-batch-mmuphin": ["Rscript", "MMUPHin"],
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_container_skeletons.py -v`
Expected: FAIL — `assert dockerfile.exists()` for `r-batch-mmuphin`.

- [ ] **Step 3: Generate the smoke dataset**

```bash
mkdir -p containers/r-batch-mmuphin/smoke
uv run python - <<'PY'
from pathlib import Path
import json
import numpy as np
import pandas as pd

rng = np.random.default_rng(0)
samples = [f"s{i}" for i in range(12)]
features = [f"f{i}" for i in range(20)]
batch = ["A"] * 6 + ["B"] * 6
group = ["case", "ctrl"] * 6

counts = rng.poisson(50, size=(len(features), len(samples))).astype(float)
# A real batch shift on half the features, so a no-op backend cannot pass.
counts[:10, 6:] *= 3.0

out = Path("containers/r-batch-mmuphin/smoke")
pd.DataFrame(counts, index=features, columns=samples).to_csv(out / "counts.tsv", sep="\t")
pd.DataFrame({"batch": batch, "group": group}, index=samples).to_csv(
    out / "metadata.tsv", sep="\t"
)
(out / "params.json").write_text(
    json.dumps({"batch": "batch", "covariates": ["group"], "target": None}, indent=2) + "\n"
)
print("smoke dataset written")
PY
```

- [ ] **Step 4: Write the Dockerfile**

`containers/r-batch-mmuphin/Dockerfile`:

```dockerfile
FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="r-batch-mmuphin"
LABEL org.opencontainers.image.description="MMUPHin adjust_batch batch-correction backend for microsuite"

ENV PATH="/opt/conda/bin:${PATH}"

# Expected commands: Rscript
# Expected R package: MMUPHin (+ jsonlite for the params.json contract)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "r-base>=4.3,<4.4" \
        r-jsonlite \
        bioconductor-mmuphin && \
    micromamba clean -a -y

COPY src/microsuite/batch/r/mmuphin.R /opt/microsuite/mmuphin.R
COPY containers/r-batch-mmuphin/smoke/ /opt/microsuite/smoke/

# Build-time smoke: run a real correction, not a package import, and fail the
# build unless it emits a table with every input feature and sample.
RUN set -eux; \
    Rscript /opt/microsuite/mmuphin.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-smoke-out.tsv; \
    test -s /tmp/microsuite-smoke-out.tsv; \
    test "$(wc -l < /tmp/microsuite-smoke-out.tsv)" -eq 21; \
    rm -f /tmp/microsuite-smoke-out.tsv

ENTRYPOINT ["Rscript"]
```

The line count is 20 features plus one header. An empty-but-present output passes `test -s`; it does not pass this.

- [ ] **Step 5: Build the image and confirm the smoke ran**

Run: `docker build -f containers/r-batch-mmuphin/Dockerfile -t r-batch-mmuphin:dev .`
Expected: build succeeds and the smoke `RUN` layer shows the `Rscript` invocation.

If Docker is unavailable on this host, mark this step blocked and say so in the task report — do not silently skip it and do not claim the image builds.

- [ ] **Step 6: Register in CI**

In `.github/workflows/docker.yml`, after the `r-ecology` matrix entry (around line 117):

```yaml
          - image: r-batch-mmuphin
            dockerfile: containers/r-batch-mmuphin/Dockerfile
            context: .
            heavy: true
```

- [ ] **Step 7: Run the tests**

Run: `uv run pytest tests/test_container_skeletons.py tests/test_ci_workflow.py -v`
Expected: all PASS. If `test_ci_workflow.py` asserts a container list, update it too.

- [ ] **Step 8: Commit**

```bash
git add containers/r-batch-mmuphin tests/test_container_skeletons.py .github/workflows/docker.yml
git commit -m "build(containers): add r-batch-mmuphin

Build-time smoke runs a real correction and asserts the output row count, so a
backend that emits an empty or truncated table fails the build."
```

---

### Task 7: `microsuite batch correct`

**Files:**
- Create: `src/microsuite/cli/batch_cmd.py`
- Modify: `src/microsuite/cli/app.py:51-64`
- Test: `tests/test_batch_cli.py`

**Interfaces:**
- Consumes: `methods.batch_correct.batch_correct` from Task 5.
- Produces: CLI `microsuite batch correct TABLE --output OUT --batch-col COL [--covariates COL]... [--target-col COL] [--backend NAME] [--runtime local|docker] [--image REF] [--run-dir DIR] [--timeout SECS] [--force]`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_batch_cli.py`:

```python
from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

from microsuite.cli.app import build_app

runner = CliRunner()


def test_batch_correct_is_registered() -> None:
    result = runner.invoke(build_app(), ["batch", "--help"])
    assert result.exit_code == 0
    assert "correct" in result.stdout


def test_options_reach_the_method(monkeypatch, tmp_path: Path) -> None:
    from microsuite.cli import batch_cmd

    captured: dict = {}
    monkeypatch.setattr(batch_cmd, "batch_correct", lambda **kw: captured.update(kw))
    table = tmp_path / "t.h5ad"
    table.write_bytes(b"")
    result = runner.invoke(
        build_app(),
        [
            "batch", "correct", str(table),
            "--output", str(tmp_path / "out.h5ad"),
            "--batch-col", "run_id",
            "--covariates", "sex",
            "--covariates", "age",
            "--backend", "mmuphin",
            "--runtime", "docker",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert captured["batch"] == "run_id"
    assert captured["covariates"] == ["sex", "age"]
    assert captured["backend"] == "mmuphin"
    assert captured["runtime"] == "docker"


def test_default_backend_is_mmuphin(monkeypatch, tmp_path: Path) -> None:
    from microsuite.cli import batch_cmd

    captured: dict = {}
    monkeypatch.setattr(batch_cmd, "batch_correct", lambda **kw: captured.update(kw))
    table = tmp_path / "t.h5ad"
    table.write_bytes(b"")
    result = runner.invoke(
        build_app(),
        ["batch", "correct", str(table), "--output", str(tmp_path / "o.h5ad"),
         "--batch-col", "run_id"],
    )
    assert result.exit_code == 0, result.stdout
    assert captured["backend"] == "mmuphin"
```

Check `src/microsuite/cli/app.py` for the actual app-factory name before writing these; if it is not `build_app`, use whatever `tests/test_cli.py` already imports.

- [ ] **Step 2: Run to verify they fail**

Run: `uv run pytest tests/test_batch_cli.py -v`
Expected: FAIL — no `batch` command.

- [ ] **Step 3: Write the CLI module**

`src/microsuite/cli/batch_cmd.py`:

```python
from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.batch.backends import SUPPORTED_BACKENDS
from microsuite.methods.batch_correct import batch_correct

app = typer.Typer(help="Batch effect correction commands.", no_args_is_help=True)


@app.command("correct")
def correct(
    table: Annotated[Path, typer.Argument(help="Input .h5ad table.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output .h5ad table.")],
    batch: Annotated[
        str, typer.Option("--batch-col", help="obs column holding the batch label.")
    ],
    backend: Annotated[
        str,
        typer.Option("--backend", help=f"One of: {', '.join(SUPPORTED_BACKENDS)}."),
    ] = "mmuphin",
    covariates: Annotated[
        list[str] | None,
        typer.Option("--covariates", help="obs column to preserve (repeatable)."),
    ] = None,
    target: Annotated[
        str | None,
        typer.Option(
            "--target-col",
            help="Outcome column. Required by supervised backends; see "
            "docs/batch_correction.md for the leakage hazard.",
        ),
    ] = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
    run_dir: Annotated[
        Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
    ] = None,
    timeout: Annotated[
        float | None, typer.Option("--timeout", help="Command timeout in seconds.")
    ] = None,
    runtime: Annotated[
        str, typer.Option("--runtime", help="R backend runtime: 'local' Rscript or 'docker'.")
    ] = "local",
    image: Annotated[
        str | None, typer.Option("--image", help="Override the r-batch-<backend> image.")
    ] = None,
) -> None:
    batch_correct(
        backend=backend,
        table=table,
        output=output,
        batch=batch,
        covariates=list(covariates or []),
        target=target,
        force=force,
        run_dir=run_dir,
        timeout=timeout,
        runtime=runtime,
        image=image,
    )
```

- [ ] **Step 4: Register the sub-app**

In `src/microsuite/cli/app.py`, add `batch_cmd` to the imports and, beside the other `add_typer` calls at lines 51-64:

```python
    app.add_typer(batch_cmd.app, name="batch")
```

- [ ] **Step 5: Run the tests**

Run: `uv run pytest tests/test_batch_cli.py tests/test_cli.py -v`
Expected: all PASS.

- [ ] **Step 6: Try it by hand**

Run: `uv run microsuite batch correct --help`
Confirm every option appears with readable help text and that `--target-col` mentions the leakage hazard.

- [ ] **Step 7: Lint and commit**

```bash
uv run ruff format src/microsuite/cli/batch_cmd.py src/microsuite/cli/app.py tests/test_batch_cli.py
uv run ruff check src tests
git add src/microsuite/cli/batch_cmd.py src/microsuite/cli/app.py tests/test_batch_cli.py
git commit -m "feat(cli): add 'microsuite batch correct'"
```

---

### Task 8: The real-execution smoke test

Every test so far mocks the subprocess, so all of them together prove only that the commands we intended get constructed. The mothur work shipped fourteen defects past exactly that kind of green suite.

This test runs the real backend and asserts on the biology: the batch effect must shrink **and** the biological signal must survive. Asserting only the first passes a backend that flattens the table to zeros.

**Files:**
- Create: `tests/integration/test_batch_correct_smoke.py`
- Test: itself

**Interfaces:**
- Consumes: `batch.correct.run_batch_correction`; `diversity.ecology.beta_significance` with `backend="vegan"`, `method="adonis2"`; `diversity.beta` for the distance matrix.
- Produces: `_two_batch_dataset(seed: int = 0) -> ad.AnnData` and `_batch_and_group_r2(adata) -> tuple[float, float]`, reused by the backend tasks that follow.

- [ ] **Step 1: Write the test**

Create `tests/integration/test_batch_correct_smoke.py`:

```python
"""Real-execution smoke test for the batch-correction backends.

Gated behind ``MICROSUITE_RUN_BATCH_SMOKE=1`` and a container engine, because
it runs the real R backends and the vegan image.

**Why this test exists.** Every other batch test mocks the subprocess, so
together they prove only that we build the commands we meant to build. The
mothur work shipped fourteen defects past a suite of exactly that kind, each
producing a complete, well-formed, wrong result.

The assertion is deliberately a *pair*. A backend that shrinks the batch effect
by flattening every difference in the table scores perfectly on the batch term
alone. Both terms are checked, always.

The dataset is generated deterministically rather than committed, so there is no
fixture to drift.
"""

from __future__ import annotations

import os
import shutil

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from microsuite.batch.correct import run_batch_correction
from microsuite.diversity.beta import beta_diversity
from microsuite.diversity.ecology import beta_significance

pytestmark = [
    pytest.mark.skipif(
        os.environ.get("MICROSUITE_RUN_BATCH_SMOKE") != "1",
        reason="set MICROSUITE_RUN_BATCH_SMOKE=1 to run the real batch backends",
    ),
    pytest.mark.skipif(shutil.which("docker") is None, reason="docker is not installed"),
]

N_PER_CELL = 8
N_FEATURES = 60


def _two_batch_dataset(seed: int = 0) -> ad.AnnData:
    """Two batches crossed with two groups, each carrying a distinct shift.

    The design is balanced: every batch contains both groups, so batch and group
    are separable. A correction that removes the group effect along with the
    batch effect is a failure, and this layout is what makes that visible.
    """
    rng = np.random.default_rng(seed)
    rows, batches, groups, names = [], [], [], []
    batch_shift = np.ones(N_FEATURES)
    batch_shift[: N_FEATURES // 2] = 4.0  # multiplicative, half the features
    group_shift = np.ones(N_FEATURES)
    group_shift[N_FEATURES // 2 :] = 2.0  # a smaller effect on the other half

    for batch in ("A", "B"):
        for group in ("case", "ctrl"):
            for replicate in range(N_PER_CELL):
                base = rng.gamma(shape=2.0, scale=30.0, size=N_FEATURES)
                if batch == "B":
                    base = base * batch_shift
                if group == "case":
                    base = base * group_shift
                rows.append(rng.poisson(base).astype(float))
                batches.append(batch)
                groups.append(group)
                names.append(f"{batch}{group}{replicate}")

    adata = ad.AnnData(np.vstack(rows))
    adata.obs_names = names
    adata.var_names = [f"f{i}" for i in range(N_FEATURES)]
    adata.obs = pd.DataFrame({"run_id": batches, "group": groups}, index=names)
    return adata


def _batch_and_group_r2(adata: ad.AnnData) -> tuple[float, float]:
    """PERMANOVA variance explained by batch and by group, via vegan adonis2."""
    distances = beta_diversity(adata, metric="braycurtis")
    result = beta_significance(
        distances,
        pd.DataFrame(adata.obs),
        method="adonis2",
        formula="run_id + group",
        backend="vegan",
        permutations=199,
        seed=0,
        runtime="docker",
    )
    indexed = result.set_index("term")["r_squared"]
    return float(indexed["run_id"]), float(indexed["group"])


def test_uncorrected_dataset_has_the_effects_the_test_assumes() -> None:
    # If the generator stops producing a batch effect, every downstream
    # assertion below becomes vacuously true. Check the premise first.
    batch_r2, group_r2 = _batch_and_group_r2(_two_batch_dataset())
    assert batch_r2 > 0.10, f"generated batch effect too small: {batch_r2}"
    assert group_r2 > 0.02, f"generated group effect too small: {group_r2}"


def test_mmuphin_shrinks_batch_and_keeps_group() -> None:
    adata = _two_batch_dataset()
    before_batch, before_group = _batch_and_group_r2(adata)

    corrected = run_batch_correction(
        adata, backend="mmuphin", batch="run_id", covariates=["group"], runtime="docker"
    )
    after_batch, after_group = _batch_and_group_r2(corrected)

    assert after_batch < before_batch * 0.5, (
        f"batch R2 did not shrink: {before_batch:.3f} -> {after_batch:.3f}"
    )
    assert after_group > before_group * 0.5, (
        f"biological signal was flattened along with the batch effect: "
        f"{before_group:.3f} -> {after_group:.3f}"
    )
    assert corrected.uns["microsuite"]["value_type"] == "relative"
```

Check `microsuite.diversity.beta` for the real function name and metric argument before writing — `tests/test_beta.py` shows the current call shape. Use whatever it uses.

- [ ] **Step 2: Run the premise test**

Run: `MICROSUITE_RUN_BATCH_SMOKE=1 uv run pytest tests/integration/test_batch_correct_smoke.py::test_uncorrected_dataset_has_the_effects_the_test_assumes -v`
Expected: PASS. If the generated batch effect is under 0.10, raise `batch_shift` until it is not — an assertion against an effect that does not exist proves nothing.

- [ ] **Step 3: Run the mmuphin smoke**

Run: `MICROSUITE_RUN_BATCH_SMOKE=1 uv run pytest tests/integration/test_batch_correct_smoke.py -v`
Expected: PASS.

If it fails, the backend is wrong, not the test. Do not loosen the thresholds to make it pass. If the correction genuinely cannot reach a 50% reduction on this dataset, record the measured numbers in the task report and raise it rather than editing the assertion quietly.

- [ ] **Step 4: Confirm the test detects a no-op**

Temporarily change `_rebuild` in `batch/correct.py` to return the input unchanged (`result.X = dense_counts(adata)`), rerun the smoke, and confirm `test_mmuphin_shrinks_batch_and_keeps_group` **fails** on the batch assertion. Revert.

- [ ] **Step 5: Confirm the test detects a flattening**

Temporarily change `_rebuild` to zero the matrix (`result.X = np.zeros_like(...)`), rerun, and confirm the test fails. Whether it fails on the batch or the group assertion, it must fail. Revert.

Both breaks matter: they are the two opposite ways a correction can be wrong, and a test that catches only one of them is half a test.

- [ ] **Step 6: Commit**

```bash
git add tests/integration/test_batch_correct_smoke.py
git commit -m "test(batch): add a real-execution smoke test for mmuphin

Asserts the pair: batch R2 shrinks and group R2 survives. Verified to fail
against both a no-op correction and a flattening one."
```

---

### Task 9: The `combat-seq` backend

Negative-binomial correction returning integer counts. This is the backend that matters for ANCOM-BC and ALDEx2 users, because it is the only one in phase 1 whose output they can consume without complaint.

**Files:**
- Create: `src/microsuite/batch/r/combat_seq.R`, `containers/r-batch-combatseq/{Dockerfile,smoke/}`
- Modify: `tests/test_container_skeletons.py`, `.github/workflows/docker.yml`, `tests/integration/test_batch_correct_smoke.py`
- Test: `tests/test_batch_correct.py` (add a case)

**Interfaces:**
- Consumes: everything from Tasks 4, 5, 6, 8. The `BATCH_BACKENDS["combat-seq"]` record already exists from Task 4.
- Produces: no new Python symbols.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_batch_correct.py`:

```python
def test_combat_seq_declares_counts(monkeypatch) -> None:
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend())
    result = run_batch_correction(_adata(), backend="combat-seq", batch="run_id")
    assert result.uns["microsuite"]["value_type"] == "counts"


def test_combat_seq_uses_its_own_script_name(monkeypatch) -> None:
    capture: dict = {}
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend(capture=capture))
    run_batch_correction(_adata(), backend="combat-seq", batch="run_id")
    assert capture["kwargs"]["script_name"] == "combat_seq"
    assert capture["kwargs"]["backend"] == "combat-seq"
```

Append to `tests/test_container_skeletons.py`'s `expected` dict:

```python
        "r-batch-combatseq": ["Rscript", "sva", "ComBat_seq"],
```

- [ ] **Step 2: Run to verify they fail**

Run: `uv run pytest tests/test_batch_correct.py tests/test_container_skeletons.py -v`
Expected: the skeleton test FAILS on the missing Dockerfile. The two Python tests may already pass from the Task 4 table — that is fine; they are regression cover for the script-name divergence.

- [ ] **Step 3: Write `combat_seq.R`**

```r
#!/usr/bin/env Rscript
# sva::ComBat_seq: negative-binomial batch adjustment returning integer counts.
# Usage: combat_seq.R counts.tsv metadata.tsv params.json corrected.tsv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: combat_seq.R counts.tsv metadata.tsv params.json corrected.tsv")
}
suppressPackageStartupMessages({
  library(jsonlite)
  library(sva)
})

counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
meta <- read.delim(args[2], row.names = 1, check.names = FALSE)
params <- fromJSON(args[3], simplifyVector = TRUE)

meta <- meta[colnames(counts), , drop = FALSE]

# ComBat_seq models counts, so a non-integer matrix means the caller fed it
# something already normalized. Fail rather than round silently.
if (any(abs(counts - round(counts)) > 1e-8)) {
  stop("ComBat_seq requires integer counts; this table holds non-integer values.")
}
storage.mode(counts) <- "integer"

batch <- factor(meta[[params$batch]])

covar_mod <- NULL
if (length(params$covariates) > 0) {
  covariates <- as.character(params$covariates)
  design <- meta[, covariates, drop = FALSE]
  for (name in covariates) {
    if (is.character(design[[name]])) design[[name]] <- factor(design[[name]])
  }
  covar_mod <- model.matrix(~., data = design)
}

adjusted <- ComBat_seq(counts = counts, batch = batch, group = NULL, covar_mod = covar_mod)

out <- data.frame(feature_id = rownames(adjusted), adjusted, check.names = FALSE)
write.table(out, args[4], sep = "\t", quote = FALSE, row.names = FALSE)
```

- [ ] **Step 4: Generate the smoke dataset and write the Dockerfile**

Reuse the generator from Task 6 Step 3, writing to `containers/r-batch-combatseq/smoke/`, with one change: `ComBat_seq` needs integer input, so drop the `* 3.0` float multiply in favour of `counts[:10, 6:] = counts[:10, 6:] * 3` on an integer array, and write with `.astype(int)`.

`containers/r-batch-combatseq/Dockerfile`:

```dockerfile
FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="r-batch-combatseq"
LABEL org.opencontainers.image.description="sva ComBat_seq batch-correction backend for microsuite"

ENV PATH="/opt/conda/bin:${PATH}"

# Expected commands: Rscript
# Expected R package: sva (ComBat_seq) (+ jsonlite for the params.json contract)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "r-base>=4.3,<4.4" \
        r-jsonlite \
        bioconductor-sva && \
    micromamba clean -a -y

COPY src/microsuite/batch/r/combat_seq.R /opt/microsuite/combat_seq.R
COPY containers/r-batch-combatseq/smoke/ /opt/microsuite/smoke/

# Build-time smoke: a real correction, with an integer-output check. ComBat_seq
# returning floats would silently break every count-requiring downstream method.
RUN set -eux; \
    Rscript /opt/microsuite/combat_seq.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-smoke-out.tsv; \
    test -s /tmp/microsuite-smoke-out.tsv; \
    test "$(wc -l < /tmp/microsuite-smoke-out.tsv)" -eq 21; \
    Rscript -e "x <- read.delim('/tmp/microsuite-smoke-out.tsv', row.names=1); \
                stopifnot(all(abs(as.matrix(x) - round(as.matrix(x))) < 1e-8))"; \
    rm -f /tmp/microsuite-smoke-out.tsv

ENTRYPOINT ["Rscript"]
```

- [ ] **Step 5: Build the image**

Run: `docker build -f containers/r-batch-combatseq/Dockerfile -t r-batch-combatseq:dev .`
Expected: success, including the integer check. If Docker is unavailable, report the step as blocked rather than skipping it silently.

- [ ] **Step 6: Register in CI**

Add to `.github/workflows/docker.yml` after the `r-batch-mmuphin` entry:

```yaml
          - image: r-batch-combatseq
            dockerfile: containers/r-batch-combatseq/Dockerfile
            context: .
            heavy: true
```

- [ ] **Step 7: Add the smoke case**

In `tests/integration/test_batch_correct_smoke.py`:

```python
def test_combat_seq_shrinks_batch_keeps_group_and_returns_counts() -> None:
    adata = _two_batch_dataset()
    before_batch, before_group = _batch_and_group_r2(adata)

    corrected = run_batch_correction(
        adata, backend="combat-seq", batch="run_id", covariates=["group"], runtime="docker"
    )
    after_batch, after_group = _batch_and_group_r2(corrected)

    assert after_batch < before_batch * 0.5, (
        f"batch R2 did not shrink: {before_batch:.3f} -> {after_batch:.3f}"
    )
    assert after_group > before_group * 0.5, (
        f"biological signal was flattened: {before_group:.3f} -> {after_group:.3f}"
    )
    assert corrected.uns["microsuite"]["value_type"] == "counts"
    # The whole point of this backend is that ANCOM-BC can consume its output.
    values = np.asarray(corrected.X)
    assert np.allclose(values, np.round(values)), "ComBat_seq returned non-integer counts"
```

- [ ] **Step 8: Run everything**

Run: `uv run pytest tests/test_batch_correct.py tests/test_container_skeletons.py -v`
Then: `MICROSUITE_RUN_BATCH_SMOKE=1 uv run pytest tests/integration/test_batch_correct_smoke.py -v`
Expected: all PASS.

- [ ] **Step 9: Commit**

```bash
uv run ruff format tests/test_batch_correct.py
uv run ruff check src tests
git add src/microsuite/batch/r/combat_seq.R containers/r-batch-combatseq tests .github/workflows/docker.yml
git commit -m "feat(batch): add the combat-seq backend

Returns integer counts, so ANCOM-BC and ALDEx2 can consume corrected tables.
The container smoke and the integration smoke both assert integrality, because
float output would pass every other check and break downstream silently."
```

---

### Task 10: The `conqur` backend

The first GitHub-sourced package. ConQuR corrects by conditional quantile regression and returns a corrected count table.

**Before writing the script:** ConQuR's exported signature must be read from the pinned commit, not assumed. The plan's script assumes `ConQuR(tax_tab, batchid, covariates, batch_ref)` where `tax_tab` is **samples × taxa** — the transpose of microsuite's on-disk orientation — and that it returns a samples × taxa data frame. Verify with `Rscript -e "args(ConQuR::ConQuR)"` inside the built image and adjust the script if it differs. Record what you found in the task report either way.

**Files:**
- Create: `src/microsuite/batch/r/conqur.R`, `containers/r-batch-conqur/{Dockerfile,smoke/}`
- Modify: `tests/test_container_skeletons.py`, `.github/workflows/docker.yml`, `tests/integration/test_batch_correct_smoke.py`

**Interfaces:**
- Consumes: Tasks 4, 5, 8.
- Produces: no new Python symbols. `batch_ref` travels through `run_batch_correction(extra_params={"batch_ref": ...})`.

- [ ] **Step 1: Add the container expectation**

In `tests/test_container_skeletons.py`:

```python
        "r-batch-conqur": ["Rscript", "ConQuR"],
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_container_skeletons.py -v`
Expected: FAIL on the missing Dockerfile.

- [ ] **Step 3: Pin the commit**

```bash
git ls-remote https://github.com/wdl2459/ConQuR HEAD
```

Record the SHA. It goes in the Dockerfile literally — a branch name would let the image change under us with no change on our side.

- [ ] **Step 4: Write `conqur.R`**

```r
#!/usr/bin/env Rscript
# ConQuR: conditional quantile regression batch removal, returning counts.
# Usage: conqur.R counts.tsv metadata.tsv params.json corrected.tsv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: conqur.R counts.tsv metadata.tsv params.json corrected.tsv")
}
suppressPackageStartupMessages({
  library(jsonlite)
  library(ConQuR)
})

counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
meta <- read.delim(args[2], row.names = 1, check.names = FALSE)
params <- fromJSON(args[3], simplifyVector = TRUE)

meta <- meta[colnames(counts), , drop = FALSE]

# ConQuR takes samples as rows; microsuite writes features as rows.
tax_tab <- as.data.frame(t(counts))
batchid <- factor(meta[[params$batch]])

# ConQuR requires a reference batch. Default to the first level so runs are
# reproducible rather than dependent on factor ordering elsewhere.
batch_ref <- if (!is.null(params$batch_ref)) as.character(params$batch_ref) else levels(batchid)[1]

if (length(params$covariates) > 0) {
  covariates <- meta[, as.character(params$covariates), drop = FALSE]
  for (name in names(covariates)) {
    if (is.character(covariates[[name]])) covariates[[name]] <- factor(covariates[[name]])
  }
} else {
  # ConQuR requires a covariate frame; an intercept-only frame is the no-covariate case.
  covariates <- data.frame(intercept_only = factor(rep("a", nrow(tax_tab))))
}

adjusted <- ConQuR(
  tax_tab = tax_tab,
  batchid = batchid,
  covariates = covariates,
  batch_ref = batch_ref
)

# Back to features-as-rows for the caller.
out_matrix <- t(as.matrix(adjusted))
out <- data.frame(feature_id = rownames(out_matrix), out_matrix, check.names = FALSE)
write.table(out, args[4], sep = "\t", quote = FALSE, row.names = FALSE)
```

- [ ] **Step 5: Write the Dockerfile**

Generate `containers/r-batch-conqur/smoke/` with the Task 6 Step 3 generator (integer counts, as for combat-seq), then:

```dockerfile
FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="r-batch-conqur"
LABEL org.opencontainers.image.description="ConQuR quantile-regression batch-correction backend for microsuite"

ENV PATH="/opt/conda/bin:${PATH}"

# Expected commands: Rscript
# Expected R package: ConQuR (GitHub, pinned by commit) (+ jsonlite)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "r-base>=4.3,<4.4" \
        r-jsonlite \
        r-remotes \
        r-quantreg \
        r-glmnet \
        r-doparallel \
        r-dplyr \
        r-ade4 \
        r-compositions \
        r-rocr \
        r-fastdummies \
        r-randomforest \
        r-gunifrac && \
    micromamba clean -a -y

# Pinned to a commit, never a branch: ConQuR has no release tags, so a branch
# reference would let this image change with no change on our side.
ARG CONQUR_SHA=<sha-from-step-3>
RUN Rscript -e "remotes::install_github('wdl2459/ConQuR', ref = '${CONQUR_SHA}', upgrade = 'never')" && \
    Rscript -e "stopifnot(requireNamespace('ConQuR', quietly = TRUE))"

COPY src/microsuite/batch/r/conqur.R /opt/microsuite/conqur.R
COPY containers/r-batch-conqur/smoke/ /opt/microsuite/smoke/

RUN set -eux; \
    Rscript /opt/microsuite/conqur.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-smoke-out.tsv; \
    test -s /tmp/microsuite-smoke-out.tsv; \
    test "$(wc -l < /tmp/microsuite-smoke-out.tsv)" -eq 21; \
    rm -f /tmp/microsuite-smoke-out.tsv

ENTRYPOINT ["Rscript"]
```

Replace `<sha-from-step-3>` with the real SHA. If any listed dependency is not on conda-forge or bioconda under that name, install it via `remotes` in the same pinned step and note the substitution in the task report — do not leave a package name in the Dockerfile that does not resolve.

- [ ] **Step 6: Build and verify the signature assumption**

```bash
docker build -f containers/r-batch-conqur/Dockerfile -t r-batch-conqur:dev .
docker run --rm --entrypoint Rscript r-batch-conqur:dev -e "print(args(ConQuR::ConQuR))"
```

Compare against the script's call. If the parameter names or the expected orientation differ, fix `conqur.R` and rebuild.

- [ ] **Step 7: Register in CI and add the smoke case**

Matrix entry as in Task 9 Step 6, with `r-batch-conqur`. Then add to `tests/integration/test_batch_correct_smoke.py` a `test_conqur_shrinks_batch_and_keeps_group` copying the body of `test_combat_seq_shrinks_batch_keeps_group_and_returns_counts`, with `backend="conqur"` and the same integrality assertion, since ConQuR also declares `counts`.

- [ ] **Step 8: Run everything**

Run: `uv run pytest tests/test_container_skeletons.py -v`
Then: `MICROSUITE_RUN_BATCH_SMOKE=1 uv run pytest tests/integration/test_batch_correct_smoke.py -v`
Expected: all PASS.

- [ ] **Step 9: Commit**

```bash
git add src/microsuite/batch/r/conqur.R containers/r-batch-conqur tests .github/workflows/docker.yml
git commit -m "feat(batch): add the conqur backend

Installed from a pinned commit, since ConQuR publishes no release tags. The
script transposes to ConQuR's samples-as-rows orientation and back."
```

---

### Task 11: The `plsda-batch` backend

The first supervised backend and the only one emitting CLR. Both facts are already encoded in the Task 4 record; this task makes them real.

**Before writing the script:** verify `PLSDAbatch::PLSDA_batch`'s signature in the built image. The plan assumes `PLSDA_batch(X, Y.trt, Y.bat, ncomp.trt, ncomp.bat)` with `X` as **samples × taxa in CLR space**, returning a list with `X.nobatch`. Confirm with `Rscript -e "args(PLSDAbatch::PLSDA_batch)"` and adjust if it differs.

**Files:**
- Create: `src/microsuite/batch/r/plsda_batch.R`, `containers/r-batch-plsdabatch/{Dockerfile,smoke/}`
- Modify: `tests/test_container_skeletons.py`, `.github/workflows/docker.yml`, `tests/integration/test_batch_correct_smoke.py`, `tests/test_batch_correct.py`

**Interfaces:**
- Consumes: Tasks 4, 5, 8.
- Produces: no new Python symbols.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_batch_correct.py`:

```python
def test_plsda_batch_requires_a_target() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="supervised"):
        run_batch_correction(_adata(), backend="plsda-batch", batch="run_id")


def test_plsda_batch_rejects_covariates() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="covariates"):
        run_batch_correction(
            _adata(), backend="plsda-batch", batch="run_id",
            covariates=["sex"], target="group",
        )


def test_plsda_batch_declares_clr_and_passes_the_target(monkeypatch) -> None:
    capture: dict = {}
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend(capture=capture))
    result = run_batch_correction(
        _adata(), backend="plsda-batch", batch="run_id", target="group"
    )
    assert capture["params"]["target"] == "group"
    assert result.uns["microsuite"]["value_type"] == "clr"
    assert result.uns["microsuite"]["batch_correct"]["target"] == "group"
```

And in `tests/test_container_skeletons.py`:

```python
        "r-batch-plsdabatch": ["Rscript", "PLSDAbatch"],
```

- [ ] **Step 2: Run to verify they fail**

Run: `uv run pytest tests/test_batch_correct.py tests/test_container_skeletons.py -v`
Expected: the skeleton test FAILS on the missing Dockerfile. The three Python tests should pass off the Task 4 record; if any does not, the capability table is wrong — fix it there, not here.

- [ ] **Step 3: Pin the commit**

```bash
git ls-remote https://github.com/EvaYiwenWang/PLSDAbatch HEAD
```

- [ ] **Step 4: Write `plsda_batch.R`**

```r
#!/usr/bin/env Rscript
# PLSDA-batch: subtracts batch-associated latent components in CLR space.
# Output is CLR log-ratios, NOT counts or relative abundances.
# Usage: plsda_batch.R counts.tsv metadata.tsv params.json corrected.tsv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: plsda_batch.R counts.tsv metadata.tsv params.json corrected.tsv")
}
suppressPackageStartupMessages({
  library(jsonlite)
  library(mixOmics)
  library(PLSDAbatch)
})

counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
meta <- read.delim(args[2], row.names = 1, check.names = FALSE)
params <- fromJSON(args[3], simplifyVector = TRUE)

meta <- meta[colnames(counts), , drop = FALSE]

if (is.null(params$target) || is.na(params$target)) {
  stop("plsda-batch is supervised and requires a target column.")
}

# PLSDA-batch operates on CLR-transformed samples-as-rows data. The offset makes
# the log defined for zero counts; it is the standard mixOmics convention.
X <- t(counts)
X <- logratio.transfo(X = X + 1, logratio = "CLR")
class(X) <- "matrix"

Y.trt <- factor(meta[[params$target]])
Y.bat <- factor(meta[[params$batch]])

ncomp_trt <- if (!is.null(params$ncomp_trt)) as.integer(params$ncomp_trt) else 2L
ncomp_bat <- if (!is.null(params$ncomp_bat)) as.integer(params$ncomp_bat) else 2L

fit <- PLSDA_batch(
  X = X,
  Y.trt = Y.trt,
  Y.bat = Y.bat,
  ncomp.trt = ncomp_trt,
  ncomp.bat = ncomp_bat
)
adjusted <- fit$X.nobatch

out_matrix <- t(as.matrix(adjusted))
out <- data.frame(feature_id = rownames(out_matrix), out_matrix, check.names = FALSE)
write.table(out, args[4], sep = "\t", quote = FALSE, row.names = FALSE)
```

- [ ] **Step 5: Write the Dockerfile**

```dockerfile
FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="r-batch-plsdabatch"
LABEL org.opencontainers.image.description="PLSDA-batch CLR-space batch-correction backend for microsuite"

ENV PATH="/opt/conda/bin:${PATH}"

# Expected commands: Rscript
# Expected R package: PLSDAbatch (GitHub, pinned by commit) (+ mixOmics, jsonlite)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "r-base>=4.3,<4.4" \
        r-jsonlite \
        r-remotes \
        bioconductor-mixomics && \
    micromamba clean -a -y

ARG PLSDABATCH_SHA=<sha-from-step-3>
RUN Rscript -e "remotes::install_github('EvaYiwenWang/PLSDAbatch', ref = '${PLSDABATCH_SHA}', upgrade = 'never')" && \
    Rscript -e "stopifnot(requireNamespace('PLSDAbatch', quietly = TRUE))"

COPY src/microsuite/batch/r/plsda_batch.R /opt/microsuite/plsda_batch.R
COPY containers/r-batch-plsdabatch/smoke/ /opt/microsuite/smoke/

# The smoke params carry a target, because this backend is supervised.
RUN set -eux; \
    Rscript /opt/microsuite/plsda_batch.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-smoke-out.tsv; \
    test -s /tmp/microsuite-smoke-out.tsv; \
    test "$(wc -l < /tmp/microsuite-smoke-out.tsv)" -eq 21; \
    rm -f /tmp/microsuite-smoke-out.tsv

ENTRYPOINT ["Rscript"]
```

Generate `containers/r-batch-plsdabatch/smoke/` with the Task 6 generator, then edit `params.json` to `{"batch": "batch", "covariates": [], "target": "group"}`.

- [ ] **Step 6: Build and verify the signature assumption**

```bash
docker build -f containers/r-batch-plsdabatch/Dockerfile -t r-batch-plsdabatch:dev .
docker run --rm --entrypoint Rscript r-batch-plsdabatch:dev -e "print(args(PLSDAbatch::PLSDA_batch))"
```

Fix `plsda_batch.R` and rebuild if the signature differs from the assumption.

- [ ] **Step 7: Register in CI and add the smoke case**

Matrix entry as before. Then in `tests/integration/test_batch_correct_smoke.py`:

```python
def test_plsda_batch_shrinks_batch_keeps_group_and_returns_clr() -> None:
    adata = _two_batch_dataset()
    before_batch, before_group = _batch_and_group_r2(adata)

    corrected = run_batch_correction(
        adata, backend="plsda-batch", batch="run_id", target="group", runtime="docker"
    )
    # CLR output is not a distance-compatible abundance table for Bray-Curtis, so
    # this backend is assessed on Euclidean distance over the CLR values, which is
    # the Aitchison distance the method itself is defined against.
    after_batch, after_group = _batch_and_group_r2_euclidean(corrected)
    before_batch_e, before_group_e = _batch_and_group_r2_euclidean(adata, clr=True)

    assert after_batch < before_batch_e * 0.5, (
        f"batch R2 did not shrink: {before_batch_e:.3f} -> {after_batch:.3f}"
    )
    assert after_group > before_group_e * 0.5, (
        f"biological signal was flattened: {before_group_e:.3f} -> {after_group:.3f}"
    )
    assert corrected.uns["microsuite"]["value_type"] == "clr"
```

Add the helper beside `_batch_and_group_r2`:

```python
def _batch_and_group_r2_euclidean(adata: ad.AnnData, *, clr: bool = False) -> tuple[float, float]:
    """Same partition, on Euclidean distance — the right geometry for CLR data.

    Bray-Curtis on CLR values is meaningless (the values are signed), so
    comparing a CLR-space correction against a Bray-Curtis baseline would
    compare two different quantities and call the difference an improvement.
    """
    from scipy.spatial.distance import pdist, squareform

    from microsuite.methods.normalize import normalize_native

    source = normalize_native(adata, method="clr") if clr else adata
    values = np.asarray(source.X, dtype=float)
    distances = pd.DataFrame(
        squareform(pdist(values, metric="euclidean")),
        index=list(map(str, adata.obs_names)),
        columns=list(map(str, adata.obs_names)),
    )
    result = beta_significance(
        distances,
        pd.DataFrame(adata.obs),
        method="adonis2",
        formula="run_id + group",
        backend="vegan",
        permutations=199,
        seed=0,
        runtime="docker",
    )
    indexed = result.set_index("term")["r_squared"]
    return float(indexed["run_id"]), float(indexed["group"])
```

- [ ] **Step 8: Run everything**

Run: `uv run pytest tests/test_batch_correct.py tests/test_container_skeletons.py -v`
Then: `MICROSUITE_RUN_BATCH_SMOKE=1 uv run pytest tests/integration/test_batch_correct_smoke.py -v`
Expected: all PASS.

- [ ] **Step 9: Commit**

```bash
uv run ruff format tests/integration/test_batch_correct_smoke.py tests/test_batch_correct.py
uv run ruff check src tests
git add src/microsuite/batch/r/plsda_batch.R containers/r-batch-plsdabatch tests .github/workflows/docker.yml
git commit -m "feat(batch): add the supervised plsda-batch backend

Emits CLR log-ratios, so its smoke test is assessed on Euclidean (Aitchison)
distance rather than Bray-Curtis; comparing signed CLR values under Bray-Curtis
would compare two different quantities."
```

---

### Task 12: The `metadict` backend

Covariate balancing plus shared dictionary learning.

**Before writing the script:** MetaDICT's repository and exported signature must both be confirmed, not assumed. Task 4 Step 5 already asked you to verify the repository path in the install hint. Here, confirm the entry-point function name and its arguments inside the built image before writing the final script. The plan's script assumes a function `MetaDICT(count, meta, batch, covariates)` returning a list whose corrected abundance matrix is `$X` — treat every part of that as unverified.

If the real API differs materially, write the script against the real one and record the difference in the task report. Do not adapt microsuite's four-argument script contract to fit the package.

**Files:**
- Create: `src/microsuite/batch/r/metadict.R`, `containers/r-batch-metadict/{Dockerfile,smoke/}`
- Modify: `tests/test_container_skeletons.py`, `.github/workflows/docker.yml`, `tests/integration/test_batch_correct_smoke.py`

**Interfaces:**
- Consumes: Tasks 4, 5, 8.
- Produces: no new Python symbols.

- [ ] **Step 1: Add the container expectation and run it red**

```python
        "r-batch-metadict": ["Rscript", "MetaDICT"],
```

Run: `uv run pytest tests/test_container_skeletons.py -v`
Expected: FAIL on the missing Dockerfile.

- [ ] **Step 2: Confirm the repository and pin the commit**

```bash
git ls-remote <confirmed-metadict-repo> HEAD
```

Use the repository confirmed in Task 4 Step 5. If it turns out `BATCH_BACKENDS["metadict"].install_hint` points somewhere wrong, correct it in `backends.py` as part of this task.

- [ ] **Step 3: Build a probe image and read the real API**

Write a minimal Dockerfile installing only R, `remotes`, and MetaDICT at the pinned SHA, build it, then:

```bash
docker run --rm --entrypoint Rscript metadict-probe:dev -e "print(ls('package:MetaDICT'))"
docker run --rm --entrypoint Rscript metadict-probe:dev -e "library(MetaDICT); print(args(MetaDICT))"
```

Write down the entry point, its arguments, its expected orientation, and the shape of its return value. Everything in the next step depends on this, and guessing here produces a script that runs and returns something wrong.

- [ ] **Step 4: Write `metadict.R` against the real signature**

Follow the shape of `mmuphin.R`: read four positional arguments, align metadata to the count columns by name, coerce the batch column to a factor, pass covariates when the params list is non-empty, and write the corrected matrix as `feature_id` plus one column per sample. Convert MetaDICT's orientation to features-as-rows before writing, as `conqur.R` does.

The output contract is fixed regardless of what the package returns: a TSV whose first column is `feature_id` and whose remaining columns are the input sample names.

- [ ] **Step 5: Write the Dockerfile**

Follow `containers/r-batch-conqur/Dockerfile` exactly in structure: micromamba base, `r-jsonlite` and `r-remotes` plus whatever MetaDICT's DESCRIPTION requires, a pinned `ARG METADICT_SHA` install, a `requireNamespace` check, the script and smoke copies, and a build-time smoke asserting 21 lines of output.

Generate the smoke dataset with the Task 6 Step 3 generator into `containers/r-batch-metadict/smoke/`.

- [ ] **Step 6: Build**

Run: `docker build -f containers/r-batch-metadict/Dockerfile -t r-batch-metadict:dev .`
Expected: success, with the smoke `RUN` producing 21 lines.

- [ ] **Step 7: Register in CI and add the smoke case**

Matrix entry as before. Then in `tests/integration/test_batch_correct_smoke.py`, add `test_metadict_shrinks_batch_and_keeps_group`, copying `test_mmuphin_shrinks_batch_and_keeps_group` with `backend="metadict"` and `value_type == "relative"`.

- [ ] **Step 8: Run everything**

Run: `uv run pytest tests/test_batch_correct.py tests/test_batch_backends.py tests/test_container_skeletons.py -v`
Then: `MICROSUITE_RUN_BATCH_SMOKE=1 uv run pytest tests/integration/test_batch_correct_smoke.py -v`
Expected: all PASS. All five backends now have a real-execution case.

- [ ] **Step 9: Commit**

```bash
git add src/microsuite/batch/r/metadict.R containers/r-batch-metadict tests .github/workflows/docker.yml src/microsuite/batch/backends.py
git commit -m "feat(batch): add the metadict backend"
```

---

### Task 13: Documentation

**Files:**
- Create: `docs/batch_correction.md`
- Modify: `docs/methods.md`, `CHANGELOG.md`, `README.md`
- Test: `tests/test_batch_docs.py` (create)

**Interfaces:**
- Consumes: `batch.backends.BATCH_BACKENDS`.
- Produces: no code symbols.

- [ ] **Step 1: Write the failing doc test**

Create `tests/test_batch_docs.py`:

```python
from __future__ import annotations

from pathlib import Path

from microsuite.batch.backends import BATCH_BACKENDS

ROOT = Path(__file__).resolve().parents[1]


def test_every_backend_is_documented() -> None:
    text = (ROOT / "docs" / "batch_correction.md").read_text(encoding="utf-8")
    for name in BATCH_BACKENDS:
        assert name in text, f"{name} is not documented"


def test_the_leakage_hazard_is_stated() -> None:
    text = (ROOT / "docs" / "batch_correction.md").read_text(encoding="utf-8")
    assert "--target-col" in text
    assert "plsda-batch" in text
    # The document must say what goes wrong, not merely that an option exists.
    assert "inflat" in text.lower()


def test_the_scale_contract_is_documented() -> None:
    text = (ROOT / "docs" / "batch_correction.md").read_text(encoding="utf-8")
    for value_type in ("counts", "relative", "clr"):
        assert value_type in text


def test_correction_is_not_presented_as_a_substitute_for_modelling() -> None:
    text = (ROOT / "docs" / "batch_correction.md").read_text(encoding="utf-8")
    assert "covariate" in text.lower()
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_batch_docs.py -v`
Expected: FAIL — `FileNotFoundError: docs/batch_correction.md`.

- [ ] **Step 3: Write `docs/batch_correction.md`**

Cover, in this order:

1. **When correction is appropriate.** Multiple runs, platforms, extraction protocols, or cohorts merged into one table, with batch confounded with nothing you care about. Include the concrete case from the audit: 21 run accessions merged with `run_id` kept only as provenance.
2. **Why it is not a substitute for modelling batch as a covariate.** Correction adjusts the table; a covariate term adjusts the inference. Where the design permits it, model batch in `diff_abundance --fix-formula` instead of, or in addition to, correcting. State plainly that correcting once and no longer thinking about batch is the failure this section exists to prevent.
3. **The backend table** — the six-column table from the spec, plus the container image name for each.
4. **The output scale contract.** What `counts`, `relative`, and `clr` mean; which downstream commands refuse which; that an unmarked table is never refused. Show the exact error a user will see and what to do about it.
5. **Supervised backends and label leakage.** `plsda-batch` (and, later, `debias-m`) fit using the outcome labels. Correcting with the outcome and then testing that same outcome inflates significance. Say what to do instead: hold the correction out of the test, or use an unsupervised backend.
6. **How to check that it worked.** Until `batch diagnose` lands, point at `microsuite diversity test --method adonis2 --formula 'run_id + group'` before and after, and state the pair rule: batch R² must fall *and* group R² must hold. A correction that improves the first by destroying the second is a failure, not a success.

- [ ] **Step 4: Update `docs/methods.md`**

Add a `batch correct` section listing the five backends, their emitted scale, and their covariate/target requirements, linking to `docs/batch_correction.md`. Match the surrounding format.

- [ ] **Step 5: Update `CHANGELOG.md`**

Under the existing unreleased `## [0.3.0] - 2026-08-05` heading's `### Added`:

```markdown
- `microsuite batch correct` — batch effect correction with five backends:
  `mmuphin` (default), `combat-seq`, `conqur`, `plsda-batch`, and `metadict`,
  each in its own container image. Corrected tables record their scale in
  `uns["microsuite"]["value_type"]` as `counts`, `relative`, or `clr`.
- Count-requiring commands (`diff_abundance --backend ancombc/aldex2`,
  `rarefy`, `normalize`) now refuse tables whose recorded scale they cannot
  consume. Tables without a recorded scale are unaffected, so no existing
  pipeline changes behaviour.
```

- [ ] **Step 6: Update `README.md`**

Add `batch correct` to whatever command list or capability table the README already carries. Match the existing format; do not introduce a new section style.

- [ ] **Step 7: Run the doc tests and the whole suite**

Run: `uv run pytest tests/test_batch_docs.py -v`
Then: `uv run pytest --junitxml=/tmp/junit.xml` and read the counts out of the XML, since `-q` suppresses the summary line.

Expected: the only failures are the eleven documented Windows-host failures (`test_metadata_models`, `test_metadata_redact`, `test_run_fastp_multiqc_script`, `test_runtime_container`, `test_system_doctor`). Any other failure is yours.

Do not treat an unfamiliar failure as pre-existing without checking it against `git stash`. In the 2026-08-02 session a real `ruff format` failure was dismissed as baseline and CI stayed red for the whole session.

- [ ] **Step 8: Full lint and type check**

```bash
uv sync --all-extras --locked
uv run ruff format --check src tests
uv run ruff check src tests
uv run ty check
```

Expected: all clean.

- [ ] **Step 9: Commit**

```bash
git add docs/batch_correction.md docs/methods.md CHANGELOG.md README.md tests/test_batch_docs.py
git commit -m "docs: document batch correction, the scale contract, and the leakage hazard

Includes the rule that a correction is judged on a pair of numbers: the batch
effect must shrink and the biological signal must survive."
```

- [ ] **Step 10: Confirm CI is green**

Push the branch and check the run. Local green is not CI green — the eleven local failures pass on Linux, and the reverse case has happened too.

```bash
git push -u origin <branch>
```

Then check the Actions run for `ci.yml`. `docker.yml`'s heavy images build only on manual dispatch with `build-heavy-containers=true`; trigger that once and confirm all five new images build, since their build-time smokes are the only real-execution check that runs in CI.

---

## Self-Review

**Spec coverage:**

| Spec section | Task |
|---|---|
| Component 1 — backends and containers | 5, 6, 9, 10, 11, 12 |
| Component 1 — R script convention | 5 (`params.json`, four positional args) |
| Component 2 — three-valued `value_type` | 2 |
| Component 2 — the five downstream guards | 3 |
| Component 2 — supervised leakage warning | 4 (raise), 7 (CLI help), 13 (docs) |
| Component 4 — runner generalization | 1 |
| Component 4 — `resolve_batch_image` | 1 |
| Code layout | 1–7 (`diagnostics.py` and `prda_pvca.R` belong to the phase-4 plan) |
| Error handling table | 4 (option rejection), 5 (design validation), 3 (scale guards) |
| Testing — unit, mocked subprocess | 1, 2, 3, 4, 5, 7 |
| Testing — per-backend integration smoke | 8, 9, 10, 11, 12 |
| Testing — the deliberate-break check | 1 (step 8), 3 (step 7), 5 (step 9), 8 (steps 4–5) |
| Documentation | 13 |

Spec items deliberately absent, both assigned to the phase 3–4 plan: `debias-m`, and `batch diagnose` with its native metrics and `r-ecology` lme4 addition.

**Two gaps closed during review:**

- The spec's guard table named `normalize --method relative` but not `--method total-sum`. They are the same operation; the plan guards both and says why in the File Structure section.
- The spec did not say what happens when a backend drops features or returns them reordered. Task 5 rebuilds by label and raises on unknown features or dropped samples, with a deliberate-break step proving the test catches a positional rebuild.

**Fixture consistency note carried into Task 5:** the two-sample `_adata()` in the first draft made `sex` perfectly confounded with `run_id`, which would have made the covariate-passing test fail against the confounding guard added in the same task. Step 1 calls for a four-sample crossed fixture plus a separate two-sample `_confounded()` helper.

**Type consistency:** `invoke_r_script`'s keyword names in Task 1 match the call in Task 5 (`backend`, `script_package`, `script_name`, `resolve_image`, `positional`, `runtime`, `image`, `engine`, `run_dir`, `timeout`, `log`, `local_missing_message`). `BatchBackend` field names in Task 4 match their uses in Task 5 (`.name`, `.script`, `.package`, `.install_hint`, `.value_type`). `record_batch_correction`'s signature in Task 2 matches its call in Task 5. `run_batch_correction`'s signature in Task 5 matches its calls in Tasks 8–12 and in `methods/batch_correct.py`.

**Unverified external APIs, flagged rather than guessed:** ConQuR, PLSDAbatch, and MetaDICT are GitHub-only with no release tags. Tasks 10, 11, and 12 each begin by reading the real signature out of the built image and adjusting the script, rather than trusting the signature written into this plan. MetaDICT is the least certain of the three and gets an explicit probe-image step.
