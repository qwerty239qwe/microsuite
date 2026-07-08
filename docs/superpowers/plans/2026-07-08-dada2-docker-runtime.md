# DADA2 Docker Runtime (P1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `--runtime local|docker` to `microsuite denoise --backend dada2-r` so the R/DADA2 step can run in the `r-dada2` container (no local R), via a reusable `runtime/container.py` helper, mirroring the dataset-verified `run_dada2_asv.sh` invocation.

**Architecture:** A new `runtime/container.py` provides a `Mount` dataclass, `build_container_command`, `PathMapper`, `require_engine`, and `resolve_dada2_image`. `denoise_dada2_r` gains `runtime`/`image`; in docker mode it mounts input (ro), output/plot parents (rw), and the packaged R script (ro), rewrites path args to container paths, and runs the same argv through the existing `run_command`. Local mode is unchanged.

**Tech Stack:** Python 3.12, `subprocess` via the repo's `run_command`, Docker CLI (invoked, never run in unit tests — argv is asserted with monkeypatched `shutil.which`/`subprocess.run`).

## Global Constraints

- `runtime="local"` (default) is byte-for-byte today's behavior; the change is non-breaking.
- Inputs mounted read-only (`:ro`), outputs/plots read-write, the host R script read-only — matching the reference `run_dada2_asv.sh`.
- Default image `ghcr.io/qwerty239qwe/microsuite/r-dada2:latest`; override precedence: `--image` → env `MICROSUITE_R_DADA2_IMAGE` → default.
- `--runtime docker` is only valid for `--backend dada2-r`; any other backend raises `MicrobiomeSuiteError`.
- All failure paths raise `MicrobiomeSuiteError` with actionable messages (missing engine; symlinked input under docker; local Rscript-missing must now mention `--runtime docker`).
- Threads are already an `int` (via `resolve_threads`) before the docker argv is built — do not resolve `auto` in the container.
- `from __future__ import annotations` at the top of every new module. External tools are argv-asserted offline; real `docker run` stays opt-in behind the external-integration marker.

---

### Task 1: `runtime/container.py` reusable helper

**Files:**
- Create: `src/microsuite/runtime/container.py`
- Test: `tests/test_runtime_container.py`

**Interfaces:**
- Produces:
  - `Mount(host: Path, container: str, mode: str = "rw")` — frozen dataclass.
  - `build_container_command(inner: list[str], image: str, mounts: list[Mount], *, engine: str = "docker") -> list[str]`.
  - `require_engine(engine: str = "docker") -> str` — `shutil.which`; raises `MicrobiomeSuiteError` if absent.
  - `resolve_dada2_image(override: str | None) -> str` — override → env `MICROSUITE_R_DADA2_IMAGE` → `DEFAULT_DADA2_IMAGE`.
  - `class PathMapper` with `add_dir(host_dir, mode, container)`, `container_dir(host_dir) -> str`, `to_container(host_path) -> str`, `mounts() -> list[Mount]`.
  - `DEFAULT_DADA2_IMAGE = "ghcr.io/qwerty239qwe/microsuite/r-dada2:latest"`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_runtime_container.py
from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.runtime.container import (
    DEFAULT_DADA2_IMAGE,
    Mount,
    PathMapper,
    build_container_command,
    require_engine,
    resolve_dada2_image,
)


def test_build_container_command_argv() -> None:
    mounts = [
        Mount(Path("/h/in"), "/work/input", "ro"),
        Mount(Path("/h/out"), "/work/out0", "rw"),
    ]
    cmd = build_container_command(
        ["/work/script/x.R", "--input-dir", "/work/input"], "img:tag", mounts
    )
    assert cmd == [
        "docker", "run", "--rm",
        "-v", "/h/in:/work/input:ro",
        "-v", "/h/out:/work/out0",
        "img:tag",
        "/work/script/x.R", "--input-dir", "/work/input",
    ]


def test_build_container_command_engine_override() -> None:
    cmd = build_container_command(["x"], "img", [], engine="podman")
    assert cmd[:3] == ["podman", "run", "--rm"]


def test_resolve_image_precedence(monkeypatch) -> None:
    monkeypatch.delenv("MICROSUITE_R_DADA2_IMAGE", raising=False)
    assert resolve_dada2_image(None) == DEFAULT_DADA2_IMAGE
    monkeypatch.setenv("MICROSUITE_R_DADA2_IMAGE", "env:img")
    assert resolve_dada2_image(None) == "env:img"
    assert resolve_dada2_image("override:img") == "override:img"


def test_require_engine_missing(monkeypatch) -> None:
    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(MicrobiomeSuiteError):
        require_engine("docker")


def test_pathmapper_dedup_and_rewrite(tmp_path: Path) -> None:
    inp = tmp_path / "in"
    out = tmp_path / "out"
    inp.mkdir()
    out.mkdir()
    mapper = PathMapper()
    mapper.add_dir(inp, "ro", "/work/input")
    mapper.add_dir(out, "rw", "/work/out0")
    mapper.add_dir(out, "rw", "/work/out0")  # duplicate -> one mount
    assert mapper.container_dir(inp) == "/work/input"
    assert mapper.to_container(out / "table.tsv") == "/work/out0/table.tsv"
    assert len(mapper.mounts()) == 2


def test_pathmapper_upgrades_ro_to_rw(tmp_path: Path) -> None:
    d = tmp_path / "d"
    d.mkdir()
    mapper = PathMapper()
    mapper.add_dir(d, "ro", "/work/d")
    mapper.add_dir(d, "rw", "/work/d")  # a writer appears -> upgrade
    assert mapper.mounts()[0].mode == "rw"
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_runtime_container.py -v`
Expected: FAIL (`ModuleNotFoundError: microsuite.runtime.container`).

- [ ] **Step 3: Create `src/microsuite/runtime/container.py`**

```python
from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError

DEFAULT_DADA2_IMAGE = "ghcr.io/qwerty239qwe/microsuite/r-dada2:latest"
_DADA2_IMAGE_ENV = "MICROSUITE_R_DADA2_IMAGE"


@dataclass(frozen=True)
class Mount:
    host: Path
    container: str
    mode: str = "rw"


def build_container_command(
    inner: list[str], image: str, mounts: list[Mount], *, engine: str = "docker"
) -> list[str]:
    command = [engine, "run", "--rm"]
    for mount in mounts:
        spec = f"{mount.host}:{mount.container}"
        if mount.mode == "ro":
            spec += ":ro"
        command.extend(["-v", spec])
    command.append(image)
    command.extend(inner)
    return command


def require_engine(engine: str = "docker") -> str:
    resolved = shutil.which(engine)
    if resolved is None:
        raise MicrobiomeSuiteError(
            f"The '{engine}' container engine is required for --runtime docker but was "
            f"not found on PATH. Install {engine}, or use --runtime local with R and the "
            "dada2 package installed."
        )
    return resolved


def resolve_dada2_image(override: str | None) -> str:
    if override:
        return override
    env = os.environ.get(_DADA2_IMAGE_ENV)
    if env:
        return env
    return DEFAULT_DADA2_IMAGE


class PathMapper:
    """Assign host directories stable container mountpoints and rewrite paths."""

    def __init__(self) -> None:
        self._dirs: dict[Path, Mount] = {}

    def add_dir(self, host_dir: Path, mode: str, container: str) -> None:
        resolved = host_dir.resolve()
        existing = self._dirs.get(resolved)
        if existing is None:
            self._dirs[resolved] = Mount(host=resolved, container=container, mode=mode)
        elif existing.mode == "ro" and mode == "rw":
            self._dirs[resolved] = Mount(
                host=resolved, container=existing.container, mode="rw"
            )

    def container_dir(self, host_dir: Path) -> str:
        return self._dirs[host_dir.resolve()].container

    def to_container(self, host_path: Path) -> str:
        resolved = host_path.resolve()
        mount = self._dirs.get(resolved.parent)
        if mount is None:
            raise MicrobiomeSuiteError(
                f"No container mount registered for {resolved.parent}"
            )
        return f"{mount.container}/{resolved.name}"

    def mounts(self) -> list[Mount]:
        return list(self._dirs.values())
```

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest tests/test_runtime_container.py -v`
Expected: PASS (6 passed).

- [ ] **Step 5: Commit**

```bash
git add src/microsuite/runtime/container.py tests/test_runtime_container.py
git commit -m "feat(runtime): reusable container-command helper (Mount/build/PathMapper/resolve)"
```

---

### Task 2: `dada2-r` docker runtime + CLI + error rework

**Files:**
- Modify: `src/microsuite/methods/denoise.py` (extract `_dada2_r_script_args`; add `runtime`/`image` to `denoise_dada2_r` and `denoise`; docker branch; guards; #7 error)
- Modify: `src/microsuite/cli/method_features_cmd.py` (add `--runtime`/`--image`, pass through)
- Test: `tests/test_denoise_cluster_methods.py` (add docker-runtime tests)

**Interfaces:**
- Consumes: `build_container_command`, `PathMapper`, `require_engine`, `resolve_dada2_image` (Task 1).
- Produces: `denoise_dada2_r(..., runtime: str = "local", image: str | None = None)`; `denoise(..., runtime: str = "local", dada2_image: str | None = None)`; CLI `--runtime`/`--image`.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_denoise_cluster_methods.py`)

```python
def test_denoise_dada2_r_docker_builds_container_command(
    tmp_path, monkeypatch
) -> None:
    import subprocess
    from microsuite.methods.denoise import denoise

    input_dir = tmp_path / "reads"
    input_dir.mkdir()
    (input_dir / "s_1.fastq.gz").write_text("x")
    (input_dir / "s_2.fastq.gz").write_text("x")
    out = tmp_path / "out"
    out.mkdir()

    monkeypatch.setattr("shutil.which", lambda name: "docker" if name == "docker" else None)
    calls: list[list[str]] = []

    def fake_run(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    denoise(
        backend="dada2-r",
        demux=input_dir,
        output_table=out / "table.tsv",
        output_rep_seqs=out / "rep.fasta",
        output_stats=out / "stats.tsv",
        output_plot_dir=out / "plots",
        mode="paired",
        threads=2,
        force=True,
        runtime="docker",
    )

    cmd = calls[0]
    assert cmd[0] == "docker" and cmd[1] == "run" and "--rm" in cmd
    assert any(tok.endswith("dada2_denoise.R") for tok in cmd)
    # input mounted ro, outputs rw
    assert any(v.endswith(":ro") and "reads" in v for v in cmd)
    # default GHCR image present
    assert any("r-dada2" in tok for tok in cmd)
    # container paths, not host paths, in the R args
    assert "--input-dir" in cmd
    idx = cmd.index("--input-dir")
    assert cmd[idx + 1].startswith("/work/")


def test_denoise_dada2_r_docker_image_override(tmp_path, monkeypatch) -> None:
    import subprocess
    from microsuite.methods.denoise import denoise

    input_dir = tmp_path / "reads"
    input_dir.mkdir()
    (input_dir / "s.fastq.gz").write_text("x")
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setattr("shutil.which", lambda name: "docker" if name == "docker" else None)
    calls: list[list[str]] = []
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kw: (calls.append(command), subprocess.CompletedProcess(command, 0, "", ""))[1],
    )
    denoise(
        backend="dada2-r", demux=input_dir,
        output_table=out / "t.tsv", output_rep_seqs=out / "r.fa", output_stats=out / "s.tsv",
        mode="single", threads=1, force=True, runtime="docker", dada2_image="myrepo/rd:1.2",
    )
    assert "myrepo/rd:1.2" in calls[0]


def test_denoise_dada2_r_docker_missing_engine(tmp_path, monkeypatch) -> None:
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import denoise

    input_dir = tmp_path / "reads"
    input_dir.mkdir()
    (input_dir / "s.fastq.gz").write_text("x")
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(MicrobiomeSuiteError, match="docker"):
        denoise(
            backend="dada2-r", demux=input_dir,
            output_table=out / "t.tsv", output_rep_seqs=out / "r.fa", output_stats=out / "s.tsv",
            mode="single", threads=1, force=True, runtime="docker",
        )


def test_denoise_dada2_r_local_error_points_to_docker(tmp_path, monkeypatch) -> None:
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import denoise

    input_dir = tmp_path / "reads"
    input_dir.mkdir()
    (input_dir / "s.fastq.gz").write_text("x")
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setattr("shutil.which", lambda name: None)  # no Rscript
    with pytest.raises(MicrobiomeSuiteError, match="--runtime docker"):
        denoise(
            backend="dada2-r", demux=input_dir,
            output_table=out / "t.tsv", output_rep_seqs=out / "r.fa", output_stats=out / "s.tsv",
            mode="single", threads=1, force=True,  # runtime defaults to local
        )


def test_denoise_docker_rejected_for_non_dada2r_backend(tmp_path) -> None:
    from microsuite._errors import MicrobiomeSuiteError
    from microsuite.methods.denoise import denoise

    with pytest.raises(MicrobiomeSuiteError, match="dada2-r"):
        denoise(
            backend="qiime2-dada2", demux=tmp_path / "d.qza",
            output_table=tmp_path / "t.qza", output_rep_seqs=tmp_path / "r.qza",
            output_stats=tmp_path / "s.qza", trunc_len=150, runtime="docker",
        )
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_denoise_cluster_methods.py -k docker -v`
Expected: FAIL (`denoise()` has no `runtime` param / behavior absent).

- [ ] **Step 3: Extract `_dada2_r_script_args` in `denoise.py`**

Add this helper (near `_append_value`/`_append_bool`). It reproduces the exact R-arg construction currently inline in `denoise_dada2_r`, but takes path **strings** so both runtimes reuse it:

```python
def _dada2_r_script_args(
    *,
    input_dir: str,
    output_table: str,
    output_rep_seqs: str,
    output_stats: str,
    output_plot_dir: str | None,
    threads: int,
    paired: bool,
    tuning: Dada2Tuning,
    max_n: int | None,
    rm_phix: bool | None,
) -> list[str]:
    args = [
        "--input-dir", input_dir,
        "--output-table", output_table,
        "--output-rep-seqs", output_rep_seqs,
        "--output-stats", output_stats,
        "--threads", str(threads),
    ]
    if output_plot_dir is not None:
        args.extend(["--output-plot-dir", output_plot_dir])
    if paired:
        args.append("--paired")
        args.extend([
            "--trim-left-f", str(tuning.trim_left_f),
            "--trunc-len-f", str(tuning.trunc_len_f),
            "--trim-left-r", str(tuning.trim_left_r),
            "--trunc-len-r", str(tuning.trunc_len_r),
        ])
        _append_value(args, "--max-ee-f", tuning.max_ee_f)
        _append_value(args, "--max-ee-r", tuning.max_ee_r)
        _append_value(args, "--min-overlap", tuning.min_overlap)
        _append_value(args, "--max-merge-mismatch", tuning.max_merge_mismatch)
        _append_bool(args, "--trim-overhang", tuning.trim_overhang)
    else:
        args.extend(["--trim-left", str(tuning.trim_left), "--trunc-len", str(tuning.trunc_len)])
        _append_value(args, "--max-ee", tuning.max_ee)
    _append_value(args, "--trunc-q", tuning.trunc_q)
    _append_value(args, "--max-n", max_n)
    _append_bool(args, "--rm-phix", rm_phix)
    _append_value(args, "--pooling-method", tuning.pooling_method)
    _append_value(args, "--chimera-method", tuning.chimera_method)
    _append_value(args, "--min-fold-parent-over-abundance", tuning.min_fold_parent_over_abundance)
    _append_bool(args, "--allow-one-off", tuning.allow_one_off)
    _append_value(args, "--n-reads-learn", tuning.n_reads_learn)
    return args
```

> Verify against the current inline construction in `denoise_dada2_r` that the flags, order, and value sources (all from `tuning`, plus `max_n`/`rm_phix`) match exactly; this helper must produce byte-for-byte the same args the current local path produces.

- [ ] **Step 4: Rewrite the body of `denoise_dada2_r`** (from the `rscript = shutil.which("Rscript")` block through the `_run(...)` call) to add `runtime`/`image` and the docker branch. Change the signature to add `runtime: str = "local"` and `image: str | None = None` (keyword-only, after `timeout`). New body:

```python
    from importlib.resources import as_file
    from microsuite.runtime.container import (
        PathMapper,
        build_container_command,
        require_engine,
        resolve_dada2_image,
    )

    if not input_dir.exists() or not input_dir.is_dir():
        raise MicrobiomeSuiteError(f"Input directory does not exist: {input_dir}")
    _prepare_outputs(output_table, output_rep_seqs, output_stats, force=force)
    if output_plot_dir is not None:
        output_plot_dir.mkdir(parents=True, exist_ok=True)

    script_res = files("microsuite.methods.r").joinpath(DADA2_R_SCRIPT)

    if runtime == "docker":
        require_engine("docker")
        # bind mounts cannot expose symlink targets inside input_dir
        for entry in input_dir.iterdir():
            if entry.is_symlink() and entry.name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
                raise MicrobiomeSuiteError(
                    f"Input FASTQ is a symlink and cannot be mounted into the container: "
                    f"{entry}. Copy/hardlink real files, or use --runtime local."
                )
        with as_file(script_res) as script_path:
            mapper = PathMapper()
            mapper.add_dir(input_dir, "ro", "/work/input")
            mapper.add_dir(script_path.parent, "ro", "/work/script")
            out_index = 0
            for out in (output_table, output_rep_seqs, output_stats):
                key = out.resolve().parent
                if key not in {m.host for m in mapper.mounts()}:
                    mapper.add_dir(key, "rw", f"/work/out{out_index}")
                    out_index += 1
            if output_plot_dir is not None:
                if output_plot_dir.resolve() not in {m.host for m in mapper.mounts()}:
                    mapper.add_dir(output_plot_dir, "rw", f"/work/out{out_index}")
            inner = [f"{mapper.container_dir(script_path.parent)}/{script_path.name}"]
            inner += _dada2_r_script_args(
                input_dir=mapper.container_dir(input_dir),
                output_table=mapper.to_container(output_table),
                output_rep_seqs=mapper.to_container(output_rep_seqs),
                output_stats=mapper.to_container(output_stats),
                output_plot_dir=(mapper.container_dir(output_plot_dir) if output_plot_dir else None),
                threads=threads, paired=paired, tuning=tuning, max_n=max_n, rm_phix=rm_phix,
            )
            command = build_container_command(
                inner, resolve_dada2_image(image), mapper.mounts(), engine="docker"
            )
            _run(command, "R/DADA2 denoising failed.", run_dir=run_dir, timeout=timeout,
                 backend="dada2-r", inputs={"input_dir": str(input_dir)},
                 outputs={"table": str(output_table), "representative_sequences": str(output_rep_seqs),
                          "denoising_stats": str(output_stats),
                          **({"plot_dir": str(output_plot_dir)} if output_plot_dir is not None else {})},
                 params=_dada2_log_params(mode="paired" if paired else "single", tuning=tuning,
                                          max_n=max_n, rm_phix=rm_phix))
        return

    rscript = shutil.which("Rscript")
    if rscript is None:
        raise MicrobiomeSuiteError(
            "R/DADA2 denoising requires the external 'Rscript' command. "
            "Install R with the dada2 package, or run it in a container with "
            "--runtime docker (uses the r-dada2 image; no local R needed). See docs/dada2.md."
        )
    with as_file(script_res) as script_path:
        command = [rscript, str(script_path)] + _dada2_r_script_args(
            input_dir=str(input_dir),
            output_table=str(output_table),
            output_rep_seqs=str(output_rep_seqs),
            output_stats=str(output_stats),
            output_plot_dir=(str(output_plot_dir) if output_plot_dir else None),
            threads=threads, paired=paired, tuning=tuning, max_n=max_n, rm_phix=rm_phix,
        )
        _run(command, "R/DADA2 denoising failed.", run_dir=run_dir, timeout=timeout,
             backend="dada2-r", inputs={"input_dir": str(input_dir)},
             outputs={"table": str(output_table), "representative_sequences": str(output_rep_seqs),
                      "denoising_stats": str(output_stats),
                      **({"plot_dir": str(output_plot_dir)} if output_plot_dir is not None else {})},
             params=_dada2_log_params(mode="paired" if paired else "single", tuning=tuning,
                                      max_n=max_n, rm_phix=rm_phix))
```

> Match the existing `_run(...)` arguments (backend/inputs/outputs/params) to whatever the current code passes — the block above shows the shape; align the `params=_dada2_log_params(...)` call to the current signature. The two `_run` calls (local/docker) are identical except for `command`; factor them into one call after the if/else if that reads cleaner.

- [ ] **Step 5: Thread `runtime`/`image` through `denoise()`** — add `runtime: str = "local"` and `dada2_image: str | None = None` to the `denoise(...)` signature (keyword-only). In the `dada2-r` dispatch branch, pass `runtime=runtime, image=dada2_image`. Add a guard near the top of `denoise` (after `backend` is known): if `runtime != "local"` and `backend != "dada2-r"`, raise `MicrobiomeSuiteError("--runtime docker is only supported for --backend dada2-r.")`.

- [ ] **Step 6: Wire the CLI** — in `src/microsuite/cli/method_features_cmd.py` `denoise_cmd`, add two options and pass them to `denoise(...)`:

```python
        runtime: Annotated[
            str, typer.Option("--runtime", help="dada2-r execution: local or docker.")
        ] = "local",
        image: Annotated[
            str | None,
            typer.Option("--image", help="Container image for --runtime docker (dada2-r)."),
        ] = None,
```
and in the `denoise(...)` call add `runtime=runtime, dada2_image=image,`.

- [ ] **Step 7: Run to verify pass + no regression**

Run: `uv run pytest tests/test_denoise_cluster_methods.py -v`
Expected: PASS (new docker tests + the existing dada2_r local tests unchanged).

- [ ] **Step 8: Commit**

```bash
git add src/microsuite/methods/denoise.py src/microsuite/cli/method_features_cmd.py tests/test_denoise_cluster_methods.py
git commit -m "feat(dada2): --runtime docker for dada2-r backend + actionable Rscript error"
```

---

## Self-Review

**Spec coverage:**
- Reusable `container.py` (Mount/build/PathMapper/require_engine/resolve) → Task 1. ✓
- `denoise_dada2_r` docker branch mirroring the reference invocation (input ro, outputs rw, host script ro, container path rewrite) → Task 2 Step 4. ✓
- Local unchanged + #7 error points to `--runtime docker` → Task 2 Step 4 (local branch) + test. ✓
- image override precedence → Task 1 `resolve_dada2_image` + Task 2 test. ✓
- guards: missing engine, symlinked input, non-dada2-r backend → Task 2 Steps 4/5 + tests. ✓
- threads pre-resolved (int) → `denoise` passes `resolved_threads`; helper stringifies. ✓
- offline argv tests; real docker opt-in → Task 2 tests monkeypatch subprocess. ✓

**Placeholder scan:** the two "match the current construction / `_run` args" notes are verification directives against existing verified code, not missing content — the new code (helper, container.py, docker branch, CLI, guards, tests) is complete. The implementer must diff `_dada2_r_script_args` against the current inline args to guarantee byte-for-byte parity.

**Consistency:** `runtime`/`image` names align: CLI `--runtime`/`--image` → `denoise(runtime=, dada2_image=)` → `denoise_dada2_r(runtime=, image=)`; helper names (`build_container_command`, `PathMapper.container_dir`/`to_container`, `resolve_dada2_image`, `require_engine`) match Task 1 exactly.
