# First-class Docker runtime for the `dada2-r` backend (P1) — Design

- **Date:** 2026-07-08
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** Codex's ERP120510 friction list, complaint #1 (and #7): the R/DADA2
  wrapper assumes a local `Rscript`, forcing users to install R locally or bypass
  the wrapper and call the R script directly in a container. This is **P1** of a
  four-part roadmap (see [[dada2-improvement-roadmap]]) to address all outstanding
  DADA2 complaints. Complaints #2/#3/#4 were already fixed by PRs #7/#8.
- **Reference:** codex's dataset-verified `run_dada2_asv.sh` implements exactly
  this docker invocation and proves it works on ERP120510; P1 mirrors its mount
  pattern inside microsuite.

## Scope

P1 only: add a first-class Docker execution path to `microsuite denoise
--backend dada2-r`, implemented through a small reusable runtime helper so P2+
and other backends can adopt it later.

### Out of scope for P1
- Output validation (#8/#9) — that is **P2**.
- DADA2 Docker docs (#5/#6) — that is **P3**.
- Naming-contract end-to-end test (#10) — that is **P4**.
- Container execution for any backend other than `dada2-r`.
- Singularity/podman engines (the helper takes an `engine` param so they are a
  later addition, but P1 ships `docker` only).

## Reference invocation (from codex's working `run_dada2_asv.sh`)

```bash
docker run --rm \
  -v "${abs_input}:/input:ro" \
  -v "${abs_asv}:/asv" \
  -v "${abs_plots}:/plots" \
  -v "${abs_script}:/dada2_denoise.R:ro" \
  "${IMAGE}" /dada2_denoise.R \
  --input-dir /input --output-table /asv/asv_table.tsv \
  --output-rep-seqs /asv/rep_seqs.fasta --output-stats /asv/denoising_stats.tsv \
  --output-plot-dir /plots --threads N [--paired ...]
```

Key facts this confirms: the `r-dada2` image's entrypoint is `Rscript`, so the
mounted script path is the first inner argument; inputs are mounted read-only,
outputs/plots read-write; the host script is mounted (not baked) so the running
R matches the installed microsuite; threads are a resolved integer before the
container starts.

## Design

### Component 1 — reusable helper `src/microsuite/runtime/container.py`

```python
@dataclass(frozen=True)
class Mount:
    host: Path        # resolved absolute host path
    container: str    # absolute path inside the container
    mode: str = "rw"  # "ro" or "rw"

def build_container_command(
    inner: list[str], image: str, mounts: list[Mount], *, engine: str = "docker"
) -> list[str]:
    """Return e.g. ['docker','run','--rm','-v','H:C:ro',...,'IMAGE', *inner]."""

def require_engine(engine: str = "docker") -> str:
    """shutil.which(engine); raise MicrobiomeSuiteError with an actionable
    install message if absent."""

def resolve_dada2_image(override: str | None) -> str:
    """override, else env MICROSUITE_R_DADA2_IMAGE, else the default GHCR ref
    'ghcr.io/qwerty239qwe/microsuite/r-dada2:latest'."""
```

`build_container_command` emits `engine run --rm` followed by one `-v
host:container[:ro]` per mount (mode omitted for `rw` to match docker's
default), then the image, then the inner argv. Mount host paths are always
resolved to absolute before mounting.

A `PathMapper` (in the same module) assigns each distinct host directory a
stable container mountpoint and rewrites a host path to its container path:

```python
class PathMapper:
    def add_dir(self, host_dir: Path, mode: str, container: str) -> None: ...
    def to_container(self, host_path: Path) -> str:  # "<mount>/<name>"
    def mounts(self) -> list[Mount]:                 # deduped
```

### Component 2 — `denoise_dada2_r` gains `runtime` + `image`

Signature adds `runtime: str = "local"` and `image: str | None = None`. When
`runtime == "local"`, behavior is byte-for-byte today's (Rscript on PATH). When
`runtime == "docker"`:

1. `require_engine("docker")` (skip the `shutil.which("Rscript")` check — R is
   in the image, not on the host).
2. Materialize the packaged R script to a real filesystem path via
   `importlib.resources.as_file(files("microsuite.methods.r").joinpath(
   DADA2_R_SCRIPT))` (handles the zip/namespace case), inside the run's temp
   scope.
3. Build a `PathMapper`:
   - `input_dir` → `/work/input` **ro**
   - script's dir → `/work/script` **ro** (script becomes `/work/script/dada2_denoise.R`)
   - the parent dir of each of `output_table`, `output_rep_seqs`,
     `output_stats`, and (if set) `output_plot_dir` → `/work/outN` **rw**
     (deduplicated: outputs sharing a parent share one mount)
4. Build the same inner Rscript argv as local mode, but with every path argument
   passed through `PathMapper.to_container(...)`. The inner argv starts with the
   container script path (the image's `Rscript` entrypoint consumes it).
5. `command = build_container_command(inner, resolve_dada2_image(image),
   mapper.mounts(), engine="docker")`.
6. Run through the existing `run_command(...)` (timeout/run_dir/log preserved),
   so provenance logging and the `CommandLog(backend="dada2-r")` are unchanged.

`_prepare_outputs(...)` still runs on the host first; because the host output
parent dirs are mounted rw and the container writes into them, outputs land at
the caller's original host paths — identical to local mode.

**Threads:** `resolve_threads(threads)` already yields an int on the host before
the argv is built, so no in-container `auto` resolution is needed.

**Symlink caveat (from the reference script):** a bind mount does not expose the
targets of symlinks that live *inside* `input_dir`. If `input_dir` contains
symlinked FASTQs, the container sees dangling links. P1 detects this: before
building the docker command, if any FASTQ directly under `input_dir` is a
symlink, raise a `MicrobiomeSuiteError` explaining that docker mode needs real
files (copy/hardlink them, or use `--runtime local`). This is a guard, not
silent breakage.

### Component 3 — CLI surface & error messages (fixes #7)

`method_features_cmd.py` `denoise` command gains, threaded through
`denoise()` → `denoise_dada2_r()`:
- `--runtime local|docker` (default `local`).
- `--image TEXT` (optional; overrides the default `r-dada2` image / env).

`--runtime docker` is only valid for `--backend dada2-r`; combined with any other
backend it raises `MicrobiomeSuiteError` (matching the existing backend-specific
option guards in `denoise.py`).

Error messages:
- Local mode, `Rscript` missing — append: *"Or run it in a container with
  `--runtime docker` (uses the r-dada2 image; no local R needed). See
  docs/dada2.md."* (fixes #7).
- Docker mode, engine missing — `require_engine` raises an actionable message
  naming `docker` and pointing to install docs.
- Docker mode, image pull/run failure — surfaced by the existing `run_command`
  failure path (`"R/DADA2 denoising failed."`), with the full docker argv in the
  run log.

## Testing

All offline (monkeypatch `shutil.which` + `subprocess.run`), matching the repo's
backend-test convention:

1. **docker argv** — `denoise_dada2_r(..., runtime="docker")` builds
   `docker run --rm -v <input>:/work/input:ro -v <script-dir>:/work/script:ro
   -v <out-parent>:/work/out0 ... <default-image> /work/script/dada2_denoise.R
   --input-dir /work/input --output-table /work/out0/<name> ...`: correct engine,
   `--rm`, ro/rw modes, deduped output mounts, and every path arg rewritten.
2. **image override** — `--image` and `MICROSUITE_R_DADA2_IMAGE` change the image
   in the argv; default is the GHCR `r-dada2` ref otherwise (`resolve_dada2_image`
   precedence: override → env → default).
3. **local unchanged** — `runtime="local"` still builds the `Rscript ...` argv and
   still errors (with the new #7 message) when Rscript is absent.
4. **guards** — docker-missing-engine message names `docker`; a symlinked FASTQ
   under `input_dir` in docker mode raises with an actionable message;
   `--runtime docker` + non-dada2-r backend raises.
5. **helper units** — `build_container_command` emits correct `-v` flags/order;
   `PathMapper` dedups shared parents and rewrites paths; `resolve_dada2_image`
   honors the precedence.

Real `docker run` execution stays opt-in behind the external-integration marker
(needs Docker + the image); offline tests fully cover argv/mount/guard/error
logic.

## Success criteria

1. `microsuite denoise --backend dada2-r --runtime docker` builds and runs the
   r-dada2 container invocation (mirroring the reference script), writing the
   ASV table / rep-seqs / stats / plots to the caller's host paths.
2. `--runtime local` (default) is byte-for-byte unchanged; the local
   Rscript-missing error now points to `--runtime docker`.
3. `--image` / `MICROSUITE_R_DADA2_IMAGE` override the image; default is the GHCR
   r-dada2 ref.
4. Symlinked-input and missing-docker cases fail with actionable errors, not
   silent breakage.
5. `container.py` is a reusable helper (Mount/build_container_command/PathMapper/
   resolve/require_engine), covered by offline unit tests; the full suite stays
   green.

## Open questions / follow-ups (not blocking P1)

- Singularity/podman engines (helper already parameterized by `engine`).
- Whether P2's output validation should run identically for local and docker
  runs (it should — validation is post-run and runtime-agnostic).
- A future `microsuite-dada2` whole-CLI-in-container path (`docker run
  microsuite-dada2 denoise ...`) is complementary and belongs to P3's docs.
