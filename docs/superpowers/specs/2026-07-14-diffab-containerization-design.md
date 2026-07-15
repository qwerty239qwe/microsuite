# Differential-abundance containerization (Round-4 I) — Design

- **Date:** 2026-07-14
- **Status:** Approved (design), pending implementation plan
- **Origin:** Round-4 complaints **#3** (all four diffab R backends shell out to a
  host `Rscript`; no `--runtime`/`--image`, no image provenance, bind-mounted
  outputs would be root-owned) and **#4** (the single `r-diffab` image is
  all-or-nothing: one broken bioconductor package breaks all four backends, and
  the build "verifies" only `requireNamespace` + `file.exists`, never a real run).
  Third diffab sub-project (**I**); follows H (model expressiveness). See
  [[microsuite-round4-roadmap]].

## Scope

Route all four R differential-abundance backends (ancombc, aldex2, maaslin2,
lefse) through the P1 container infra so each can run either on a host `Rscript`
(`--runtime local`, default, unchanged behavior) or inside a per-backend Docker
image (`--runtime docker`), with the container executed as the caller's UID/GID
(writable bind-mounted outputs) and the resolved image + digest recorded beside
the output. Replace the brittle combined `r-diffab` image with four per-backend
images, each verified at build time by running a tiny real input through its
`.R` script (not merely importing the package).

### Out of scope for I
- Standardized cross-backend result contract (#5) — **J**.
- Formula parity for aldex2/maaslin2/lefse (repeated measures) — later follow-up.
- Containerizing the QIIME 2 ANCOM-BC path (`qiime2-ancombc`) — already a
  separate qiime image; unaffected here.
- Exposing a `--engine podman` switch — engine stays `docker` (the infra already
  parameterizes `engine`; not surfaced as a user option under YAGNI).

## Verified current state

- `runtime/container.py`: `build_container_command(inner, image, mounts, *, engine)`
  (no `--user`), `PathMapper`, `require_engine` (message is dada2-specific),
  `resolve_dada2_image`. No UID/GID helper, no diffab image resolver, no digest
  helper.
- `diffab/ancombc.py` (`run_ancombc`) and `diffab/r_backends.py`
  (`run_r_diffab_backend`) each `shutil.which("Rscript")`, write TSVs (+
  `params.json` for ancombc) into a `TemporaryDirectory`, and `run_command`
  `[rscript, script, ...args, output]`. No runtime/image parameters.
- `methods/diff_abundance.py` dispatches to those two wrappers; CLIs
  `cli/diffab_cmd.py` (`ancombc`) and `cli/method_stats_cmd.py` (`diff_abundance`)
  expose no runtime/image options.
- `containers/r-diffab/Dockerfile`: installs all four bioconductor packages in one
  image, `COPY`s all four `.R`, `ENTRYPOINT ["Rscript"]`; build check is
  `requireNamespace` only. CI (`docker.yml`) builds it as one heavy image; its
  smoke step runs `requireNamespace` + `file.exists`.
- `tests/test_container_skeletons.py` enumerates `r-diffab` and asserts
  `docs/methods.md` links `[R diffab](../containers/r-diffab/Dockerfile)`.

## Design

### Component 1 — `runtime/container.py` extensions

- `build_container_command(inner, image, mounts, *, engine="docker", user=None)`:
  when `user` is truthy, insert `["--user", user]` immediately after `--rm`.
  Existing callers (dada2) pass no `user` → argv unchanged.
- `current_user_spec() -> str | None`: return `f"{os.getuid()}:{os.getgid()}"`
  when both `os.getuid`/`os.getgid` exist, else `None` (Windows → no `--user`).
- `resolve_diffab_image(backend, override) -> str`: `override` wins; else env
  `MICROSUITE_R_DIFFAB_<BACKEND>_IMAGE` (backend upper-cased); else default
  `ghcr.io/qwerty239qwe/microsuite/r-diffab-<backend>:latest`. Add
  `DEFAULT_DIFFAB_IMAGE_PREFIX` + `_DIFFAB_IMAGE_ENV_PREFIX` constants.
- `resolve_image_digest(engine, image) -> str | None`: `shutil.which(engine)`;
  if absent return `None`. Run `<engine> inspect --format {{index .RepoDigests 0}}`
  then fall back to `--format {{.Id}}`; return the first non-empty stripped
  stdout, else `None`. Swallow `OSError`/`SubprocessError`/non-zero exit →
  `None` (provenance is best-effort, never fatal).
- Generalize `require_engine`'s message: drop "the dada2 package", say "use
  `--runtime local` with R and the required package installed" (backend-agnostic;
  the existing test only asserts it raises).

### Component 2 — shared runner (`diffab/_runner.py`, new)

`invoke_r_backend(*, backend, positional, runtime="local", image=None,
engine="docker", run_dir=None, timeout=None, log, local_missing_message) -> None`

- `positional: list[str | Path]` is the exact ordered argument list the `.R`
  script expects. **Convention: the last `Path` is the output (rw-mounted parent);
  every earlier `Path` is an input (ro-mounted parent); `str` items pass through
  verbatim.** Both wrappers already end with the output Path, so this rule holds.
- `runtime` must be `"local"` or `"docker"` (else `MicrobiomeSuiteError`).
- **local:** `shutil.which("Rscript")`; if `None`, raise
  `MicrobiomeSuiteError(local_missing_message)`. Resolve the script via
  `files("microsuite.diffab.r").joinpath(f"{backend}.R")`. Run
  `[rscript, str(script), *[str(a) for a in positional]]`.
- **docker:** `resolve_diffab_image(backend, image)`; `require_engine(engine)`.
  Build a `PathMapper`: register each unique input-parent dir `ro` and the
  output-parent dir `rw` (a shared dir is upgraded to `rw` by `PathMapper`),
  assigning deterministic mountpoints (`/mnt/d0`, `/mnt/d1`, …). Inner command is
  `[f"/opt/microsuite/{backend}.R", *(mapper.to_container(a) if Path else a)]`
  (image `ENTRYPOINT` is `Rscript`). Assemble via `build_container_command(...,
  engine=engine, user=current_user_spec())`.
- Run through `run_command(command, failure_message=f"{backend} failed.",
  run_dir=run_dir, log=log, timeout=timeout)`.
- **docker only, after the run:** write `<output-parent>/<backend>_container.json`
  = `{"runtime": "docker", "engine": engine, "image": resolved_image,
  "digest": resolve_image_digest(engine, resolved_image)}` (uniform provenance
  for all four backends; #3). `local` writes no sidecar.

### Component 3 — wrappers route through the runner

- `diffab/ancombc.py` `run_ancombc(...)`: add `runtime="local"`, `image=None`,
  `engine="docker"` params. Keep the obs-column pre-check. Drop the inline
  `shutil.which` (the runner owns it). After writing counts/metadata/params.json,
  call `invoke_r_backend(backend="ancombc", positional=[counts_path,
  metadata_path, params_path, output], runtime=runtime, image=image,
  engine=engine, run_dir=run_dir, timeout=timeout, log=<current CommandLog>,
  local_missing_message=<current ANCOM-BC Rscript message>)`.
- `diffab/r_backends.py` `run_r_diffab_backend(...)`: add the same three params.
  Keep the group-column check. Drop the inline `shutil.which`. Call
  `invoke_r_backend(backend=backend, positional=[counts_path, metadata_path,
  group, output], runtime=runtime, image=image, engine=engine, ...,
  local_missing_message=<current per-package message>)`.

### Component 4 — method + CLI threading

- `methods/diff_abundance.py` `diff_abundance(...)`: add `runtime="local"`,
  `image=None`, `engine="docker"`; forward to `run_r_diffab_backend` and
  `run_ancombc`. The `qiime2-ancombc` branch is unaffected.
- `cli/method_stats_cmd.py` `diff_abundance_cmd`: add `--runtime` (default
  `"local"`) and `--image` options; forward to `diff_abundance`.
- `cli/diffab_cmd.py` `ancombc`: add `--runtime` (default `"local"`) and
  `--image`; forward to `run_ancombc`.

### Component 5 — four per-backend images (replace the combined one)

Delete `containers/r-diffab/`. Create `containers/r-diffab-<backend>/Dockerfile`
for `ancombc`, `aldex2`, `maaslin2`, `lefse`, each:
- `FROM mambaorg/micromamba:1.5.10`, the two OCI labels, `# Expected commands:`
  comment, `ENV PATH`.
- Install pinned `"r-base>=4.3,<4.4"` + **only that backend's** bioconductor
  package (`bioconductor-ancombc` / `-aldex2` / `-maaslin2` / `-lefser`), then
  `micromamba clean -a -y`.
- `COPY src/microsuite/diffab/r/<backend>.R /opt/microsuite/<backend>.R`.
- `COPY containers/r-diffab-<backend>/smoke/ /opt/microsuite/smoke/` — a tiny
  fixture (counts.tsv + metadata.tsv, plus params.json for ancombc) sized so the
  method actually **converges** (ANCOM-BC2 in particular needs enough
  samples/taxa; the smoke fixture is a real, if minimal, dataset, not a 1-row
  stub).
- **Build-time smoke:** `RUN Rscript /opt/microsuite/<backend>.R
  /opt/microsuite/smoke/... /tmp/microsuite-smoke-out.tsv && test -s
  /tmp/microsuite-smoke-out.tsv` (echo diagnostics before any cleanup). The
  smoke output belongs under `/tmp` because the micromamba base runs as an
  unprivileged user and `/opt/microsuite` is root-owned after `COPY`. A broken
  package fails the build, per-image, and only for that backend (#4).
- `ENTRYPOINT ["Rscript"]`.

Because the smoke run is a `RUN` step, a successful `docker build` already proves
the backend runs end-to-end; CI needs no separate smoke step.

### Component 6 — CI, docs, skeleton test

- `.github/workflows/docker.yml`: replace the single `r-diffab` matrix entry with
  four (`r-diffab-ancombc`, `r-diffab-aldex2`, `r-diffab-maaslin2`,
  `r-diffab-lefse`), each `heavy: true`; remove the old "Smoke test R diffab
  image" step (build-time `RUN` smoke replaces it).
- `docs/methods.md`: replace the `[R diffab](../containers/r-diffab/Dockerfile)`
  link with the four per-backend links; note each backend is one image and runs
  via `--runtime docker --image`.
- `tests/test_container_skeletons.py`: replace the `r-diffab` entry with the four
  `r-diffab-<backend>` entries (each with its own expected tokens — `Rscript` +
  the one package) and update the `docs/methods.md` link assertions to match.

## Testing (all offline / unit)

1. **`tests/test_runtime_container.py`** (extend): `build_container_command(...,
   user="1000:1000")` puts `--user 1000:1000` right after `--rm`; omitting `user`
   leaves argv unchanged. `current_user_spec()` returns `"<uid>:<gid>"`
   (monkeypatch `os.getuid`/`os.getgid`) and `None` when they're absent.
   `resolve_diffab_image`: default per backend, env override, explicit override.
   `resolve_image_digest`: monkeypatch `subprocess.run` → `RepoDigests[0]`, the
   `.Id` fallback when RepoDigests is empty, and `None` on missing engine /
   non-zero / raised error.
2. **`tests/test_diffab_runner.py`** (new): `invoke_r_backend` with a stubbed
   `run_command` (capture argv). `runtime="local"` → `[rscript, script, *args]`
   and no sidecar. `runtime="docker"` → argv is `docker run --rm --user u:g -v
   ...:ro -v ...:rw <image> /opt/microsuite/<backend>.R <mapped args>`, inputs
   `ro`, output parent `rw`, and `<backend>_container.json` written with the
   resolved image + digest (monkeypatch `resolve_image_digest`). Invalid
   `runtime` raises; missing local `Rscript` raises with the supplied message.
3. **`tests/test_diffab_ancombc.py`** / **`test_diff_abundance_method.py`**
   (update): the params.json contract is unchanged, but the command assertion now
   goes through the runner — assert `run_ancombc`/`run_r_diffab_backend` still
   produce the same local argv (default `runtime="local"`), and that
   `--runtime docker --image X` reaches the runner. Add a CLI smoke
   (`--runtime docker --image ...`) via a monkeypatched backend for both
   `diffab ancombc` and `method diff_abundance`.
4. **`tests/test_container_skeletons.py`**: the four new Dockerfiles exist with
   the right labels/tokens; the combined `r-diffab` path no longer referenced.

**Not offline-testable (build-time / CI):** the four Dockerfile smoke runs need
`docker build`; verified in CI, not in the pytest suite.

## Success criteria

1. All four backends accept `--runtime local|docker` and `--image`; `local`
   (default) is byte-for-byte the current host-`Rscript` behavior.
2. `--runtime docker` runs the per-backend image as the caller's UID/GID
   (bind-mounted output stays caller-owned) and writes
   `<backend>_container.json` with the resolved image reference + digest.
3. The combined `r-diffab` image is gone; four per-backend images each fail their
   own build if their bioconductor package can't actually run a tiny real input.
4. Full suite + `ty check` + ruff + format green; container-skeleton + methods
   docs tests updated to the four images.

## Open questions / follow-ups (not blocking I)
- Publishing the four image tags to GHCR (push workflow) mirrors the existing
  dada2 image publication; wiring the release/publish job is a deploy follow-up.
- Once I lands, `diffab`'s `_container.json` digest feeds naturally into J's
  standardized result contract (provenance columns).
- The oral pipeline's bespoke MaAsLin2 Docker detour can retire once H+I+J land
  (see [[oral-pipeline-microsuite-refactor]]).
