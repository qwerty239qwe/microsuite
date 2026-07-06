# Reference-DB Subsystem (`refdb`) — Design

- **Date:** 2026-07-04
- **Status:** Approved (design), pending implementation plan
- **Author:** YT Lin (with Claude)
- **Origin:** Assessment of whether microsuite can replicate the
  [FOMC 16S rRNA Taxonomy Assignment Algorithm](https://github.com/tsute/FOMC_16S_rRNA_Taxonomy_Assignment_Algorithm).
  The assessment exposed that microsuite has **no reference-database management
  layer** — the single highest-leverage gap. This spec covers the foundational
  sub-project that unblocks the rest.

## Background & problem

microsuite's `tax_classify` accepts a `--classifier` argument that is a raw
filesystem path handed straight to a backend (`method_taxonomy_cmd.py:34`). There
is no acquisition, no build/format step, no versioning, no caching, and no
validation of reference databases. The `qiime2` backend is additionally
hardcoded to `classify-sklearn` (`methods/tax_classify.py:349`), so even the
reference-based consensus methods QIIME2 offers are unreachable.

Consequences observed against the FOMC pipeline:
- Cannot **acquire** a reference DB (no SILVA / GreenGenes2 / NCBI 16S / HOMD /
  MOMD download).
- Cannot **build** a classifier artifact from FASTA + taxonomy (what FOMC's
  combined HOMD+MOMD+NCBI reference requires).
- Cannot **register / version / cache / validate** a DB for reproducible reuse.

This is a missing subsystem, not a patch.

## Scope

This spec covers **sub-project A** of a four-part roadmap. Only A is fully
specified here; B–D are stubs to be brainstormed when reached.

| # | Sub-project | Depends on | Summary |
|---|---|---|---|
| **A** | **`refdb` subsystem** | — | Provider interface + registry to fetch, build, cache, and validate reference DBs. **This spec.** |
| B | Reference-based consensus taxonomy backend | A | Expose QIIME2 `classify-consensus-vsearch` / `classify-consensus-blast` at configurable identity/coverage (e.g. 98/98), the FOMC-signature method. |
| C | Workflow composability | A, B | Branchable pipeline steps (assign → unassigned → cluster → reassign) so FOMC-shaped custom pipelines can be assembled, not just fixed named workflows. |
| D | Robustness & testing pass | A, B | Backend-dispatch hardening and external-tool integration tests across the suite. The DB-relevant slice is folded into A; the broader pass is deferred. |

### Out of scope for A
- The consensus-taxonomy backend itself (sub-project B).
- Any change to how workflows are composed (sub-project C).
- Fetchers for DBs beyond the target amplicon set (SILVA, NCBI 16S, GTDB, UNITE,
  GG2, HOMD, MOMD); anything further follows the same provider pattern later.
- Non-amplicon biodbs work (its existing gene/protein/pathway coverage is
  untouched).

## Design

### Acquisition strategy (load-bearing decision)

Base QIIME2 has no fetch functions. The QIIME2 ecosystem's RESCRIPt plugin
(`q2-rescript`) can fetch SILVA/NCBI/GTDB/UNITE, but **biodbs is the default
acquisition provider** and will be extended to cover *all* reference DBs —
including SILVA, NCBI 16S, GTDB, UNITE, GreenGenes2, HOMD, and MOMD. RESCRIPt is
kept as an **optional alternate provider** for users who prefer it or already
have a QIIME2 environment, but it is not the default and is not required.

| Provider | Role | DBs | Notes |
|---|---|---|---|
| **`biodbs`** (default) | Primary acquisition for every DB | SILVA, NCBI 16S, GTDB, UNITE, GreenGenes2, **HOMD, MOMD** | Upstream biodbs gains amplicon-reference fetchers (this sub-project). Returns raw sequences + taxonomy; microsuite builds artifacts. |
| `rescript` (optional) | Alternate for RESCRIPt-covered DBs | SILVA, NCBI 16S, GTDB, UNITE, GG2 | Wrapped only when explicitly selected. Already emits `.qza`. Requires a QIIME2 env. |

> **As-built reconciliation (2026-07-06).** biodbs **v0.4.0** ships fetchers for
> **HOMD, SILVA, GTDB, GreenGenes, UNITE, and PR2** — this is what the delivered
> `biodbs` provider supports. It does **not** provide **MOMD** or an **NCBI 16S
> Microbial sequence** download (its NCBI functions are gene/taxonomy REST, not
> 16S RefSeq sequences). The provider's `_DB_ADAPTERS` therefore has exactly
> those six DBs, and an unsupported name (including `momd`/`ncbi-16s`) raises a
> clear error. Consequence: the flagship **FOMC "HOMD+MOMD+NCBI combined"
> reference cannot yet be produced through this provider** — the HOMD arm works;
> MOMD and NCBI 16S need a separate source (a future sub-project). Treat the
> MOMD/NCBI-16S entries above as aspirational, not delivered.

**biodbs division of labor.** biodbs
([github.com/qwerty239qwe/biodbs](https://github.com/qwerty239qwe/biodbs)) is a
REST-API fetch framework with a clean `fetch/<db>/` base-fetcher pattern
(`biodbs/fetch/_base.py`), today covering gene/protein/pathway/chemical DBs
(BioMart, KEGG, QuickGO, HPA, Ensembl, ChEMBL, UniProt, PubChem, …). This
sub-project extends biodbs upstream with amplicon-reference fetchers for the
full DB set above, following that existing pattern. biodbs stays
**acquisition-only** — it returns raw sequence + taxonomy files. The **build**
step (FASTA+taxonomy → vsearch reference / BLAST db / QIIME2 `.qza`) lives in
microsuite, because it is bioinformatics-tool invocation (`makeblastdb`,
`qiime tools import`), not REST fetching. Because biodbs is acquisition-only,
the microsuite `build.py` layer is the **primary** build path (the RESCRIPt
provider is the only one that can short-circuit to a pre-built `.qza`).

microsuite must stay usable when a provider's optional dependency is absent: a
raw `--classifier` path always works; selecting `--provider rescript` without a
QIIME2 env yields a clear error rather than a crash.

### Package layout

New package `src/microsuite/refdb/`:

```
refdb/
  __init__.py
  spec.py        # RefDbSpec: name, version, provider, target (16S/ITS/...), build_targets, sources[]
  registry.py    # JSON manifest cache: resolve/record built DBs by (name, version, build_target)
  providers/
    _base.py     # RefDbProvider: fetch(spec) -> RawRefDb ; build(raw, build_target) -> BuiltArtifact
    biodbs.py    # DEFAULT: all DBs via biodbs fetchers; normalize to FASTA + taxonomy TSV
    rescript.py  # OPTIONAL alternate: wraps q2-rescript get-silva-data/get-ncbi-data/... (+ classifier download)
  build.py       # RawRefDb -> vsearch reference | BLAST db (makeblastdb) | qiime2 .qza (qiime tools import); concat+dedup for merges
```

### Core types

- `RawRefDb` — local paths to a sequences FASTA + taxonomy TSV (RESCRIPt may
  hand back `.qza`; the RESCRIPt provider records those directly).
- `BuiltArtifact` — a built, backend-ready file plus its `build_target` tag and
  checksum.
- `RefDbSpec` — declarative description of a DB: `name`, `version`, `provider`,
  `target`, `build_targets`, and a `sources` list (>1 source ⇒ a merge, e.g.
  FOMC's HOMD+MOMD+NCBI).

### Provider interface

Both providers implement the same two methods (`providers/_base.py`):

- `fetch(spec) -> RawRefDb` — obtain sequences + taxonomy locally.
- `build(raw, build_target) -> BuiltArtifact` — produce a backend-ready artifact.
  The default `biodbs` provider always routes through `build.py`; the optional
  `rescript` provider already emits `.qza` and may short-circuit.

### Registry & caching

A JSON manifest under a cache directory (`$MICROSUITE_REFDB_DIR`, default
`~/.cache/microsuite/refdb/`). Entries keyed by `(name, version, build_target)`
with a stored checksum. Behavior:

- First fetch/build records the entry and artifact path.
- Subsequent resolves reuse the artifact; nothing re-downloads or rebuilds.
- Checksum mismatch triggers a rebuild.
- A corrupt/unreadable manifest yields a clear `MicrobiomeSuiteError`, not a
  crash.

This is the versioning/caching/validation the repo lacks today.

### CLI surface

```
microsuite refdb fetch <name> --version <v> \
    [--provider biodbs|rescript] --build <vsearch|blast|qiime2>
```

`--provider` defaults to `biodbs`; `rescript` is opt-in.

Prints the built artifact path and records it in the registry.
`tax_classify --classifier` accepts **either**:
- a raw filesystem path (today's behavior — unchanged, unbroken), **or**
- a registry reference such as `refdb:homd@15.22`, resolved via the registry.

### FOMC merge case

A `RefDbSpec` with multiple `sources` (HOMD v15.22 + MOMD v5.1 + NCBI 16S).
Providers fetch each source; `build.py` concatenates and de-duplicates into one
artifact registered as `refdb:fomc-combined@20221029`. This provides the
reference that sub-project B's consensus-taxonomy backend will align against at
98% identity / 98% coverage.

> **As-built (2026-07-06):** the multi-source merge *machinery* is delivered and
> tested (fixture-scale), but the specific FOMC combined DB is **not yet
> producible** because the delivered biodbs provider supplies HOMD but not MOMD
> or NCBI 16S sequences (see the reconciliation note under Acquisition
> strategy). The HOMD arm can be fetched today; wiring MOMD + NCBI 16S sources is
> a follow-up.

## Testing & robustness

Real fetches are large and network-bound, so tests must not pull real DBs.

- **Fake provider** — `FakeProvider` returning tiny bundled FASTA+taxonomy
  fixtures (~10 sequences). Exercises fetch → build → register → resolve
  end-to-end with no network. This is the bulk of coverage.
- **Tiny fixture reference DB** committed under
  `src/microsuite/data/fixtures/` (a ~10-sequence mock "HOMD-like" set) so
  `build.py` (concat/dedup + `makeblastdb` path) and the merge case get real
  assertions.
- **RESCRIPt / biodbs wrappers** — unit-test command construction (assert the
  argv we would invoke) without running the external tool, matching the existing
  backend pattern (`shutil.which` guard + `run_command`). Real-tool runs go
  behind the existing external-integration-test marker (opt-in).
- **Registry** — round-trip tests: record → resolve; checksum mismatch →
  rebuild; corrupt manifest → clear error.
- **Robustness** — every new failure path raises `MicrobiomeSuiteError` with an
  actionable message (missing tool, missing DB, bad version), matching the
  repo's current convention. This folds priority D's DB-relevant slice into A;
  the broader dispatch-hardening pass remains sub-project D.

## Success criteria

1. `microsuite refdb fetch` (default `--provider biodbs`) can build and register
   a DB from the biodbs lane, proven with the fixture set.
2. The optional `rescript` provider is proven by argv-assertion + opt-in
   integration; selecting it without a QIIME2 env gives a clear error.
3. The FOMC combined DB can be produced as a single registered artifact from a
   multi-source spec (fixture-scale in CI).
4. `tax_classify --classifier` resolves both a raw path and a `refdb:` reference,
   with the raw-path behavior unchanged.
5. Re-running a fetch reuses the cached artifact (no re-download/rebuild) and a
   checksum mismatch forces a rebuild.
6. Full unit suite runs offline; real-tool paths are opt-in.
7. biodbs gains amplicon-reference fetchers for the full DB set (SILVA, NCBI 16S,
   GTDB, UNITE, GG2, HOMD, MOMD) following its existing `fetch/<db>/` pattern.

## Open questions / dependencies

- biodbs needs upstream extension to fetch every amplicon reference DB (SILVA,
  NCBI 16S, GTDB, UNITE, GG2, HOMD, MOMD). That work is part of this sub-project
  and lands in the biodbs repo. Since it is the default provider, the microsuite
  side can be developed against the `FakeProvider` and fixtures while the biodbs
  fetchers land in parallel.
- Exact RESCRIPt subcommand flags (optional provider) to be pinned during
  implementation against the installed `q2-rescript` version.
