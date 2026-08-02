# microsuite audit — evidence from `oral_microbiome_gingivitis`

**Date:** 2026-08-02
**Subject repo:** `C:\Users\qwert\PycharmProjects\oral_microbiome_gingivitis` (private; 39 shell wrappers, 73 Python helpers, 74 TOML configs, 29 docs, ~21 study accessions)
**Method:** every workaround script, container, and written complaint in the consumer repo traced back to the microsuite code or doc it works around.

Paths beginning `scripts/`, `docs/`, `configs/`, `projects/`, `nextflow/` are in the **gingivitis** repo. Paths beginning `src/`, `containers/`, `docs/methods.md` are in **microsuite**.

---

## 1. Executive summary — five highest-value findings

1. **`diff_abundance --backend maaslin2` runs MaAsLin2 with `normalization = "NONE"` on a raw count matrix.** Library-size variation is never removed, so every coefficient is depth-contaminated. It returns a complete, well-formed, wrong result. The project silently replaced it with `normalization = "TSS"`. (`src/microsuite/diffab/r/maaslin2.R:31` vs `scripts/r/run_maaslin2_repeated_measures.R:67`)
2. **The R differential-abundance family accepts exactly one fixed effect and no random effects**, so a repeated-measures study cannot use it at all. The project wrote a 76-line R script, a dedicated Dockerfile, and four Python helpers to replace it. Meanwhile `microsuite diffab ancombc` *does* expose `--fix-formula`/`--rand-formula` — the capability exists but is unreachable from `diff_abundance`. (`src/microsuite/methods/diff_abundance.py:17-29`, `src/microsuite/diffab/r_backends.py:22-33`)
3. **`simpson` means Gini-Simpson (1 − D) and `simpson_d` means D, and no microsuite doc says so.** `grep -rn simpson docs/` in microsuite returns zero hits. vegan and mothur use the opposite convention. A user who reads the metric name gets a silently inverted diversity ordering. The project maintains a hand-written remap table. (`src/microsuite/diversity/alpha.py:454,458` vs `scripts/python/finalize_microsuite_diversity.py:13-14`)
4. **Provenance (`--run-dir`) is wired into 9 of ~22 CLI modules, and none of the native ones.** `diversity_cmd`, `table_cmd`, `viz_cmd`, `import_cmd`, `ordination_cmd`, `ml_cmd` emit nothing. The project rebuilt the entire run-bundle layer itself. (`src/microsuite/cli/`, `docs/workflow_reporting.md:1-31`)
5. **Native beta diversity is Bray-Curtis and Jaccard only.** Sørensen, Aitchison, and tree-based UniFrac all had to be sourced elsewhere; UniFrac exists in microsuite only behind QIIME 2 artifacts. The project built a whole container to get them from scikit-bio. (`src/microsuite/diversity/beta.py:9`, `docs/ORDINATION_ANALYSIS.md:5-7`)

**Finding counts:** Correctness 4 · Missing capability 6 · Ergonomics 6 · Reproducibility 3 (19 total, plus 3 already-resolved items recorded in §6).

---

## 2. Correctness

> Weighted highest per the brief. Every item here produces output that looks valid.

### C1. `maaslin2` backend disables normalization on raw counts
- **Evidence:** `src/microsuite/diffab/r/maaslin2.R:31` — `normalization = "NONE"`, with `transform = "LOG"` on line 32. The matrix it receives is raw counts: `src/microsuite/diffab/r_backends.py:41-45` writes `dense_counts(adata)` straight to `counts.tsv` with no transformation.
- **Project contradiction:** `scripts/r/run_maaslin2_repeated_measures.R:67` — `normalization = "TSS"`.
- **Label:** Inferred. The project chose TSS deliberately and never wrote down "microsuite's setting is wrong", but the divergence is unambiguous and TSS is MaAsLin2's own default.
- **What microsuite does:** log-transforms unnormalized counts. Sequencing depth enters every feature's model as an unmodelled additive term. On a study where depth correlates with any covariate — batch, run, extraction date — the result is confounded, and MaAsLin2 will not warn.
- **What the project needed:** TSS normalization, plus the ability to choose `normalization`, `transform`, `analysis_method`, `min_prevalence`, `min_abundance` at all (their script exposes all five; microsuite exposes none).
- **Reach:** Every user of `--backend maaslin2` who does not pre-rarefy. High.

### C2. `simpson` / `simpson_d` naming collision, undocumented
- **Evidence:** `src/microsuite/diversity/alpha.py:454` `simpson = 1.0 - dominance(...)`; `:458` `simpson_d = dominance(...)`. `scripts/python/finalize_microsuite_diversity.py:13-14` remaps microsuite `simpson_d` → project `simpson`, and microsuite `simpson` → project `gini_simpson`.
- **Contradicts:** nothing — that is the problem. `docs/methods.md` lists the native alpha backend (line 136) without enumerating metrics, and no microsuite doc mentions Simpson at all.
- **Label:** Documented (the remap exists in code); the ambiguity itself is Inferred.
- **What microsuite does:** follows the scikit-bio convention (`simpson` = 1 − D) while vegan, mothur, and most microbial-ecology papers use `simpson` = D. Both are "correct" in some lineage; neither is written down. The failure is silent and directional — the two quantities are monotonically inverted, so a misread flips the entire conclusion about which group is more diverse.
- **What the project needed:** unambiguous names (`gini_simpson`, `simpson_d`) or a documented convention.
- **Reach:** Everyone who requests `simpson`. High, and the cost of being wrong is high.

### C3. Zero-depth samples produce well-formed NaN/0 rows
- **Evidence:** `src/microsuite/diversity/alpha.py:64` applies `_clean_counts` per sample; `:497-504` drops all non-positive entries, leaving an empty vector for an all-zero sample. `shannon` then returns `np.nan` (`:443-444`) while `sobs` returns `0.0` (`:472-474`). No warning, no error.
- **Project workaround:** `scripts/python/prepare_microsuite_metadata.py:23-46` — `filter_empty_samples()` strips zero-total samples before microsuite sees them and writes an `excluded_samples.tsv` recording `reason = zero_total_count`.
- **Label:** Inferred. The helper's existence and its dedicated exclusion-reason column are the evidence.
- **What microsuite does:** emits a row per sample mixing NaN and 0 across metrics, which flows into downstream means, plots, and tests as either a silent NaN-propagation or a real-looking zero.
- **What the project needed:** microsuite to refuse, or at minimum flag, zero-depth samples.
- **Reach:** Anyone subsetting by sample list, filtering aggressively, or rarefying. Medium.

### C4. R-backend metadata alignment can yield all-NA rows instead of an error
- **Evidence:** `src/microsuite/diffab/r/maaslin2.R:21` — `metadata <- metadata[colnames(counts), , drop = FALSE]`. R row-name indexing with an unmatched name returns an all-`NA` row rather than raising.
- **Project contradiction:** `scripts/r/run_maaslin2_repeated_measures.R:35-36` uses `match()` and then `stop("Metadata does not cover every count-table sample.")`.
- **Label:** Inferred.
- **Mitigating context:** in the normal path `r_backends.py:41-46` writes both files from the same AnnData, so they always align. The exposure is for anyone invoking the bundled R script directly, or any future path where the two inputs diverge.
- **Reach:** Low today, but it is a latent silent-failure mode in shipped code, and the same defensive check is missing across the aldex2/lefse scripts.

---

## 3. Missing capability

### M1. No multiple fixed effects or random effects in `diff_abundance`
- **Evidence (documented, verbatim):** `docs/oral_microbiome_work_report_2026-07-10.md:276-278` — *"The general microsuite differential-abundance CLI did not expose multiple fixed effects or subject-level random effects. A dedicated MaAsLin2 Docker image was used to support the longitudinal design without requiring host R installation."*
- **microsuite side:** `src/microsuite/methods/diff_abundance.py:17-29` — the signature carries `group: str` and nothing else. `src/microsuite/diffab/r_backends.py:22-33` — identical. `docs/methods.md:156` concedes it: *"current wrapper exposes the selected group as a single fixed effect."*
- **The asymmetry that matters:** `src/microsuite/cli/diffab_cmd.py:25` — the dedicated `ancombc` command *does* take `--fix-formula`, `--rand-formula`, `--reference`, `--prv-cut`, `--global`, `--pairwise`, `--trend`, `--dunnet`. So microsuite already knows how to express these models; the capability is simply absent from `diff_abundance` and from every R backend. The project uses both surfaces side by side in one script (`scripts/run_differential_abundance.sh:272-296` for hand-rolled MaAsLin2, `:311-340` for microsuite ANCOM-BC2).
- **What the project needed:** `abundance ~ time_day + (1 | subject_code)` and `abundance ~ phase_code + (1 | subject_code)`, at four ranks, times two models — eight runs (`docs/oral_microbiome_work_report_2026-07-10.md:283-309`). Plus, notably, **design-matrix diagnostics**: their script had to add its own rank-deficiency guard (`scripts/r/run_maaslin2_repeated_measures.R:55-58`) after the combined time+phase model turned out to be structurally confounded. microsuite offers no such check.
- **Also hand-rolled around it:** `scripts/python/list_diffab_models.py`, `list_diffab_methods.py`, `list_ancombc_models.py`, `summarize_diffab_results.py`, `summarize_ancombc_results.py`, `plot_diffab_publication.py`, `plot_differential_abundance.py`, `plot_ancombc_results.py`, `scripts/containers/Dockerfile.maaslin2`, `configs/ERP120510/ERP120510_differential_abundance_relaxed_human_oral_weighted.toml`.
- **Reach:** Every longitudinal, paired, or multi-covariate design. Very high — this is the single most common real-world microbiome design microsuite cannot serve.

### M2. No cross-run / cross-study feature-table merge
- **Evidence:** `scripts/python/combine_ordination_tables.py` (unions tables by feature ID, refuses duplicate sample IDs rather than auto-prefixing, records SHA-256 per source) and `nextflow/ordination.nf`. `docs/NEXTFLOW_ORDINATION.md:82-84` describes combining **21 run tables** into one bundle.
- **microsuite side:** no merge/concat entry exists in the Table Transforms section of `docs/methods.md:121-130`.
- **Label:** Documented (by the doc and the dedicated workflow).
- **What the project needed:** union by feature ID, sample-ID collision detection, a manifest of which run contributed each sample, and source hashing.
- **Reach:** Any meta-analysis or multi-accession study. High.

### M3. No batch-effect correction of any kind
- **Evidence:** `docs/HOMD_ABUNDANCE_QC.md:179-181` — a high-distance QC flag *"can result from disease, cohort biology, primer region, DNA extraction, sequencing, taxonomy naming, or other batch effects"*; `:187` — *"Primer, platform, cohort, site, and taxonomy-version differences remain important confounders."*
- **What they actually did about it: nothing.** A repo-wide grep for `combat|MMUPHin|batch.?correct|harmoniz` returns zero hits outside those two advisory sentences. They merge 21 runs (M2) and correct none of it; `run_id` survives only as a provenance column in `inputs/combined_sample_manifest.tsv`, usable as a plot grouping for visual inspection.
- **Label:** Documented that the confounding exists; Inferred that microsuite offers no remedy (confirmed absent from `docs/methods.md`).
- **What the project needed:** at minimum, batch as a covariate in the differential-abundance model (blocked by M1); ideally a correction method (MMUPHin/ComBat-seq) or a restricted-permutation beta test blocking on study.
- **Reach:** High for multi-study work; this is the gap that most limits what their 21-run merge can conclude.

### M4. Native beta metrics are Bray-Curtis and Jaccard only
- **Evidence:** `src/microsuite/diversity/beta.py:9` — `BETA_METRICS = {"bray-curtis", "jaccard"}`; `:12-24` raises on anything else.
- **Project consequences, three separate workarounds:**
  - Sørensen: available only as long-form pairs from `beta-turnover`; `scripts/python/finalize_microsuite_diversity.py:49-61` pivots it into a square matrix and fills the zero diagonal. The refactor spec calls this out explicitly — `docs/superpowers/specs/2026-07-12-scripts-microsuite-refactor-design.md:86-87`: *"**[keep]** (microsuite `BETA_METRICS` = bray-curtis, jaccard only; Sørensen stays bespoke)."*
  - Aitchison and UniFrac: `scripts/python/run_ordination_analysis.py:25-26` imports `skbio.diversity.beta_diversity` and `scipy.spatial.distance` directly.
  - `docs/ORDINATION_ANALYSIS.md:5-7`: *"It uses the MicroSuite API for Bray-Curtis, Jaccard, and PCoA, and scikit-learn for PCA and t-SNE. When a compatible rooted Newick tree is supplied, scikit-bio calculates weighted and unweighted UniFrac."*
- **UniFrac specifically:** `docs/methods.md:137-138` offers UniFrac only through `qiime2-diversity-lib` / `qiime2-core-metrics-phylogenetic`, i.e. requiring QIIME artifacts. The project wanted the opposite and said so — `docs/ORDINATION_ANALYSIS.md:81`: *"This helper uses the Newick tree directly and does not require QIIME 2 artifacts for UniFrac."*
- **Label:** Documented.
- **Reach:** High. Aitchison is the default for compositional workflows; UniFrac is standard in 16S.

### M5. No PCA, no t-SNE
- **Evidence:** `scripts/python/run_ordination_analysis.py:26-27` imports `sklearn.decomposition.PCA` and `sklearn.manifold.TSNE`. `docs/ORDINATION_ANALYSIS.md:5-6` and `:68-69` (outputs `ordination/pca.tsv`, `ordination/tsne.tsv`). `scripts/containers/Dockerfile.ordination:29` installs `scikit-learn` specifically for this.
- **microsuite side:** `docs/methods.md:142-144` offers constrained ordination (RDA/CCA/db-RDA) and PCoA, but no unconstrained PCA and no non-linear embedding.
- **Label:** Documented.
- **Reach:** Medium-high — PCA on CLR coordinates is the standard compositional ordination.

### M6. No alpha-diversity plots and no distance heatmap
- **Evidence (documented, explicit):** `docs/superpowers/specs/2026-07-12-scripts-microsuite-refactor-design.md:99-101` — `alpha_diversity_summary.png`, `alpha_depth_scatter.png`, `alpha_shannon_by_time.png`, `alpha_observed_asvs_by_time.png`, `alpha_metrics_by_phase.png` → *"**[keep]** (no microsuite alpha-plot command)"*; `:102-103` — `beta_bray_curtis_heatmap.png` and the baseline-dissimilarity boxplots → *"**[keep]** (no distance-heatmap / bespoke baseline-boxplot command)."*
- **Notable:** this spec is a deliberately *aggressive* migration ("replace both the table/transform logic AND the hand-drawn plots with new microsuite commands where an equivalent exists", `:5-8`). These four are what survived the aggression — they are the residue of genuine absence, not of reluctance. `README.md:210-211` confirms the outcome: *"The bespoke alpha plots, the Bray-Curtis distance heatmap, and the baseline-dissimilarity boxplots are kept (no microsuite equivalent)."*
- **Reach:** Medium. Alpha-by-group is the most-produced figure in the field.

---

## 4. Ergonomics / friction

### E1. Alpha diversity is one subprocess and one file per metric
- **Evidence:** `scripts/run_diversity_metrics.sh:175-178` loops eight metrics, each a separate `microsuite diversity alpha` invocation writing `microsuite_raw/<metric>.tsv`. `scripts/python/finalize_microsuite_diversity.py` then exists solely to merge those eight files — its module docstring (`:2`) reads *"Merge microsuite diversity outputs into the workflow's stable TSV contract."*
- **Compounding:** microsuite emits no sample-depth column, so the finalizer recomputes `total_reads` from the count table itself (`:22-37`, `:94`), and renames all eight metric columns on the way out (`:10-19`).
- **microsuite side:** `src/microsuite/cli/diversity_cmd.py` — `alpha` takes a single `--metric` and writes one file. Eight invocations means eight h5ad loads.
- **Label:** Documented (docstring + the file's existence).
- **Reach:** High — every user producing more than one alpha metric.

### E2. `beta-turnover` output shape differs from every other beta output
- **Evidence:** `scripts/python/finalize_microsuite_diversity.py:49-61` — `write_sorensen()` reads long-form `sample_a`/`sample_b`/`sorensen_dissimilarity` rows and pivots them into the square matrix that `diversity beta` would have produced natively.
- **Label:** Inferred.
- **Reach:** Medium — anyone combining turnover with the other beta outputs.

### E3. `import tsv` requires a pre-aligned, `sample_id`-keyed metadata TSV
- **Evidence:** `scripts/python/prepare_microsuite_metadata.py` — ~200 lines whose entire job is producing `metadata_for_microsuite.tsv`. It resolves the sample-ID column (with fallback and a warning when the configured column is wrong), rejects duplicate IDs, hard-errors on samples missing from metadata, drops zero-total samples (C3), and renames QIIME-reserved headers. That last one carries its own comment: *"QIIME 2 reserves these headers for the ID column; a *non-ID* metadata column that uses one of them makes `qiime` reject the entire metadata file"* (`:48-52`, sets at `:53-58`).
- **Called from:** `scripts/run_diversity_metrics.sh:151-154`, `scripts/run_differential_abundance.sh:317-319`, and elsewhere — it is a standing prerequisite for touching microsuite at all.
- **Label:** Inferred, with the QIIME-reserved-name portion documented in-code.
- **What microsuite could do:** accept a sample-ID column name, align by it, and give a targeted error naming the missing samples — instead of requiring the caller to guarantee the contract.
- **Reach:** High. Every project bridging SRA metadata to microsuite hits this.

### E4. `--front` accepts only one literal primer
- **Evidence (documented):** `docs/reports/pipeline_refactor_and_srp090878_report_20260710.md:88` — *"Because different samples have the V7-V9 amplicon read from either end, and cutadapt/microsuite's `--front` only accepts one literal primer, used cutadapt's `file:` adapter syntax (`configs/SRP090878_v7v9_front_primers.fasta`, both primers) so cutadapt tries both per read and uses whichever matches."*
- **Note:** they got through by smuggling cutadapt's own `file:` syntax through microsuite's string option. That works but is undiscoverable and undocumented.
- **Reach:** Medium — mixed-orientation and multi-primer libraries; more common in older 454 and pooled-amplicon deposits than people expect.

### E5. ANCOM-BC2 output is passed through raw, with version-dependent column names
- **Evidence:** `scripts/python/summarize_ancombc_results.py:17-27` — `effect_columns()` carries the docstring *"Map ANCOM-BC2 effects to q-value columns across package versions"* and handles both `q_val_*` and `q_*` prefixes; `:30-33` `result_value()` handles `modern`/`legacy` field-name pairs.
- **What microsuite does:** writes ANCOM-BC2's wide output verbatim, with no schema normalization and no pinned ANCOMBC version. `docs/methods.md:31` puts ANCOMBC in "External runtimes and databases remain user supplied."
- **Label:** Documented (docstring).
- **Reach:** Medium, but it turns any ANCOMBC upgrade into a silent downstream breakage.

### E6. R backends split inconsistently between host `Rscript` and Docker
- **Evidence:** `docs/SRP090878_END_TO_END.md:104-107` — *"The default differential-abundance stage includes ANCOM-BC2. With the current MicroSuite API that backend still requires host `Rscript` plus the R packages `ANCOMBC` and `jsonlite`; MaAsLin2 itself runs in Docker. This requirement can be removed when MicroSuite's planned per-backend Docker runtime is released."* Repeated at `docs/ERP120510_END_TO_END.md:92-94`.
- **microsuite side:** `docs/methods.md:153` advertises the ancombc backend as *"per-backend R/Bioconductor image, run via `--runtime docker`"*, and `containers/r-diffab-ancombc/` exists. The project's own script has `--ancombc-runtime` and `--ancombc-image` flags (`scripts/run_differential_abundance.sh:42-43`).
- **Assessment:** the capability landed, but the field experience and the doc disagree, and the project still documents the host-Rscript requirement as current. Either the docker path is incomplete for ancombc or the docs oversell it. Worth verifying end to end.
- **Reach:** Medium.

---

## 5. Reproducibility

### R1. `--run-dir` provenance covers 9 of ~22 CLI modules, and none of the native ones
- **Evidence:** modules referencing `run_dir` in `src/microsuite/cli/`: `diffab_cmd`, `method_diversity_cmd`, `method_features_cmd`, `method_preprocess_cmd`, `method_qiime_io_cmd`, `method_stats_cmd`, `method_tables_cmd`, `method_taxonomy_cmd`, `network_cmd`. **Absent from:** `diversity_cmd.py`, `table_cmd.py`, `viz_cmd.py`, `import_cmd.py`, `ordination_cmd.py`, `ml_cmd.py`, `data_cmd.py`, `workflow_cmd.py` — which is precisely the set this project runs most (`import tsv`, `table normalize`, `diversity alpha/beta/beta-turnover/beta-significance`, `viz`).
- **Project consequence:** they rebuilt the layer. `docs/workflow_reporting.md:1-16` defines their own immutable bundle (`workflow.json`, `resolved_config.json`, `run_manifest.json`, `stages/*.json`), backed by `scripts/python/workflow_reporting.py` (319 lines) and `scripts/python/workflow_registry.py`.
- **Documented, verbatim:** `docs/workflow_reporting.md:22-31` — *"Once MicroSuite metadata-A is implemented, its commands will attach their per-command envelopes to the same logical run without changing analysis output paths… Custom stages (FOMC, SMDI, plots, and MaAsLin2) are therefore represented before equivalent MicroSuite APIs exist."*
- **Also:** `docs/methods.md:189` advertises `microsuite report --backend native` as producing *"HTML provenance reports from run metadata"* — but for the native commands there is no run metadata to report on.
- **Reach:** Very high. This is the difference between a reproducible run and a directory of loose TSVs.

### R2. Required APIs live on unreleased branches; the project ships a version-detection shim
- **Evidence:** `scripts/python/microsuite_compat.py:10-18` declares a `REQUIRED_RUN_BUNDLE_API` of seven symbols and, when any is missing, raises (`:35-41`): *"MicroSuite checkout lacks the Nextflow run-bundle API: {…}. Update MicroSuite to a revision containing 'feat(metadata): add versioned run bundle contracts' (currently origin/dev/nextflow-run-bundles at ad86c9a or later)."*
- **Corroborating:** `docs/SRP051201_END_TO_END.md:19-21` — pyrosequencing settings *"require a MicroSuite version containing the corresponding CLI options."*
- **microsuite side:** `pyproject.toml:3` — `version = "0.2.0"`. There is no release cadence or feature-version gate a consumer can pin against, so the consumer pins a **git SHA on a dev branch** and writes runtime `hasattr` checks.
- **Label:** Documented.
- **Reach:** High for anyone tracking microsuite closely; it makes "which microsuite produced this result" unanswerable from the version string.

### R3. microsuite is consumed as a source checkout, never as a package or published image
- **Evidence:** `README.md:44-45` clones the microsuite repo alongside the project. Every wrapper needs `microsuite_dir` and invokes via `uv run --project "${MICROSUITE_DIR}" microsuite …` (`scripts/run_diversity_metrics.sh:120-129`, `:162`; `scripts/run_differential_abundance.sh:320-323`). `scripts/containers/Dockerfile.ordination:4-6` — *"The local MicroSuite checkout is supplied as a named BuildKit context, so this image does not require access to a private GHCR image"* — with `:22-31` copying `/pyproject.toml`, `/src`, `/README.md`, `/uv.lock` from that context.
- **Nuance:** GHCR images *do* exist and are used where available (`scripts/run_differential_abundance.sh:75` defaults to `ghcr.io/qwerty239qwe/microsuite/r-diffab-maaslin2:latest`), but they are private, so any downstream container build falls back to source.
- **Label:** Documented.
- **Reach:** High. No PyPI install, no public image, no lockable version.

---

## 6. What works well — do not break these

The project is a heavy, largely satisfied consumer. The 2026-07-12 refactor spec was explicitly *aggressive* — it set out to replace bespoke code with microsuite commands wherever one existed, accepting changed figures — and most of it succeeded. What it adopted is what works.

- **`microsuite trim` (cutadapt, fastp) and `microsuite qc` (fastqc, multiqc).** Left untouched by the aggressive refactor: *"**[keep]** already thin wrappers over `microsuite trim`/`qc` + external tools; no new command applies"* (`docs/superpowers/specs/2026-07-12-scripts-microsuite-refactor-design.md:122-124`). Used across all ~21 accessions.
- **`microsuite denoise --backend dada2-r`.** Complaint items 1-6 are all marked MITIGATED (`docs/microsuite_complaints.txt:8-45`): paired-FASTQ suffix parsing, sample-ID preservation in ASV tables, diagnostic plots, single-sample robustness, and the combined `microsuite-dada2` image. The follow-on provenance/QC work also landed — `README.md:221-225` confirms `dada2_denoise_manifest.json`, `dada2_qc_summary.tsv`, and the `--amplicon-length` overlap check now run under both local and docker runtimes. This is microsuite's best-received subsystem; the complaint file's own summary (`:154-157`) says so.
- **`microsuite tax_classify --backend qiime2`** with SILVA classifiers — kept as the taxonomy path (`refactor-design.md:66`), driving every accession.
- **The `feature_table.h5ad` hub: `import tsv` + `abundance` + `table normalize --method clr`.** `README.md:198-212`. This resolved complaints 12-15 outright, including the AnnData rank-column-name conflict that used to force a manual `feature_id` rewrite (`refactor-design.md:41-45`). CLR is now one command per rank instead of four blocks.
- **`microsuite tax_assignment_summary`.** Complaint #17 asked for exactly this (*"taxonomy_unassigned_summary.tsv is a useful concept and should be native"*, `microsuite_complaints.txt:126-130`); it shipped and is now invoked in four places.
- **`microsuite viz`** — `barplot`, `braycurtis-ordination` (with `--subject`/trajectory styling), `taxa-by-group`, `clr-by-group`. Adopted despite the figures visibly changing from the bespoke ones (`README.md:205-219`), which is a strong signal.
- **`microsuite diffab ancombc`** with `--fix-formula`/`--rand-formula`/`--reference` — adopted as the second differential-abundance method and wired with full option pass-through (`scripts/run_differential_abundance.sh:328-340`). This is the model the R backends should be brought up to.
- **The alpha and beta numerics themselves.** Across 29 docs and a dedicated complaints file, no one ever disputes a computed value. Every diversity finding in this report is about naming, packaging, output shape, or coverage — never arithmetic.
- **`microsuite phylogeny`, `microsuite functional_profile`, `microsuite diversity beta-significance`** — used as-is.

**Already resolved since the complaint file was written** (recorded so they are not re-litigated): per-backend R images now exist, superseding `scripts/containers/Dockerfile.maaslin2`'s stated reason for being (*"Keep this separate from microsuite's all-backend r-diffab image so optional ANCOM-BC/ALDEx2/LEfSe package availability cannot block this analysis"*, `:1-3`, after the incident at `docs/oral_microbiome_work_report_2026-07-10.md:437-438`); DADA2 parameter provenance and low-retention warnings (complaints 7-9, 11); TSV→CLR ergonomics and h5ad→TSV export (complaints 12-13).

---

## 7. Ranked fix list

Ordered by (impact × evidence strength) ÷ effort.

| # | Fix | Cat | Impact | Evidence | Effort |
|---|---|---|---|---|---|
| 1 | maaslin2 `normalization`/`transform` — default TSS, expose both | Correctness | High | Strong | Trivial |
| 2 | Document + disambiguate Simpson; add `gini_simpson` alias | Correctness | High | Strong | Trivial |
| 3 | `diversity alpha` accepts repeated `--metric`, emits one tidy table with `total_reads` | Ergonomics | High | Strong | Small |
| 4 | Multi-fixed-effect + random-effect support on the R diffab backends | Missing | Very high | Strong | Medium |
| 5 | `--run-dir` provenance on the native commands | Reproducibility | Very high | Strong | Medium |
| 6 | Beta metrics: square-matrix Sørensen, Aitchison, Newick-based UniFrac | Missing | High | Strong | Medium |
| 7 | Metadata alignment inside `import tsv` (ID column, missing-sample errors, reserved names) | Ergonomics | High | Medium | Small |
| 8 | Error (or flag) on zero-depth samples in alpha/normalize | Correctness | Medium | Medium | Trivial |
| 9 | Cross-run feature-table merge command | Missing | High | Strong | Medium |
| 10 | Batch as a modelled covariate; restricted permutations for beta tests | Missing | High | Medium | Large |
| 11 | PCA / t-SNE ordination | Missing | Medium | Strong | Small |
| 12 | Alpha-by-group plots + distance heatmap | Missing | Medium | Strong | Small |
| 13 | Normalize ANCOM-BC2 output schema; pin ANCOMBC | Ergonomics | Medium | Strong | Small |
| 14 | Multi-primer `--front` (accept a FASTA or repeated flag), documented | Ergonomics | Medium | Strong | Small |
| 15 | Publish to PyPI + public images with real version gating | Reproducibility | High | Strong | Large |

### Reasoning for the top five

**1 — maaslin2 normalization.** A one-line default change (`"NONE"` → `"TSS"`) plus two pass-through parameters closes the only finding in this audit where microsuite silently produces statistically invalid numbers under its documented, advertised usage. The consumer already demonstrated the correct setting, so there is no design question to resolve — only a decision about whether to treat the change as a bug fix or a breaking default. It should be a bug fix, with a note in the changelog: results produced by the old default are not comparable.

**2 — Simpson disambiguation.** Also near-zero effort, and it addresses the defect class the brief singles out: complete, well-formed, silently wrong. The two quantities are monotonically inverted, so the failure does not degrade a result, it reverses it. Adding a `gini_simpson` alias, documenting the convention in `docs/methods.md`, and emitting a one-time note when bare `simpson` is requested costs an afternoon. That the only consumer we have evidence from felt compelled to maintain a private remap table is the argument.

**3 — Alpha multi-metric output.** This single change deletes `finalize_microsuite_diversity.py` entirely — the eight-subprocess loop, the eight-file merge, the column renaming, and the recomputation of sample depth from the count table. It is the highest ratio of consumer code eliminated to microsuite code written in the whole list. Bundle it with fix 2, since the column-naming decision is the same decision, and ship `total_reads` in the same table while you are there.

**4 — Random and multiple fixed effects.** The most-cited gap in the repo, stated verbatim in a work report, and the reason for an entire parallel R + Docker + six-Python-helper stack. Effort is genuinely medium — it means threading a formula through `diff_abundance`, `r_backends.py`, and the three R scripts — but the design is already settled: `cli/diffab_cmd.py:25` shows exactly what the option surface should look like, and the consumer's `run_maaslin2_repeated_measures.R` is a working reference implementation of the MaAsLin2 half, rank-deficiency guard included. Copy the rank-deficiency guard too; they needed it on their first real model.

**5 — Native-command provenance.** Ranked fifth rather than higher only because the effort is real and the payoff is invisible until something goes wrong. But note what the absence caused: a sophisticated consumer built a complete parallel provenance system (`workflow.json`, `run_manifest.json`, `stage-result.v1`, a registry, 319 lines of reporting code) and then wrote a compatibility shim pinning a dev-branch SHA to get at the half-built microsuite version. The metadata module already exists (`src/microsuite/metadata/`, twelve files); the gap is that the native CLI commands do not call it. Wiring `--run-dir` through `diversity_cmd`, `table_cmd`, `viz_cmd`, and `import_cmd` is mostly mechanical and would let the consumer delete their shim and adopt the real contract. Doing this without also doing fix 15 (real versions) leaves them still pinning SHAs, so consider sequencing those together.
