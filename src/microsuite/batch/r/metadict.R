#!/usr/bin/env Rscript
# MetaDICT: covariate balancing plus shared dictionary learning, returning
# relative abundances (per the capability record in backends.py).
#
# UNVERIFIED AT RUNTIME: this signature is taken from the published man page
# (man/metadict.Rd at the pinned commit,
# https://github.com/BoYuan07/MetaDICT) and, because that man page alone
# left open questions (see below), from the package source itself
# (R/MetaDICT.R and R/data_check.R at the same commit) -- there is no
# container engine available in the environment this script was written in,
# so it has never been executed against the real MetaDICT package. The
# container's build-time smoke test (see
# containers/r-batch-metadict/Dockerfile) is this script's first real
# execution. If that smoke fails, re-check this script against
# `Rscript -e "args(MetaDICT::MetaDICT)"` inside the built image.
#
# man/metadict.Rd / R/MetaDICT.R confirm:
#   MetaDICT(count, meta, covariates = "all", tree = NULL, taxonomy = NULL,
#            distance_matrix = NULL, tax_level = NULL,
#            customize_parameter = FALSE, alpha = 0.1, beta = 0.01,
#            normalization = "uq", max_iter = 10000, imputation = FALSE,
#            verbose = TRUE, optim_trace = FALSE)
# - `count` is documented as "taxa-by-sample" -- features-as-rows, the same
#   orientation microsuite writes counts.tsv in, so (unlike conqur.R /
#   plsda_batch.R) NO transpose is needed on the way in.
# - `meta` "must include a column named 'batch'" -- a literal column name,
#   not a name passed as an argument. There is no `batch=` parameter.
# - `covariates` is a character vector of column names in `meta`, or the
#   string "all" (its default) meaning every non-"batch" column.
#   R/data_check.R: `if (length(covariates)==1 && covariates=="all") ...
#   else if (is.vector(covariates)) meta_filtered <- meta[, covariates, ...]`
#   -- a zero-length vector is a valid `is.vector` value and selects zero
#   columns, so "no covariates requested" is expressed natively as
#   `character(0)`, unlike ConQuR, which has no such native empty case.
# - Return value: `list(count = <data.frame, taxa x samples>, D=..., R=...,
#   w=..., meta=..., dist_mat=...)`. `$count` matches the plan's assumed
#   `$X` in spirit (a features x samples corrected table) but not in name.
#
# ONE MORE THING R/data_check.R REVEALS THAT THE MAN PAGE DOES NOT SAY:
# MetaDICT *requires* one of `distance_matrix`, `tree`, or `taxonomy` --
# with all three NULL it calls `stop("Neither a phylogenetic tree nor a
# taxonomy table is provided...")`. None of those three is expressible in
# microsuite's four-argument contract (counts.tsv, metadata.tsv,
# params.json, corrected.tsv) as currently defined by any other backend.
# This script therefore defines its OWN, UNDOCUMENTED-BY-METADICT
# extension of the params.json contract: an optional
# `distance_matrix_path`, a TSV of taxon-by-taxon dissimilarities with row
# and column names matching counts.tsv's feature ids. When absent (as in
# the build-time smoke below, which supplies no taxonomy/tree/distance
# data), this script falls back to a uniform distance matrix (every
# off-diagonal entry equal, diagonal zero) -- i.e., "no information about
# which taxa are similar," the most neutral assumption available. This
# fallback is an INFERENCE, not something metadict.Rd or data_check.R
# recommends; MetaDICT's own examples always supply a real dist_mat.
#
# A SECOND UNDOCUMENTED BEHAVIOUR FOUND ONLY BY READING R/MetaDICT.R:
# the function builds its corrected matrix `X` by cbind-ing per-batch
# blocks in the order `unique(meta$batch)` enumerates them, then does
# `colnames(res.metadict) <- colnames(count)` -- reassigning the ORIGINAL
# input's column names onto a matrix whose actual column order is grouped
# by batch. If a caller's samples are not already contiguous by batch,
# this silently mislabels columns (right names, wrong data underneath).
# To neutralize this without patching the vendored package, this script
# sorts `count`/`meta` by batch before calling MetaDICT, which makes the
# per-batch block order equal to the (now-sorted) input order, so the
# reassigned colnames are correct by construction; the output is then
# restored to the caller's original sample order before writing.
#
# Usage: metadict.R counts.tsv metadata.tsv params.json corrected.tsv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: metadict.R counts.tsv metadata.tsv params.json corrected.tsv")
}
suppressPackageStartupMessages({
  library(jsonlite)
  library(MetaDICT)
})

counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
meta <- read.delim(args[2], row.names = 1, check.names = FALSE)
params <- fromJSON(args[3], simplifyVector = TRUE)

meta <- meta[colnames(counts), , drop = FALSE]

# MetaDICT reads the batch id from a column literally named "batch" (see
# above), not from a name passed as an argument, so the caller's chosen
# batch column is copied into one named exactly that.
meta[["batch"]] <- factor(meta[[params$batch]])

sample_order <- colnames(counts)

# Sort samples into contiguous batch blocks before calling MetaDICT, to
# neutralize the column-mislabeling behaviour documented above in
# R/MetaDICT.R. Restored to `sample_order` before writing.
ord <- order(meta[["batch"]])
counts_sorted <- counts[, ord, drop = FALSE]
meta_sorted <- meta[ord, , drop = FALSE]

covariates <- if (length(params$covariates) > 0) as.character(params$covariates) else character(0)

# OUR OWN CONTRACT EXTENSION (see header comment): an optional taxon
# dissimilarity matrix supplied via params.json, since none of
# distance_matrix/tree/taxonomy has a home in the four positional
# arguments. Falls back to an uninformative uniform matrix when absent.
if (!is.null(params$distance_matrix_path) && !is.na(params$distance_matrix_path)) {
  dist_raw <- as.matrix(read.delim(params$distance_matrix_path, row.names = 1, check.names = FALSE))
  dist_mat <- dist_raw[rownames(counts_sorted), rownames(counts_sorted)]
} else {
  # INFERRED, UNDOCUMENTED FALLBACK: MetaDICT (R/data_check.R) requires a
  # square taxa-by-taxa distance matrix with zero on the diagonal when
  # `distance_matrix` is supplied directly (as opposed to being derived
  # from `tree`/`taxonomy`); nothing in the docs says what a "no
  # information available" distance matrix should look like. A uniform
  # off-diagonal value treats every pair of taxa as equally (dis)similar,
  # which is this script's least-committal choice, not a package default.
  n_features <- nrow(counts_sorted)
  dist_mat <- matrix(1, nrow = n_features, ncol = n_features)
  diag(dist_mat) <- 0
  rownames(dist_mat) <- rownames(counts_sorted)
  colnames(dist_mat) <- rownames(counts_sorted)
}

fit <- MetaDICT(
  count = counts_sorted,
  meta = meta_sorted,
  covariates = covariates,
  distance_matrix = dist_mat
)
adjusted <- as.matrix(fit$count)

# Restore the caller's original sample order (see header comment on the
# batch-sort workaround), then features-as-rows with feature_id first.
adjusted <- adjusted[, sample_order, drop = FALSE]
out <- data.frame(feature_id = rownames(adjusted), adjusted, check.names = FALSE)
write.table(out, args[4], sep = "\t", quote = FALSE, row.names = FALSE)
