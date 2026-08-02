#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
    stop("usage: capture_reference.R /path/to/SpiecEasi [fixture-output-dir]")
}

expected_commit <- "faed6a4476fe0a8dc701ea15cbdfe98d56ce6704"
source_dir <- normalizePath(args[[1L]], mustWork = TRUE)
fixture_dir <- if (length(args) == 2L) {
    normalizePath(args[[2L]], mustWork = TRUE)
} else {
    normalizePath(dirname(sys.frame(1)$ofile), mustWork = TRUE)
}

actual_commit <- system2(
    "git",
    c("-C", shQuote(source_dir), "rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = TRUE
)
if (length(actual_commit) != 1L || !identical(actual_commit, expected_commit)) {
    stop(
        sprintf(
            "SpiecEasi checkout must be exactly %s; found %s",
            expected_commit,
            paste(actual_commit, collapse = " ")
        )
    )
}

description <- read.dcf(file.path(source_dir, "DESCRIPTION"))
if (!identical(unname(description[1L, "Version"]), "1.99.0")) {
    stop("unexpected SpiecEasi DESCRIPTION version")
}
if (!requireNamespace("VGAM", quietly = TRUE)) {
    stop("VGAM is required to capture SpiecEasi reference matrices")
}
if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("MASS is required for the SpiecEasi generalized-inverse fallback")
}

# Keep the oracle boundary explicit: load only normalization, the file that
# provides cor2cov(), and the SparCC implementation from the pinned checkout.
source(file.path(source_dir, "R", "normalization.R"), local = globalenv())
source(file.path(source_dir, "R", "mvdistributions.R"), local = globalenv())
source(file.path(source_dir, "R", "spaRcc.R"), local = globalenv())

options(digits = 17, scipen = 999)
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

read_fixture <- function(filename) {
    value <- as.matrix(
        read.delim(
            file.path(fixture_dir, filename),
            row.names = 1L,
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
    )
    storage.mode(value) <- "double"
    value
}

write_matrix <- function(value, filename, feature_names) {
    if (!all(is.finite(value))) stop(sprintf("non-finite output for %s", filename))
    dimnames(value) <- list(feature_names, feature_names)
    write.table(
        value,
        file = file.path(fixture_dir, filename),
        sep = "\t",
        quote = FALSE,
        row.names = TRUE,
        col.names = NA
    )
}

seeds <- c(10010L, 10011L, 10012L)
for (dataset in c("dense", "zero")) {
    counts <- read_fixture(sprintf("%s_counts.tsv", dataset))
    for (seed in seeds) {
        set.seed(seed)
        result <- sparcc(counts)
        write_matrix(
            result$Cor,
            sprintf("%s_reference_cor_seed_%d.tsv", dataset, seed),
            colnames(counts)
        )
    }
}

inner <- read_fixture("inner_compositions.tsv")
inner_variation <- av(inner)
inner_basis <- basis_var(inner_variation)
inner_initial <- C_from_V(inner_variation, inner_basis$Vbase)$Cor
inner_completed <- sparccinner(inner)$Cor
write_matrix(
    inner_initial,
    "inner_initial_reference_cor.tsv",
    colnames(inner)
)
write_matrix(
    inner_completed,
    "inner_reference_cor.tsv",
    colnames(inner)
)

cat(sprintf("SpiecEasi commit: %s\n", actual_commit))
cat(sprintf("SpiecEasi version: %s\n", description[1L, "Version"]))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("VGAM version: %s\n", as.character(utils::packageVersion("VGAM"))))
cat(sprintf("MASS version: %s\n", as.character(utils::packageVersion("MASS"))))
