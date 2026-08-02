#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
    stop(
        paste(
            "usage: run_spieceasi_reference.R",
            "/path/to/SpiecEasi /path/to/fixtures /path/to/output"
        )
    )
}

expected_commit <- "faed6a4476fe0a8dc701ea15cbdfe98d56ce6704"
source_dir <- normalizePath(args[[1L]], mustWork = TRUE)
fixture_dir <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = TRUE)

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
spieceasi_version <- unname(description[1L, "Version"])
if (!identical(spieceasi_version, "1.99.0")) {
    stop(sprintf("unexpected SpiecEasi DESCRIPTION version: %s", spieceasi_version))
}
if (!requireNamespace("VGAM", quietly = TRUE)) {
    stop("VGAM is required to run the live SpiecEasi reference")
}
if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("MASS is required for the SpiecEasi generalized-inverse fallback")
}

cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("VGAM version: %s\n", as.character(utils::packageVersion("VGAM"))))
cat(sprintf("SpiecEasi version: %s\n", spieceasi_version))
cat(sprintf("SpiecEasi commit: %s\n", actual_commit))

# Keep the GPL oracle outside this repository and load only the three files
# required to execute SparCC against the committed count fixtures.
source(file.path(source_dir, "R", "normalization.R"), local = globalenv())
source(file.path(source_dir, "R", "mvdistributions.R"), local = globalenv())
source(file.path(source_dir, "R", "spaRcc.R"), local = globalenv())

options(digits = 17, scipen = 999)
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

read_counts <- function(dataset) {
    value <- as.matrix(
        read.delim(
            file.path(fixture_dir, sprintf("%s_counts.tsv", dataset)),
            row.names = 1L,
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
    )
    storage.mode(value) <- "double"
    value
}

write_correlation <- function(value, dataset, feature_names) {
    if (!all(is.finite(value))) {
        stop(sprintf("non-finite SpiecEasi correlation for %s", dataset))
    }
    dimnames(value) <- list(feature_names, feature_names)
    write.table(
        value,
        file = file.path(output_dir, sprintf("%s_reference_cor.tsv", dataset)),
        sep = "\t",
        quote = FALSE,
        row.names = TRUE,
        col.names = NA
    )
}

for (dataset in c("dense", "zero")) {
    counts <- read_counts(dataset)
    set.seed(10010L)
    result <- sparcc(counts)
    write_correlation(result$Cor, dataset, colnames(counts))
}
