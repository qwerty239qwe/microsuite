args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: lefse.R counts.tsv metadata.tsv group_col output.tsv")
}

counts_path <- args[[1]]
metadata_path <- args[[2]]
group_col <- args[[3]]
output_path <- args[[4]]

if (!requireNamespace("lefser", quietly = TRUE)) {
  stop("R package 'lefser' is required. Install it with BiocManager::install('lefser').")
}
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  stop("R package 'SummarizedExperiment' is required.")
}

counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_path, row.names = 1, check.names = FALSE)
if (!(group_col %in% colnames(metadata))) {
  stop(paste("Group column not found:", group_col))
}

metadata <- metadata[colnames(counts), , drop = FALSE]
conditions <- factor(metadata[[group_col]])
if (nlevels(conditions) != 2) {
  stop("LEfSe requires exactly two groups.")
}

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = as.matrix(counts)),
  colData = metadata
)
# LEfSe expects relative abundances, not raw counts. Given counts, lefser only
# emits a warning and carries on, producing LDA scores confounded by sequencing
# depth -- a complete, plausible, wrong result. relativeAb converts to counts
# per million, which is what the method assumes.
se <- lefser::relativeAb(se)
# lefser renamed this parameter across releases: `groupCol` in the versions
# bioconda ships for R 4.3, `classCol` in current devel (1.23.0). Passing the
# wrong one is not a hard error -- it lands in `...` while the real parameter
# takes its default -- so lefser fails with "must refer to a valid dichotomous
# (two-level) variable" even though the column is a clean two-level factor.
# Ask the installed function which name it takes rather than pinning a guess.
lefser_formals <- names(formals(lefser::lefser))
class_arg <- if ("classCol" %in% lefser_formals) {
  "classCol"
} else if ("groupCol" %in% lefser_formals) {
  "groupCol"
} else {
  stop(
    "Unsupported lefser version: lefser() takes neither 'classCol' nor ",
    "'groupCol'. Found: ", paste(lefser_formals, collapse = ", ")
  )
}
lefser_args <- list(se)
lefser_args[[class_arg]] <- group_col
result <- do.call(lefser::lefser, lefser_args)
write.table(result, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
