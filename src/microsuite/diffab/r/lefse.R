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
result <- lefser::lefser(se, classCol = group_col)
write.table(result, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
