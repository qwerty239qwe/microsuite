args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: aldex2.R counts.tsv metadata.tsv group_col output.tsv")
}

counts_path <- args[[1]]
metadata_path <- args[[2]]
group_col <- args[[3]]
output_path <- args[[4]]

if (!requireNamespace("ALDEx2", quietly = TRUE)) {
  stop("R package 'ALDEx2' is required. Install it with BiocManager::install('ALDEx2').")
}

counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_path, row.names = 1, check.names = FALSE)
if (!(group_col %in% colnames(metadata))) {
  stop(paste("Group column not found:", group_col))
}

metadata <- metadata[colnames(counts), , drop = FALSE]
conditions <- as.character(metadata[[group_col]])
if (length(unique(conditions)) < 2) {
  stop("ALDEx2 requires at least two groups.")
}

clr <- ALDEx2::aldex.clr(
  counts,
  conditions,
  mc.samples = 128,
  denom = "all",
  verbose = FALSE
)

if (length(unique(conditions)) == 2) {
  test <- ALDEx2::aldex.ttest(clr, paired.test = FALSE, verbose = FALSE)
  effect <- ALDEx2::aldex.effect(clr, verbose = FALSE)
  result <- cbind(feature = rownames(test), test, effect)
} else {
  test <- ALDEx2::aldex.kw(clr)
  result <- cbind(feature = rownames(test), test)
}

write.table(result, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
