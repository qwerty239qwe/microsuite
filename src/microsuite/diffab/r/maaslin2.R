args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: maaslin2.R counts.tsv metadata.tsv group_col output.tsv")
}

counts_path <- args[[1]]
metadata_path <- args[[2]]
group_col <- args[[3]]
output_path <- args[[4]]

if (!requireNamespace("Maaslin2", quietly = TRUE)) {
  stop("R package 'Maaslin2' is required. Install it with BiocManager::install('Maaslin2').")
}

counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_path, row.names = 1, check.names = FALSE)
if (!(group_col %in% colnames(metadata))) {
  stop(paste("Group column not found:", group_col))
}

metadata <- metadata[colnames(counts), , drop = FALSE]
features <- as.data.frame(t(counts), check.names = FALSE)

run_dir <- tempfile("maaslin2-")
dir.create(run_dir)
# normalization must not be "NONE" here. The caller passes a raw count table,
# so without library-size normalization every feature in a deeply sequenced
# sample carries a systematically larger value, and the model attributes that
# depth difference to whichever covariate happens to correlate with it. The
# output stays well formed and the p-values look plausible, so the error is
# invisible. TSS is MaAsLin 2's own default and pairs with transform = "LOG".
Maaslin2::Maaslin2(
  input_data = features,
  input_metadata = metadata,
  output = run_dir,
  fixed_effects = c(group_col),
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  standardize = FALSE,
  plot_heatmap = FALSE,
  plot_scatter = FALSE,
  cores = 1
)

results_path <- file.path(run_dir, "all_results.tsv")
if (!file.exists(results_path)) {
  stop("MaAsLin2 did not produce all_results.tsv")
}
results <- read.delim(results_path, check.names = FALSE)
write.table(results, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
