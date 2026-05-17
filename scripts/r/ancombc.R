args <- commandArgs(trailingOnly = TRUE)
counts_path <- args[[1]]
metadata_path <- args[[2]]
group_col <- args[[3]]
output_path <- args[[4]]

if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
  stop("R package 'ANCOMBC' is required. Install it with BiocManager::install('ANCOMBC').")
}

counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_path, row.names = 1, check.names = FALSE)
if (!(group_col %in% colnames(metadata))) {
  stop(paste("Group column not found:", group_col))
}

if ("ancombc2" %in% getNamespaceExports("ANCOMBC")) {
  fit <- ANCOMBC::ancombc2(
    data = counts,
    meta_data = metadata,
    fix_formula = group_col,
    group = group_col,
    p_adj_method = "BH",
    prv_cut = 0,
    lib_cut = 0,
    struc_zero = FALSE,
    neg_lb = FALSE,
    global = FALSE
  )
  result <- fit$res
} else {
  fit <- ANCOMBC::ancombc(
    phyloseq = NULL,
    feature_table = counts,
    meta_data = metadata,
    formula = group_col,
    p_adj_method = "BH",
    prv_cut = 0,
    lib_cut = 0
  )
  result <- fit$res
}

write.table(result, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
