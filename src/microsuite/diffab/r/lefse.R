args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: lefse.R counts.tsv metadata.tsv params.json output.tsv")
}

counts_path <- args[[1]]
metadata_path <- args[[2]]
params_path <- args[[3]]
output_path <- args[[4]]

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("R package 'jsonlite' is required.")
}
if (!requireNamespace("lefser", quietly = TRUE)) {
  stop("R package 'lefser' is required. Install it with BiocManager::install('lefser').")
}
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  stop("R package 'SummarizedExperiment' is required.")
}

params <- jsonlite::fromJSON(params_path, simplifyVector = TRUE)
required <- c(
  "group", "reference", "comparison", "seed", "kruskal_threshold",
  "wilcoxon_threshold", "lda_threshold", "p_adjust_method", "trim_names"
)
missing_params <- setdiff(required, names(params))
if (length(missing_params) > 0) {
  stop(paste("Missing LEfSe parameters:", paste(missing_params, collapse = ", ")))
}

counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_path, row.names = 1, check.names = FALSE)
group_col <- params$group
if (!(group_col %in% colnames(metadata))) {
  stop(paste("Group column not found:", group_col))
}

missing_metadata <- setdiff(colnames(counts), rownames(metadata))
if (length(missing_metadata) > 0) {
  stop(paste("Samples missing from metadata:", paste(missing_metadata, collapse = ", ")))
}
metadata <- metadata[colnames(counts), , drop = FALSE]
metadata[[group_col]] <- factor(
  as.character(metadata[[group_col]]),
  levels = c(params$reference, params$comparison)
)
if (anyNA(metadata[[group_col]]) || nlevels(metadata[[group_col]]) != 2) {
  stop("LEfSe requires exactly two groups.")
}
if (!is.null(params$subclass)) {
  subclass_col <- params$subclass
  if (!(subclass_col %in% colnames(metadata))) {
    stop(paste("Subclass column not found:", subclass_col))
  }
  metadata[[subclass_col]] <- factor(as.character(metadata[[subclass_col]]))
  combinations <- table(metadata[[subclass_col]], metadata[[group_col]])
  if (ncol(combinations) != 2 || any(combinations == 0)) {
    stop(
      "LEfSe subclass levels must be represented in both groups; subclass is a ",
      "crossed blocking/replicate factor, not a nested or random-effect term."
    )
  }
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
lefser_args$kruskal.threshold <- params$kruskal_threshold
lefser_args$wilcox.threshold <- params$wilcoxon_threshold
lefser_args$lda.threshold <- params$lda_threshold
lefser_args$trim.names <- params$trim_names
if ("checkAbundances" %in% lefser_formals) {
  lefser_args$checkAbundances <- TRUE
}
if ("method" %in% lefser_formals) {
  lefser_args$method <- params$p_adjust_method
} else if (!identical(params$p_adjust_method, "none")) {
  stop("Installed lefser does not support p-value adjustment via the 'method' argument.")
}
if (!is.null(params$subclass)) {
  subclass_arg <- if ("subclassCol" %in% lefser_formals) {
    "subclassCol"
  } else if ("blockCol" %in% lefser_formals) {
    "blockCol"
  } else {
    stop("Installed lefser does not support a subclass/block column.")
  }
  lefser_args[[subclass_arg]] <- params$subclass
}
set.seed(as.integer(params$seed))
if (is.null(params$subclass)) {
  result <- do.call(lefser::lefser, lefser_args)
} else {
  # lefser 1.22.0's fillPmatZmat() builds per-comparison sample indices with
  # apply(). For a balanced class x subclass table, apply() simplifies the
  # equal-length vectors to a matrix; seq_along() then sends one sample at a
  # time to coin::wilcox_test(), which fails with an IndependenceProblem error.
  # Keep the upstream algorithm but force the indices to remain a list.
  internal_names <- c("filterKruskal", "wilcox_pstats", "ldaFunction", ".trunc")
  internals <- lapply(internal_names, function(name) {
    tryCatch(getFromNamespace(name, "lefser"), error = function(error) NULL)
  })
  names(internals) <- internal_names
  if (any(vapply(internals, is.null, logical(1)))) {
    stop(
      "Installed lefser does not expose the internals required for reliable ",
      "subclass analysis. Use the pinned r-diffab-lefse image."
    )
  }

  relab_data <- SummarizedExperiment::assay(se, i = 1L)
  class_factor <- as.factor(SummarizedExperiment::colData(se)[[group_col]])
  subclass_factor <- droplevels(
    as.factor(SummarizedExperiment::colData(se)[[params$subclass]])
  )
  relab_sub <- internals$filterKruskal(
    relab = relab_data,
    class = class_factor,
    p.value = params$kruskal_threshold,
    method = params$p_adjust_method
  )

  if (nrow(relab_sub) > 0L) {
    subclass_indices <- seq_along(levels(subclass_factor))
    comparisons <- expand.grid(
      subclass_indices,
      subclass_indices
    )
    indices <- lapply(seq_len(nrow(comparisons)), function(index) {
      first_subclass <- levels(subclass_factor)[comparisons[index, 1L]]
      second_subclass <- levels(subclass_factor)[comparisons[index, 2L]]
      c(
        which(
          class_factor == levels(class_factor)[1L] &
            subclass_factor == first_subclass
        ),
        which(
          class_factor == levels(class_factor)[2L] &
            subclass_factor == second_subclass
        )
      )
    })
    z_matrix <- pvalue_matrix <- matrix(
      NA_real_,
      nrow = nrow(relab_sub),
      ncol = length(indices),
      dimnames = list(rownames(relab_sub), NULL)
    )
    for (index in seq_along(indices)) {
      tests <- internals$wilcox_pstats(
        relab_sub,
        class = class_factor,
        index = indices[[index]]
      )
      pvalue_matrix[, index] <- stats::p.adjust(
        tests["pvalue", ], method = params$p_adjust_method
      )
      z_matrix[, index] <- tests["statistic", ]
    }
    passes_pvalue <- !is.na(pvalue_matrix) &
      pvalue_matrix <= params$wilcoxon_threshold * 2
    keep <- rowSums(passes_pvalue) == ncol(passes_pvalue)
    z_kept <- z_matrix[keep, , drop = FALSE]
    consistent_direction <- abs(rowSums(z_kept)) == rowSums(abs(z_kept))
    relab_sub <- relab_sub[names(consistent_direction[consistent_direction]), , drop = FALSE]
  }

  if (nrow(relab_sub) == 0L) {
    result <- data.frame(features = character(), scores = numeric())
  } else {
    lda_input <- as.data.frame(t(relab_sub))
    lda_input$class <- class_factor
    raw_scores <- internals$ldaFunction(lda_input, levels(class_factor))
    scores <- sign(raw_scores) * log10(1 + abs(raw_scores))
    result <- data.frame(
      features = names(scores),
      scores = as.vector(scores),
      stringsAsFactors = FALSE
    )
    result <- internals$.trunc(result, params$trim_names)
    result <- result[abs(result$scores) >= params$lda_threshold, , drop = FALSE]
  }
}
if (!is.data.frame(result) || ncol(result) != 2) {
  stop("lefser returned an unexpected result schema; expected a two-column data.frame.")
}
colnames(result) <- c("features", "scores")
result$features <- as.character(result$features)
result$scores <- as.numeric(result$scores)
result <- result[order(-abs(result$scores), result$features, method = "radix"), , drop = FALSE]
write.table(result, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
