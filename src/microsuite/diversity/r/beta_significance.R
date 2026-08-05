#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
    stop(
        paste(
            "Usage:",
            "beta_significance.R DISTANCE.tsv METADATA.tsv METHOD FORMULA STRATA PERMUTATIONS SEED OUTPUT.tsv"
        )
    )
}

suppressPackageStartupMessages(library(vegan))

distance_path <- args[[1]]
metadata_path <- args[[2]]
method <- tolower(args[[3]])
formula_text <- trimws(args[[4]])
strata_name <- trimws(args[[5]])
permutations <- as.integer(args[[6]])
seed <- as.integer(args[[7]])
output_path <- args[[8]]

if (!method %in% c("adonis2", "anosim2")) {
    stop("method must be adonis2 or anosim2")
}
if (!nzchar(formula_text)) {
    stop("formula must not be empty")
}
if (is.na(permutations) || permutations < 0L) {
    stop("permutations must be a non-negative integer")
}
if (is.na(seed)) {
    stop("seed must be an integer")
}

read_matrix <- function(path) {
    table <- read.delim(
        path,
        header = TRUE,
        row.names = 1,
        check.names = FALSE,
        comment.char = "",
        quote = ""
    )
    matrix <- as.matrix(table)
    storage.mode(matrix) <- "double"
    matrix
}

distance_matrix <- read_matrix(distance_path)
metadata <- read.delim(
    metadata_path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    comment.char = "",
    quote = ""
)

if (is.null(rownames(distance_matrix)) || is.null(colnames(distance_matrix))) {
    stop("distance matrix must have sample IDs in both dimensions")
}
if (!all(rownames(distance_matrix) %in% rownames(metadata))) {
    stop("distance matrix contains sample IDs missing from metadata")
}
metadata <- metadata[rownames(distance_matrix), , drop = FALSE]
distance <- stats::as.dist(distance_matrix)
strata_values <- NULL
if (nzchar(strata_name)) {
    if (!strata_name %in% colnames(metadata)) {
        stop(sprintf("strata column not found: %s", strata_name))
    }
    strata_values <- metadata[[strata_name]]
}

set.seed(seed)
permutation_scheme <- if (permutations == 0L) {
    "none"
} else if (is.null(strata_values)) {
    "unrestricted"
} else {
    "blocked"
}

if (method == "adonis2") {
    model_formula <- stats::as.formula(paste("distance ~", formula_text))
    fit <- vegan::adonis2(
        model_formula,
        data = metadata,
        permutations = permutations,
        strata = strata_values,
        by = "terms"
    )
    fit <- as.data.frame(fit)
    model_rows <- !tolower(rownames(fit)) %in% c("residual", "total")
    fit <- fit[model_rows, , drop = FALSE]
    get_column <- function(name) {
        if (name %in% colnames(fit)) fit[[name]] else rep(NA_real_, nrow(fit))
    }
    output <- data.frame(
        backend = "vegan",
        method = "adonis2",
        term = rownames(fit),
        df = get_column("Df"),
        sum_of_squares = get_column("SumOfSqs"),
        r_squared = get_column("R2"),
        f_value = get_column("F"),
        p_value = get_column("Pr(>F)"),
        formula = formula_text,
        strata = if (nzchar(strata_name)) strata_name else NA_character_,
        permutations = permutations,
        permutation_scheme = permutation_scheme,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
} else {
    if (!formula_text %in% colnames(metadata)) {
        stop("anosim2 formula must be one metadata grouping column")
    }
    grouping <- metadata[[formula_text]]
    if (length(unique(grouping)) < 2L) {
        stop("ANOSIM requires at least two groups")
    }
    fit <- vegan::anosim(
        distance,
        grouping = grouping,
        permutations = permutations,
        strata = strata_values
    )
    statistic <- unname(fit$statistic)
    p_value <- unname(fit$signif)
    output <- data.frame(
        backend = "vegan",
        method = "anosim2",
        term = formula_text,
        statistic = statistic,
        r = statistic,
        p_value = p_value,
        n_groups = length(unique(grouping)),
        formula = formula_text,
        strata = if (nzchar(strata_name)) strata_name else NA_character_,
        permutations = permutations,
        permutation_scheme = permutation_scheme,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

write.table(
    output,
    file = output_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)
