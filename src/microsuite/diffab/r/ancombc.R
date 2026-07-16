args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: ancombc.R counts.tsv metadata.tsv params.json output.tsv")
}
counts_path <- args[[1]]
metadata_path <- args[[2]]
params_path <- args[[3]]
output_path <- args[[4]]

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("R package 'jsonlite' is required to read the ANCOM-BC parameters.")
}
params <- jsonlite::fromJSON(params_path)

counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_path, row.names = 1, check.names = FALSE)

fix_formula <- params$fix_formula
rand_formula <- if (is.null(params$rand_formula)) NULL else params$rand_formula
group_col <- if (is.null(params$group)) NULL else params$group

# Apply reference levels (base R; before the ANCOMBC requirement so validation is
# reachable without the heavy package installed).
reference <- params$reference
if (length(reference)) {
  for (col in names(reference)) {
    if (!(col %in% colnames(metadata))) stop(sprintf("Reference column not found: %s", col))
    lvl <- reference[[col]]
    metadata[[col]] <- factor(metadata[[col]])
    if (!(lvl %in% levels(metadata[[col]]))) {
      stop(sprintf("Reference level '%s' not found in column '%s'.", lvl, col))
    }
    metadata[[col]] <- stats::relevel(metadata[[col]], ref = lvl)
  }
}

check_formula_columns <- function(fml) {
  if (is.null(fml) || !nzchar(fml)) return(invisible(NULL))
  vars <- all.vars(stats::as.formula(paste("~", fml)))
  missing <- setdiff(vars, colnames(metadata))
  if (length(missing)) {
    stop(sprintf("Formula references unknown metadata columns: %s", paste(missing, collapse = ", ")))
  }
}
check_formula_columns(fix_formula)
check_formula_columns(rand_formula)

# Full-rank fixed-effects design (rejects confounded models, e.g. phase + time).
mm <- stats::model.matrix(stats::as.formula(paste("~", fix_formula)), data = metadata)
rank <- qr(mm)$rank
if (rank < ncol(mm)) {
  stop(sprintf(
    "Fixed-effects design is rank deficient (rank %d < %d columns): the model is confounded. Simplify fix_formula or drop a collinear term.",
    rank, ncol(mm)
  ))
}

# Group-dependent tests need a group column.
group_tests <- c(params$global, params$pairwise, params$trend, params$dunnet)
if (any(unlist(group_tests)) && is.null(group_col)) {
  stop("--global/--pairwise/--trend/--dunnet require --group.")
}

tryCatch(
  loadNamespace("ANCOMBC"),
  error = function(exc) {
    stop(sprintf("R package 'ANCOMBC' could not be loaded: %s", conditionMessage(exc)))
  }
)
if (!("ancombc2" %in% getNamespaceExports("ANCOMBC"))) {
  stop("This command requires ANCOM-BC2; update the ANCOMBC package.")
}

fit <- ANCOMBC::ancombc2(
  data = counts,
  meta_data = metadata,
  fix_formula = fix_formula,
  rand_formula = rand_formula,
  group = group_col,
  p_adj_method = params$p_adj_method,
  prv_cut = params$prv_cut,
  lib_cut = params$lib_cut,
  struc_zero = isTRUE(params$struc_zero),
  neg_lb = isTRUE(params$neg_lb),
  global = isTRUE(params$global),
  pairwise = isTRUE(params$pairwise),
  trend = isTRUE(params$trend),
  dunnet = isTRUE(params$dunnet),
  pseudo_sens = isTRUE(params$pseudo_sens),
  n_cl = params$n_cl
)
write.table(fit$res, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)

# Resolved-config provenance beside the output.
model_factors <- Filter(is.factor, metadata)
provenance <- list(
  fix_formula = fix_formula,
  rand_formula = if (is.null(rand_formula)) NA else rand_formula,
  group = if (is.null(group_col)) NA else group_col,
  reference = reference,
  controls = list(
    p_adj_method = params$p_adj_method, prv_cut = params$prv_cut, lib_cut = params$lib_cut,
    struc_zero = isTRUE(params$struc_zero), neg_lb = isTRUE(params$neg_lb),
    global = isTRUE(params$global), pairwise = isTRUE(params$pairwise),
    trend = isTRUE(params$trend), dunnet = isTRUE(params$dunnet),
    pseudo_sens = isTRUE(params$pseudo_sens), n_cl = params$n_cl
  ),
  factor_reference_levels = lapply(model_factors, function(x) levels(x)[1]),
  ancombc_version = as.character(utils::packageVersion("ANCOMBC")),
  r_version = R.version.string
)
prov_path <- file.path(dirname(output_path), "ancombc_provenance.json")
writeLines(jsonlite::toJSON(provenance, auto_unbox = TRUE, null = "null", pretty = TRUE), prov_path)
