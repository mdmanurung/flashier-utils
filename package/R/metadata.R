# Independent descriptive summaries; ImmGenT and jamboree group-ranking workflow credit.
.metadata_activity <- function(fit, metadata, required, id_col, na_action) {
  view <- .resolve_view(fit)
  match <- .match_ids(rownames(view$activity), metadata, id_col, na_action = na_action,
                      required_columns = unique(c(id_col, required)))
  list(view = view, metadata = match$metadata, A = view$activity[match$matching$matrix_row, , drop = FALSE], matching = match$matching, exclusions = match$exclusions)
}

.weights <- function(weights, ids) {
  if (is.null(weights)) return(rep(1, length(ids)))
  if (!is.numeric(weights)) stop("weights must be numeric and named by sample")
  .validate_ids(names(weights), "Weight sample")
  if (any(!is.finite(weights)) || any(weights < 0)) stop("weights must be finite and nonnegative")
  index <- .select_ids(ids, names(weights), "weight sample")
  weights[index]
}

#' Summarize activities by explicit metadata groups
#'
#' Per-observation weighting is the default; cells are not donor-balanced.
#' With weights, only means are supported: SD and quantiles are NA and labelled
#' unavailable. Descriptive SD never estimates uncertainty of the group mean.
#' @param fit Native fit or view.
#' @param metadata Explicit sample metadata.
#' @param by One or more grouping columns.
#' @param id_col Sample ID column.
#' @param unit_col Optional independent-unit ID for counts.
#' @param weights Optional nonnegative named sample weights.
#' @param na_action fail or explicit omit with recorded exclusions.
#' @return Group/factor/n/n_units/mean/sd/median/q25/q75 table; no p-values.
#' @export
#' @details Origin: NEW_QOL / none.
#'   Credits metadata grouping workflow; use provenance("summarize_factor_activity")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examples
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' metadata <- data.frame(sample=rownames(X),donor=rownames(X),
#'   group=rep(c("A","B"),each=12))
#' summarize_factor_activity(view,metadata,"group",unit_col="donor")
summarize_factor_activity <- function(fit, metadata, by, id_col = "sample", unit_col = NULL, weights = NULL, na_action = "fail") {
  .validate_ids(by, "Grouping column")
  data <- .metadata_activity(fit, metadata, c(by, unit_col), id_col, na_action)
  md <- data$metadata
  w <- .weights(weights, as.character(md[[id_col]]))
  if (!nrow(md)) stop("No samples remain after metadata omissions")
  groups <- do.call(interaction, c(md[by], list(drop = TRUE, lex.order = TRUE)))
  parts <- lapply(split(seq_len(nrow(md)), groups), function(rows) {
    if (sum(w[rows]) <= 0) stop("A group has zero total weight")
    stats <- lapply(seq_len(ncol(data$A)), function(k) {
      values <- data$A[rows, k]
      quantiles <- if (is.null(weights)) stats::quantile(values, c(0.25, 0.5, 0.75), names = FALSE) else rep(NA_real_, 3)
      data.frame(factor = colnames(data$A)[k], n = length(rows),
        n_units = if (is.null(unit_col)) NA_integer_ else length(unique(md[[unit_col]][rows])),
        mean = stats::weighted.mean(values, w[rows]), sd = if (is.null(weights)) stats::sd(values) else NA_real_,
        median = quantiles[2], q25 = quantiles[1], q75 = quantiles[3])
    })
    table <- do.call(rbind, stats)
    cbind(md[rep(rows[1], nrow(table)), by, drop = FALSE], table)
  })
  out <- do.call(rbind, parts)
  out$basis_id <- rep(data$view$manifest$basis_id, nrow(out))
  rownames(out) <- NULL
  attr(out, "analysis_metadata") <- c(data$view$manifest, list(weighting = if (is.null(weights)) "per_observation" else "named_weights_means_only", weighted_quantiles = "unsupported", matching = data$matching, exclusions = data$exclusions, n_definition = "observations", unit_col = unit_col))
  .attach_provenance(out, "summarize_factor_activity", match.call(), .resolved_parameters("summarize_factor_activity", environment()))
}

#' Rank factors by descriptive in-sample metadata R squared
#'
#' Independently written stats::lm calculation, conceptually from jamboree
#' ANOVA_factors. No arbitrary top-percent filter or inferential p-values.
#' @param fit Native fit or view.
#' @param metadata Sample metadata.
#' @param variable One explanatory metadata column.
#' @param statistic Only r2 is supported.
#' @param adjust Optional adjustment columns.
#' @param id_col Sample ID column.
#' @param na_action fail or explicit omit.
#' @return Factor rank, statistic, value and sample count with basis metadata.
#' @export
#' @details Origin: ADAPTED / none.
#'   Credits ANOVA_factors / stats::lm workflow; use provenance("rank_factors_by_metadata")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examples
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' metadata <- data.frame(sample=rownames(X),donor=rownames(X),
#'   group=rep(c("A","B"),each=12))
#' rank_factors_by_metadata(view,metadata,"group")
rank_factors_by_metadata <- function(fit, metadata, variable, statistic = "r2", adjust = NULL, id_col = "sample", na_action = "fail") {
  if (statistic != "r2" || length(variable) != 1L) stop("One variable and descriptive r2 are supported")
  data <- .metadata_activity(fit, metadata, c(variable, adjust), id_col, na_action)
  .validate_ids(c(variable, adjust), "Predictor column")
  formula <- stats::reformulate(c(variable, adjust))
  X <- stats::model.matrix(formula, data$metadata)
  if (qr(X)$rank != ncol(X)) stop("Descriptive design is rank deficient")
  out <- lapply(seq_len(ncol(data$A)), function(k) {
    y <- data$A[, k]
    value <- if (sum((y - mean(y))^2) == 0) NA_real_ else {
      model <- stats::lm(y ~ X - 1)
      1 - sum(stats::residuals(model)^2) / sum((y - mean(y))^2)
    }
    data.frame(factor = colnames(data$A)[k], statistic = "r2", value = value, n_samples = length(y), basis_id = data$view$manifest$basis_id)
  })
  out <- do.call(rbind, out)
  out$rank <- rank(-out$value, ties.method = "average", na.last = "keep")
  attr(out, "analysis_metadata") <- c(data$view$manifest, list(formula = formula, matching = data$matching, exclusions = data$exclusions, inference = "descriptive_in_sample"))
  .attach_provenance(out, "rank_factors_by_metadata", match.call(), .resolved_parameters("rank_factors_by_metadata", environment()))
}

#' Compare two groups descriptively in matched factor and feature space
#'
#' Contrast c(A,B) always means A minus B. Workflow adapted independently from
#' ZemmourLib FlashierDGE. No p-values or automatic conversion to fold changes.
#' @param fit Native fit or view.
#' @param metadata Explicit sample metadata.
#' @param group Grouping column.
#' @param contrast Two distinct group labels, first minus second.
#' @param id_col Sample ID column.
#' @param unit_col Optional independent-unit column for counts.
#' @param weights Optional named observation weights.
#' @param na_action fail or explicit omit.
#' @return factor_contrast, feature_contrast, counts and matching metadata.
#' @export
#' @details Origin: ADAPTED / none.
#'   Credits FlashierDGE group contrast workflow; use provenance("compare_factor_groups")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examples
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' metadata <- data.frame(sample=rownames(X),donor=rownames(X),
#'   group=rep(c("A","B"),each=12))
#' compare_factor_groups(view,metadata,"group",c("B","A"),unit_col="donor")$factor_contrast
compare_factor_groups <- function(fit, metadata, group, contrast, id_col = "sample", unit_col = NULL, weights = NULL, na_action = "fail") {
  .validate_ids(contrast, "Contrast group")
  if (length(contrast) != 2L) stop("contrast requires exactly two different groups: first minus second")
  data <- .metadata_activity(fit, metadata, c(group, unit_col), id_col, na_action)
  md <- data$metadata
  w <- .weights(weights, as.character(md[[id_col]]))
  means <- lapply(contrast, function(label) {
    rows <- which(as.character(md[[group]]) == label)
    if (!length(rows) || sum(w[rows]) <= 0) stop("Contrast group has no observations or zero total weight: ", label)
    colSums(data$A[rows, , drop = FALSE] * w[rows]) / sum(w[rows])
  })
  beta <- means[[1]] - means[[2]]
  label <- paste(contrast, collapse = " - ")
  factor_table <- data.frame(factor = names(beta), contrast = label, estimate = as.numeric(beta), basis_id = data$view$manifest$basis_id)
  attr(factor_table, "analysis_metadata") <- data$view$manifest
  projected <- backproject_contrast(data$view, beta, "coefficients", factor_basis = data$view$manifest, contrasts = label)
  counts <- do.call(rbind, lapply(contrast, function(label) {
    rows <- which(as.character(md[[group]]) == label)
    data.frame(group = label, n_samples = length(rows), n_units = if (is.null(unit_col)) NA_integer_ else length(unique(md[[unit_col]][rows])), total_weight = sum(w[rows]))
  }))
  .attach_provenance(list(factor_contrast = factor_table, feature_contrast = projected$effects, counts = counts,
    analysis_metadata = c(data$view$manifest, list(contrast_direction = label, weighting = if (is.null(weights)) "per_observation" else "named_weights", matching = data$matching, exclusions = data$exclusions, unit_col = unit_col, inference = "descriptive"))), "compare_factor_groups", match.call(), .resolved_parameters("compare_factor_groups", environment()))
}
