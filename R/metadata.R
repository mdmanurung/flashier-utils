# Independent descriptive summaries; ImmGenT and jamboree group-ranking workflow credit.
.metadata_activity <- function(fit, metadata, required, id_col, na_action) {
  view <- .resolve_view(fit)
  match <- .match_ids(rownames(view$activity), metadata, id_col, extra = "drop", na_action = na_action,
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
#' Group identity uses metadata tuples, so punctuation cannot merge groups.
#' A zero-factor fit returns a typed empty summary table.
#' With an activity_threshold, raw active_n and active_proportion summaries are
#' added per factor and group using a strict greater-than comparison on the
#' resolved representation coordinates; with weights, active_weighted_proportion
#' is added as the weighted mean of the binary active indicator, while
#' active_n and active_proportion remain raw observation summaries.
#' @param fit Native fit or view.
#' @param metadata Explicit sample metadata.
#' @param by One or more grouping columns.
#' @param id_col Sample ID column.
#' @param unit_col Optional independent-unit ID for counts.
#' @param weights Optional nonnegative named sample weights.
#' @param na_action fail or explicit omit with recorded exclusions.
#' @param activity_threshold Optional single nonnegative finite threshold on the
#'   resolved representation coordinates; no display scaling is performed.
#' @param activity_direction positive, negative or absolute coordinate side
#'   compared against the threshold with a strict > operator.
#' @return Group/factor/n/n_units/mean/sd/median/q25/q75 table, plus active_n,
#'   active_proportion and (with weights) active_weighted_proportion when a
#'   threshold is supplied; no p-values.
#' @export
#' @author David Zemmour / Zemmour Lab (ImmGenT workflows).
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
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
#' summarize_factor_activity(view,metadata,"group",activity_threshold=0.5)
summarize_factor_activity <- function(fit, metadata, by, id_col = "sample", unit_col = NULL, weights = NULL, na_action = "fail", activity_threshold = NULL, activity_direction = c("positive", "negative", "absolute")) {
  .validate_ids(by, "Grouping column")
  reserved <- c("factor", "n", "n_units", "mean", "sd", "median", "q25", "q75",
    "basis_id", "active_n", "active_proportion", "active_weighted_proportion", "baseline", "centered")
  if (any(by %in% reserved)) stop("Grouping columns conflict with reserved output names; rename: ",
    paste(intersect(by, reserved), collapse = ", "))
  activity_direction <- match.arg(activity_direction)
  prevalence <- !is.null(activity_threshold)
  if (prevalence && (!is.numeric(activity_threshold) || length(activity_threshold) != 1L ||
      !is.finite(activity_threshold) || activity_threshold < 0)) {
    stop("activity_threshold must be a single nonnegative finite number")
  }
  data <- .metadata_activity(fit, metadata, c(by, unit_col), id_col, na_action)
  md <- data$metadata
  w <- .weights(weights, as.character(md[[id_col]]))
  if (!nrow(md)) stop("No samples remain after metadata omissions")
  groups <- .tuple_groups(md[by])
  parts <- lapply(split(seq_len(nrow(md)), groups), function(rows) {
    if (sum(w[rows]) <= 0) stop("A group has zero total weight")
    stats <- lapply(seq_len(ncol(data$A)), function(k) {
      values <- data$A[rows, k]
      quantiles <- if (is.null(weights)) stats::quantile(values, c(0.25, 0.5, 0.75), names = FALSE) else rep(NA_real_, 3)
      row <- data.frame(factor = colnames(data$A)[k], n = length(rows),
        n_units = if (is.null(unit_col)) NA_integer_ else length(unique(md[[unit_col]][rows])),
        mean = stats::weighted.mean(values, w[rows]), sd = if (is.null(weights)) stats::sd(values) else NA_real_,
        median = quantiles[2], q25 = quantiles[1], q75 = quantiles[3])
      if (prevalence) {
        active <- switch(activity_direction, positive = values, negative = -values,
          absolute = abs(values)) > activity_threshold
        row$active_n <- sum(active)
        row$active_proportion <- mean(active)
        if (!is.null(weights)) row$active_weighted_proportion <- sum(w[rows] * active) / sum(w[rows])
      }
      row
    })
    table <- if (length(stats)) do.call(rbind, stats) else {
      empty <- data.frame(factor = character(), n = integer(), n_units = integer(), mean = numeric(),
        sd = numeric(), median = numeric(), q25 = numeric(), q75 = numeric())
      if (prevalence) {
        empty$active_n <- integer()
        empty$active_proportion <- numeric()
        if (!is.null(weights)) empty$active_weighted_proportion <- numeric()
      }
      empty
    }
    cbind(md[rep(rows[1], nrow(table)), by, drop = FALSE], table)
  })
  out <- do.call(rbind, parts)
  out$basis_id <- rep(data$view$manifest$basis_id, nrow(out))
  rownames(out) <- NULL
  attr(out, "analysis_metadata") <- c(data$view$manifest, list(weighting = if (is.null(weights)) "per_observation" else "named_weights_means_only", weighted_quantiles = "unsupported", matching = data$matching, exclusions = data$exclusions, n_definition = "observations", unit_col = unit_col, activity_threshold = activity_threshold, activity_direction = if (prevalence) activity_direction else NULL, prevalence_comparison = if (prevalence) "strict_greater_than_resolved_coordinates" else NULL))
  .attach_provenance(out, "summarize_factor_activity", match.call(), .resolved_parameters("summarize_factor_activity", environment()))
}

#' Center grouped factor-activity means for display
#'
#' Consumes [summarize_factor_activity()] (no competing calculator) and
#' subtracts a declared per-factor baseline from each group mean:
#' `centered = mean - baseline[factor]`. The baseline is either the
#' **equal-group** mean of group means (each group counts once) or the
#' **observation-weighted** mean (group means weighted by group size `n`);
#' the choice is declared in the output metadata. Raw means, baselines,
#' centered values and group counts are all returned. Display-only
#' transformation: no inferential labeling, and centering never changes
#' ordering within a factor.
#'
#' @param fit Native fit or representation.
#' @param metadata Sample metadata.
#' @param by Grouping column(s).
#' @param baseline `equal_group` or `observation_weighted`.
#' @param id_col Sample ID column.
#' @param unit_col Optional independent-unit column (carried through).
#' @param weights Optional named observation weights (passed to
#'   [summarize_factor_activity()]; the baseline itself uses group sizes `n`).
#' @param na_action fail or explicit omit.
#' @return The [summarize_factor_activity()] table plus `baseline` and
#'   `centered` columns, with centering declared in `analysis_metadata`.
#' @author Mikhael Manurung. Group-summary workflow credits David Zemmour / Zemmour Lab.
#' @keywords internal
center_factor_activity <- function(fit, metadata, by, baseline = c("equal_group", "observation_weighted"),
                                   id_col = "sample", unit_col = NULL, weights = NULL, na_action = "fail") {
  baseline <- match.arg(baseline)
  tab <- summarize_factor_activity(fit, metadata, by, id_col = id_col, unit_col = unit_col,
                                   weights = weights, na_action = na_action)
  if (nrow(tab)) {
    spl <- split(seq_len(nrow(tab)), tab$factor)
    base <- vapply(spl, function(idx) {
      if (baseline == "equal_group") mean(tab$mean[idx]) else stats::weighted.mean(tab$mean[idx], tab$n[idx])
    }, numeric(1))
    tab$baseline <- as.numeric(base[tab$factor])
    tab$centered <- tab$mean - tab$baseline
  } else {
    tab$baseline <- numeric(0)
    tab$centered <- numeric(0)
  }
  md <- attr(tab, "analysis_metadata")
  md$centering <- list(baseline = baseline, definition = "centered = mean - baseline[factor]",
                       baseline_weighting = if (baseline == "observation_weighted") "group size n" else "none (each group once)",
                       inferential = FALSE)
  attr(tab, "analysis_metadata") <- md
  .attach_provenance(tab, "center_factor_activity", match.call(), .resolved_parameters("center_factor_activity", environment()))
}

#' Rank factors by descriptive in-sample metadata R squared
#'
#' Independently written stats::lm calculation, conceptually from jamboree
#' ANOVA_factors. No arbitrary top-percent filter or inferential p-values.
#' A zero-factor fit returns a typed empty ranking table.
#' @param fit Native fit or view.
#' @param metadata Sample metadata.
#' @param variable One explanatory metadata column.
#' @param statistic Only r2 is supported.
#' @param adjust Optional adjustment columns.
#' @param id_col Sample ID column.
#' @param na_action fail or explicit omit.
#' @return Factor rank, statistic, value and sample count with basis metadata.
#' @keywords internal
#' @author Stephens Lab; individual script authors not verified.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: ADAPTED / none.
#'   Credits ANOVA_factors / stats::lm workflow; use provenance("rank_factors_by_metadata")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examplesIf FALSE
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
  predictors <- lapply(c(variable, adjust), as.name)
  formula <- stats::as.formula(call("~", Reduce(function(x, y) call("+", x, y), predictors)))
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
  out <- if (length(out)) do.call(rbind, out) else data.frame(factor = character(),
    statistic = character(), value = numeric(), n_samples = integer(), basis_id = character())
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
#' @author David Zemmour.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
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
  .require_biological_factors(colnames(data$A))
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
