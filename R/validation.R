# Newly written diagnostics in explicitly matching feature units.
#' Compare reconstructed and directly observed feature contrasts
#'
#' Both estimates must use the same contrast, covariates, weighting and declared
#' units. Prediction R squared is distinct from squared correlation and can be
#' negative. Regression slope and top-feature overlap are descriptive.
#' @param reconstructed Backprojection result or feature table.
#' @param observed Direct feature-effect table.
#' @param id_col Feature ID column.
#' @param estimate_col Numerical effect column.
#' @param contrast Optional single contrast label.
#' @param feature_scale Explicit common non-unknown feature scale.
#' @param top_n Positive requested overlap sizes.
#' @param zero_tolerance Nonnegative sign eligibility threshold.
#' @param plot Include scatter plot.
#' @return Metrics, matched_features, residuals, plot, exclusions, analysis_metadata.
#' @export
#' @author Mikhael Manurung.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: NEW_ORCHESTRATION / none.
#'   Credits new reconstruction diagnostics; no exact upstream function claimed; use provenance("validate_backprojection")
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
#' beta <- setNames(rep(1,ncol(view$activity)),colnames(view$activity))
#' result <- backproject_contrast(view,beta,"coefficients",factor_basis=view$manifest)
#' # Illustrative same-scale feature vector; not a calibrated statistical truth.
#' observed <- result$effects[,c("feature","estimate")]
#' validate_backprojection(result,observed,feature_scale="simulated_centered_intensity",plot=FALSE)$metrics
validate_backprojection <- function(reconstructed, observed, id_col = "feature", estimate_col = "estimate", contrast = NULL, feature_scale,
                                    top_n = c(20L, 50L, 100L), zero_tolerance = 0, plot = TRUE) {
  if (missing(feature_scale) || identical(feature_scale, "unknown")) stop("Supply matching declared feature_scale")
  if (inherits(reconstructed, "backprojection_result")) reconstructed <- reconstructed$effects
  prepare <- function(table) {
    if (!is.data.frame(table) || !all(c(id_col, estimate_col) %in% names(table))) stop("Input requires feature IDs and estimates")
    scale <- attr(table, "analysis_metadata")$feature_scale
    if (!is.null(scale) && !identical(scale, feature_scale)) stop("Feature scales differ")
    if ("contrast" %in% names(table)) {
      labels <- unique(table$contrast)
      if (is.null(contrast) && length(labels) != 1L) stop("Choose one contrast explicitly")
      if (!is.null(contrast)) table <- table[table$contrast == contrast, , drop = FALSE]
    }
    .validate_ids(as.character(table[[id_col]]), "Feature")
    if (!nrow(table) || !is.numeric(table[[estimate_col]]) || any(!is.finite(table[[estimate_col]]))) stop("Feature estimates must be finite and nonempty")
    table
  }
  r <- prepare(reconstructed); o <- prepare(observed)
  if ("contrast" %in% names(r) && "contrast" %in% names(o) && !identical(unique(as.character(r$contrast)), unique(as.character(o$contrast)))) stop("Contrast directions/labels differ")
  ids <- intersect(as.character(r[[id_col]]), as.character(o[[id_col]]))
  if (!length(ids)) stop("No common feature IDs")
  matched <- data.frame(feature = ids, reconstructed = r[[estimate_col]][match(ids, r[[id_col]])], observed = o[[estimate_col]][match(ids, o[[id_col]])])
  exclusions <- data.frame(feature = c(setdiff(r[[id_col]], ids), setdiff(o[[id_col]], ids)), reason = c(rep("missing_observed", length(setdiff(r[[id_col]], ids))), rep("missing_reconstructed", length(setdiff(o[[id_col]], ids)))))
  matched$residual <- matched$observed - matched$reconstructed
  variable <- nrow(matched) > 1L && stats::sd(matched$observed) > 0 && stats::sd(matched$reconstructed) > 0
  correlation <- if (variable) stats::cor(matched$observed, matched$reconstructed) else NA_real_
  denominator <- sum((matched$observed - mean(matched$observed))^2)
  r2 <- if (denominator > 0) 1 - sum(matched$residual^2) / denominator else NA_real_
  slope <- intercept <- NA_real_
  if (nrow(matched) > 1L && stats::sd(matched$observed) > 0) {
    model <- stats::lm(reconstructed ~ observed, data = matched)
    intercept <- unname(stats::coef(model)[1]); slope <- unname(stats::coef(model)[2])
  }
  if (length(zero_tolerance) != 1L || !is.finite(zero_tolerance) || zero_tolerance < 0) stop("zero_tolerance must be nonnegative")
  eligible <- abs(matched$observed) > zero_tolerance & abs(matched$reconstructed) > zero_tolerance
  sign <- if (any(eligible)) mean(sign(matched$observed[eligible]) == sign(matched$reconstructed[eligible])) else NA_real_
  metrics <- data.frame(metric = c("pearson_correlation", "prediction_r2", "correlation_squared", "slope", "intercept", "sign_concordance"), value = c(correlation, r2, correlation^2, slope, intercept, sign), denominator = c(rep(nrow(matched), 5L), sum(eligible)))
  metrics$reason <- ifelse(is.na(metrics$value), "constant estimates or no eligible features", NA_character_)
  if (!is.numeric(top_n) || any(!is.finite(top_n)) || any(top_n < 1) || any(top_n != as.integer(top_n))) stop("top_n must contain positive integers")
  overlap <- do.call(rbind, lapply(unique(top_n), function(n) {
    size <- min(n, nrow(matched))
    a <- utils::head(matched$feature[order(-abs(matched$observed), matched$feature)], size)
    b <- utils::head(matched$feature[order(-abs(matched$reconstructed), matched$feature)], size)
    data.frame(requested_n = n, effective_n = size, overlap_n = length(intersect(a, b)), overlap_fraction = length(intersect(a, b)) / size)
  }))
  p <- if (plot) ggplot2::ggplot(matched, ggplot2::aes(x = .data$observed, y = .data$reconstructed)) + ggplot2::geom_point() + ggplot2::geom_abline(slope = 1, intercept = 0) else NULL
  .attach_provenance(list(metrics = metrics, top_n_overlap = overlap, matched_features = matched,
       residuals = matched[, c("feature", "residual")], plot = p, exclusions = exclusions,
       analysis_metadata = list(feature_scale = feature_scale, contrast = contrast, estimand_compatibility = "user_declared_same_design_units_and_weighting", tie_rule = "absolute_effect_then_feature_ID", zero_tolerance = zero_tolerance)), "validate_backprojection", match.call(), .resolved_parameters("validate_backprojection", environment()))
}
