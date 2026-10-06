# Public stats correlations and audited CorShrinkData interface; no p-values.
#' Describe correlations between estimated factor activities
#'
#' Repeated cells/visits do not establish donor-level inference. Posterior score
#' uncertainty is ignored. CorShrink delegates pairwise estimation to its native
#' API and reports symmetry/PSD properties rather than claiming them universally.
#' @param fit Native fit or view.
#' @param engine Required pearson, spearman or corshrink.
#' @param factors Optional unique factor IDs.
#' @param metadata Optional aligned sample metadata retained for diagnostics.
#' @param subset Optional unique sample IDs.
#' @param use complete.obs, pairwise.complete.obs or everything for base methods.
#' @param control CorShrink native settings; bootstrap and inverse modes are unsupported.
#' @return Correlation matrix, pair counts, native result and descriptive metadata.
#' @export
#' @details Origin: WRAPPER / call_only.
#'   Credits stats::cor; CorShrink entrypoints pending audit; use provenance("factor_cor")
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
#' factor_cor(view,"pearson")$correlation
factor_cor <- function(fit, engine, factors = NULL, metadata = NULL, subset = NULL, use = "complete.obs", control = list()) {
  if (missing(engine)) stop("Choose correlation engine explicitly")
  engine <- match.arg(engine, c("pearson", "spearman", "corshrink"))
  use <- match.arg(use, c("complete.obs", "pairwise.complete.obs", "everything"))
  view <- .resolve_view(fit)
  A <- factor_activity(view, factors = factors, samples = subset)
  if (!ncol(A) || nrow(A) < 3L) stop("Correlation needs selected factors and at least three observations")
  if (!is.null(metadata)) .match_ids(rownames(view$activity), metadata)
  constant <- vapply(seq_len(ncol(A)), function(k) stats::sd(A[, k], na.rm = TRUE) == 0, logical(1))
  counts <- crossprod(!is.na(A))
  if (engine != "corshrink") {
    if (length(control)) stop("Base correlations have no additional control settings")
    correlation <- matrix(NA_real_, ncol(A), ncol(A), dimnames = list(colnames(A), colnames(A)))
    valid <- which(!constant)
    if (length(valid)) correlation[valid, valid] <- stats::cor(A[, valid, drop = FALSE], method = engine, use = use)
    native <- NULL
  } else {
    .check_engine("CorShrink", "CorShrinkData")
    if (any(constant)) stop("CorShrink constant columns are unsupported; remove them explicitly")
    if (use != "pairwise.complete.obs" && use != "complete.obs") stop("CorShrink uses complete or explicitly pairwise observations")
    if (use == "complete.obs") A <- A[stats::complete.cases(A), , drop = FALSE]
    if (!is.list(control) || (length(control) && is.null(names(control))) || any(!names(control) %in% c("cor_method", "thresh_up", "thresh_down", "tol", "dosym", "maxiter", "ash.control"))) stop("Unsupported CorShrink setting")
    if (isTRUE(control$dosym)) stop("CorShrink dosym=TRUE is unsupported: audited native branch fails on Matrix transpose")
    previous_warn <- getOption("warn")
    on.exit(options(warn = previous_warn), add = TRUE)
    native <- do.call(CorShrink::CorShrinkData, c(list(data = A, sd_boot = FALSE, type = "cor", image = "null"), control))
    if (!is.list(native) || !is.matrix(native$cor) || !identical(dim(native$cor), c(ncol(A),ncol(A)))) stop("Unprobed CorShrink result schema")
    correlation <- native$cor
    dimnames(correlation) <- list(colnames(A), colnames(A))
    counts <- crossprod(!is.na(A))
  }
  properties <- list(symmetry_error = if (all(is.finite(correlation))) max(abs(correlation-t(correlation))) else NA_real_,
    minimum_eigenvalue_of_symmetric_part = if (all(is.finite(correlation))) min(eigen((correlation+t(correlation))/2,symmetric=TRUE,only.values=TRUE)$values) else NA_real_)
  .attach_provenance(list(correlation = correlation, pair_counts = counts, constant_factors = colnames(A)[constant], native = native,
       analysis_metadata = c(view$manifest, list(engine = engine, use = use, selected_sample_ids = rownames(A),
         n_samples = nrow(A), posterior_score_uncertainty = "ignored", interpretation = "descriptive_not_donor_level_inference", matrix_properties = properties, native_warning_policy = if (engine == "corshrink") "upstream suppresses warnings internally; caller warn option restored" else "unchanged", control = control))), "factor_cor", match.call(), .resolved_parameters("factor_cor", environment()))
}
