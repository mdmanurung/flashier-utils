# Independent linear algebra adapted from ZemmourLib::FlashierDGE / ImmGenT.
# No upstream source copied; no unconditional joint posterior claim.
.validate_basis <- function(source, target) {
  if (inherits(source, "factor_representation")) source <- .validate_view(source)$manifest
  if (!is.list(source) || is.null(source$fit_id) || is.null(source$basis_id)) stop("Supply factor_basis as the activity representation manifest")
  if (!identical(source$fit_id, target$fit_id)) stop("Coefficients belong to a different fit")
  if (!identical(source$basis_id, target$basis_id)) stop("Incompatible coefficient basis: scale, sign, order or feature units differ")
  if (!identical(source$factor_ids, target$factor_ids) || !identical(source$feature_scale, target$feature_scale)) stop("Coefficient basis manifest is inconsistent")
  invisible(source)
}

#' Back-project a declared factor contrast into feature units
#'
#' Computes B beta using the effects matching the activity model. Offsets
#' cancel for contrasts. Point inputs produce no p-values or feature SEs.
#' Workflow origin: ZemmourLib FlashierDGE and ImmGenT; independently written.
#' @param fit Native fit or matched view.
#' @param factor_effects Named coefficients, coefficient table or draws.
#' @param input Required coefficients or draws (never inferred from shape).
#' @param representation Matched view or settings.
#' @param factor_basis Source activity manifest required for detached coefficients.
#' @param term,contrasts,conditions Explicit table selections when applicable.
#' @param allow_partial Report a restricted program contribution.
#' @param coefficient_covariance Full named coefficient covariance, when supported.
#' @param uncertainty conditional; pipeline_bootstrap requires a separate workflow.
#' @param level Interval coverage when supported.
#' @param keep_draws Retain projected draws when supported.
#' @param feature_block_size Positive block size for feature operations.
#' @return A backprojection_result with effects and analysis_metadata.
#' @export
#' @details Origin: ADAPTED / none.
#'   Credits FlashierDGE reconstruction workflow; use provenance("backproject_contrast")
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
#' backproject_contrast(view,beta,"coefficients",factor_basis=view$manifest)$effects
backproject_contrast <- function(fit, factor_effects, input, representation = NULL,
                                factor_basis = NULL, term = NULL, contrasts = NULL,
                                conditions = NULL, allow_partial = FALSE,
                                coefficient_covariance = NULL, uncertainty = "conditional",
                                level = 0.95, keep_draws = FALSE, feature_block_size = 2000L) {
  if (missing(input)) stop("input must explicitly declare coefficients or draws")
  input <- match.arg(input, c("coefficients", "draws"))
  if (uncertainty != "conditional") stop("Only fixed-program conditional uncertainty is supported")
  if (length(feature_block_size) != 1L || !is.finite(feature_block_size) || feature_block_size < 1 || feature_block_size != as.integer(feature_block_size)) stop("feature_block_size must be a positive integer")
  view <- .resolve_view(fit, representation)
  basis <- if (is.null(factor_basis)) attr(factor_effects, "analysis_metadata") else factor_basis
  .validate_basis(basis, view$manifest)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) stop("level must be between zero and one")
  if (input == "draws") return(.attach_provenance(.backproject_draws(view, factor_effects, term, contrasts, conditions, allow_partial, level, keep_draws, feature_block_size), "backproject_contrast", match.call(), .resolved_parameters("backproject_contrast", environment())))
  coefficients <- factor_effects
  if (is.data.frame(coefficients)) {
    if (!all(c("factor", "estimate") %in% names(coefficients))) stop("Coefficient table requires factor and estimate")
    for (name in c("term", "contrast", "condition")) {
      chosen <- switch(name, term = term, contrast = contrasts, condition = conditions)
      if (name %in% names(coefficients)) {
        if (is.null(chosen) && length(unique(coefficients[[name]])) > 1L) stop("Explicit ", name, " selection is required")
        if (!is.null(chosen)) coefficients <- coefficients[coefficients[[name]] %in% chosen, , drop = FALSE]
      } else if (!is.null(chosen)) stop("Coefficient table has no ", name, " column")
    }
    if (!nrow(coefficients)) stop("No coefficients match requested selection")
    label <- if ("contrast" %in% names(coefficients)) as.character(coefficients$contrast) else if ("term" %in% names(coefficients)) as.character(coefficients$term) else rep("contrast", nrow(coefficients))
    if ("condition" %in% names(coefficients)) label <- paste(label, coefficients$condition, sep = ":")
    if (anyDuplicated(paste(coefficients$factor, label, sep = "\r"))) stop("Duplicate factor/contrast coefficients")
    C <- matrix(NA_real_, length(unique(coefficients$factor)), length(unique(label)),
                dimnames = list(unique(as.character(coefficients$factor)), unique(label)))
    C[cbind(match(coefficients$factor, rownames(C)), match(label, colnames(C)))] <- coefficients$estimate
  } else if (is.numeric(coefficients) && is.null(dim(coefficients))) {
    .validate_ids(names(coefficients), "Coefficient factor")
    C <- matrix(coefficients, ncol = 1L, dimnames = list(names(coefficients), if (is.null(contrasts)) "contrast" else contrasts))
  } else C <- coefficients
  if (!is.matrix(C) || !is.numeric(C) || any(!is.finite(C))) stop("Coefficients must be a finite, complete numeric vector or matrix")
  .validate_ids(rownames(C), "Coefficient factor")
  .validate_ids(colnames(C), "Contrast")
  ids <- view$manifest$factor_ids
  if (any(!rownames(C) %in% ids)) stop("Unknown coefficient factors")
  omitted <- setdiff(ids, rownames(C))
  if (length(omitted) && !allow_partial) stop("Missing factors; use allow_partial=TRUE for a restricted contribution")
  if (!nrow(C)) stop("No factor coefficients supplied")
  B <- view$effects[, match(rownames(C), ids), drop = FALSE]
  delta <- B %*% C
  table <- data.frame(feature = rep(rownames(B), ncol(C)), contrast = rep(colnames(C), each = nrow(B)), estimate = as.vector(delta))
  metadata <- c(view$manifest, list(input_type = input, factor_basis = basis, contrast_direction = "user_declared", coefficient_basis_declaration = if (is.null(factor_basis)) "attached" else "user_supplied",
                included_factors = rownames(C), omitted_factors = omitted,
                estimand = if (length(omitted)) "partial_program_contribution" else "represented_feature_contrast",
                uncertainty_type = "none_point_coefficients", workflow_origin = "ZemmourLib::FlashierDGE / ImmGenT"))
  if (!is.null(coefficient_covariance)) {
    covariances <- if (is.matrix(coefficient_covariance)) list(coefficient_covariance) else coefficient_covariance
    if (!is.list(covariances) || length(covariances) != ncol(C)) stop("Supply one full covariance per contrast; marginal SEs are insufficient")
    if (length(covariances) > 1L) {
      .validate_ids(names(covariances), "Covariance contrast")
      covariances <- covariances[.select_ids(colnames(C), names(covariances), "covariance contrast")]
    }
    errors <- lapply(covariances, function(Sigma) {
      if (!is.matrix(Sigma) || !is.numeric(Sigma) || any(!is.finite(Sigma)) || !identical(dim(Sigma), c(nrow(C), nrow(C)))) stop("Full coefficient covariance must be a finite K by K matrix")
      .validate_ids(rownames(Sigma), "Covariance factor")
      .validate_ids(colnames(Sigma), "Covariance factor")
      Sigma <- Sigma[.select_ids(rownames(C), rownames(Sigma), "covariance factor"), .select_ids(rownames(C), colnames(Sigma), "covariance factor"), drop = FALSE]
      tolerance <- 1e-10 * max(1, max(abs(Sigma)))
      if (max(abs(Sigma - t(Sigma))) > tolerance || min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values) < -tolerance) stop("Coefficient covariance must be symmetric positive semidefinite")
      sqrt(pmax(rowSums((B %*% Sigma) * B), 0))
    })
    table$std.error <- unlist(errors, use.names = FALSE)
    critical <- stats::qnorm((1 + level) / 2)
    table$lower <- table$estimate - critical * table$std.error
    table$upper <- table$estimate + critical * table$std.error
    metadata$uncertainty_type <- "normal_approximation_coefficient_covariance_fixed_programs"
    metadata$level <- level
  }
  table$basis_id <- rep(view$manifest$basis_id, nrow(table))
  attr(table, "analysis_metadata") <- metadata
  .attach_provenance(structure(list(schema_version = "1.0.0", effects = table,
                 factor_mapping = data.frame(factor = rownames(C), native_index = match(rownames(C), ids)),
                 analysis_metadata = metadata), class = "backprojection_result"), "backproject_contrast", match.call(), .resolved_parameters("backproject_contrast", environment()))
}

.backproject_draws <- function(view, draws, term, contrasts, conditions, allow_partial, level, keep_draws, block_size) {
  label <- if (is.null(contrasts)) "contrast" else contrasts
  if (length(label) != 1L) stop("Draw mode supports one selected contrast per call")
  if (is.data.frame(draws)) {
    if (!all(c("draw", "factor", "estimate") %in% names(draws))) stop("Draw table requires draw, factor and estimate")
    for (name in c("term", "contrast", "condition")) {
      selected <- switch(name, term = term, contrast = contrasts, condition = conditions)
      if (name %in% names(draws)) {
        if (is.null(selected) && length(unique(draws[[name]])) > 1L) stop("Explicit ", name, " selection is required")
        if (!is.null(selected)) draws <- draws[draws[[name]] %in% selected, , drop = FALSE]
        if (length(unique(draws[[name]])) > 1L) stop("Draw mode supports one selected coefficient set per call")
        if (name == "contrast" && nrow(draws)) label <- as.character(draws$contrast[1])
      }
    }
    if (!nrow(draws) || anyNA(draws$draw) || anyNA(draws$factor) || anyDuplicated(paste(draws$draw, draws$factor, sep = "\r"))) stop("Draw IDs/factors must form unique complete vectors")
    D <- matrix(NA_real_, length(unique(draws$draw)), length(unique(draws$factor)), dimnames = list(as.character(unique(draws$draw)), unique(as.character(draws$factor))))
    D[cbind(match(draws$draw, rownames(D)), match(draws$factor, colnames(D)))] <- draws$estimate
  } else D <- draws
  if (!is.matrix(D) || !is.numeric(D) || nrow(D) < 2L || !ncol(D) || any(!is.finite(D))) stop("Draws must be a finite M by K matrix with at least two draws")
  .validate_ids(colnames(D), "Draw factor")
  if (is.null(rownames(D))) rownames(D) <- paste0("draw", seq_len(nrow(D)))
  .validate_ids(rownames(D), "Draw")
  .validate_ids(label, "Contrast")
  selected <- .select_ids(colnames(D), view$manifest$factor_ids, "factor")
  omitted <- setdiff(view$manifest$factor_ids, colnames(D))
  if (length(omitted) && !allow_partial) stop("Missing factors; allow_partial is required for restricted draw contributions")
  B <- view$effects[, selected, drop = FALSE]
  if (keep_draws && 8 * nrow(D) * nrow(B) > 256 * 1024^2) stop("Requested retained feature draws exceed 256 MiB; use fewer draws/features or an explicit file backend")
  tables <- list()
  retained <- if (keep_draws) matrix(NA_real_, nrow(D), nrow(B), dimnames = list(rownames(D), rownames(B))) else NULL
  for (start in seq.int(1L, nrow(B), by = block_size)) {
    rows <- start:min(nrow(B), start + block_size - 1L)
    block <- D %*% t(B[rows, , drop = FALSE])
    summary <- .draw_summary(block, level)
    tables[[length(tables) + 1L]] <- cbind(data.frame(feature = rownames(B)[rows], contrast = label, estimate = summary$mean, basis_id = view$manifest$basis_id), summary)
    if (keep_draws) retained[, rows] <- block
  }
  metadata <- c(view$manifest, list(input_type = "draws", factor_basis = view$manifest,
    included_factors = colnames(D), omitted_factors = omitted, contrast_direction = "user_declared",
    estimand = if (length(omitted)) "partial_program_contribution" else "represented_feature_contrast",
    uncertainty_type = "coefficient_draws_fixed_programs", nsamp = nrow(D), draw_ids = rownames(D), level = level,
    prob_positive_definition = "P(theta>0)", quantile_type = 7L))
  out <- do.call(rbind, tables)
  attr(out, "analysis_metadata") <- metadata
  result <- list(schema_version = "1.0.0", effects = out, factor_mapping = data.frame(factor = colnames(D), native_index = selected), analysis_metadata = metadata)
  if (keep_draws) result$draws <- retained
  structure(result, class = "backprojection_result")
}
