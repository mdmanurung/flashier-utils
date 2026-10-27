# New comparison orchestration. Linear algebra: R Core Team; native EBMF:
# Willwerscheid, Carbonetto, Wang and Stephens. No optimizer source is copied.
.program_unit <- function(x) {
  scale <- max(abs(x), 0)
  if (scale == 0) return(x)
  z <- x / scale
  z / sqrt(sum(z^2))
}

.program_log_norm <- function(x) {
  scale <- max(abs(x), 0)
  if (scale == 0) return(-Inf)
  log(scale) + log(sum((x / scale)^2)) / 2
}

.program_source <- function(x) {
  if (inherits(x, "gbcd_program_result")) x <- x$representation
  if (inherits(x, c("flash", "factor_representation"))) {
    view <- .resolve_view(x)
    records <- .factor_support(view$manifest)
    out <- list(activity = view$activity, effects = view$effects, offset = view$offset,
      support = records[, c("factor", "role", "activity_support", "effect_support"), drop = FALSE],
      feature_scale = view$manifest$feature_scale, preprocessing = view$manifest$preprocessing,
      basis_id = view$manifest$basis_id, fit_id = view$manifest$fit_id,
      independent_unit = view$manifest$independent_unit, engine = view$manifest$engine,
      manifest = view$manifest)
  } else {
    if (!is.list(x) || !all(c("activity", "effects", "support", "feature_scale") %in% names(x)))
      stop("Supply a native fit, validated representation, or declared activity/effects/support/feature_scale matrices")
    A <- x$activity; B <- x$effects
    if (!is.matrix(A) || !is.numeric(A) || !is.matrix(B) || !is.numeric(B) ||
        nrow(A) < 1L || nrow(B) < 1L || ncol(A) != ncol(B) ||
        !identical(colnames(A), colnames(B)) || any(!is.finite(A)) || any(!is.finite(B)))
      stop("Declared program matrices need finite matching factor axes")
    for (ids in list(rownames(A), rownames(B), if (ncol(A)) colnames(A) else character())) .validate_ids(ids, "Program source")
    if (!identical(sort(names(x$support)), c("activity", "effects")))
      stop("support must declare activity and effects")
    expand <- function(value) {
      if (!is.character(value) || anyNA(value) || any(!value %in% c("signed", "nonnegative", "unknown")))
        stop("Invalid declared support")
      if (length(value) == 1L) return(rep(value, ncol(A)))
      if (length(value) != ncol(A) || !setequal(names(value), colnames(A)))
        stop("Per-factor support must be named by factor")
      unname(value[colnames(A)])
    }
    role <- x$roles
    if (is.null(role)) role <- stats::setNames(rep("factor", ncol(A)), colnames(A))
    if (!is.character(role) || length(role) != ncol(A) || anyNA(names(role)) || anyDuplicated(names(role)) ||
        !setequal(names(role), colnames(A)) || anyNA(role) || any(!role %in% c("factor", "initialized")))
      stop("roles must name factor or initialized for each factor")
    offset <- x$offset
    if (is.null(offset)) offset <- stats::setNames(rep(0, nrow(B)), rownames(B))
    if (!is.numeric(offset) || !identical(names(offset), rownames(B)) || any(!is.finite(offset)))
      stop("offset must be finite and named in feature order")
    records <- data.frame(factor = if (ncol(A)) colnames(A) else character(), role = unname(role[colnames(A)]),
      activity_support = expand(x$support[["activity"]]), effect_support = expand(x$support[["effects"]]))
    out <- list(activity = A, effects = B, offset = offset, support = records,
      feature_scale = x$feature_scale, preprocessing = if (is.null(x$preprocessing)) list() else x$preprocessing,
      independent_unit = x$independent_unit, engine = if (is.null(x$engine)) "declared_matrices" else x$engine)
    out$basis_id <- digest::digest(out, algo = "sha256")
    out$fit_id <- if (is.null(x$fit_id)) out$basis_id else x$fit_id
  }
  .measurement_scale(out$feature_scale)
  if (!is.list(out$preprocessing)) stop("preprocessing must be a declared list")
  if (!is.character(out$fit_id) || length(out$fit_id) != 1L || is.na(out$fit_id) || !nzchar(out$fit_id))
    stop("Invalid source fit identity")
  for (side in c("activity", "effects")) {
    support <- out$support[[if (side == "activity") "activity_support" else "effect_support"]]
    for (j in which(support == "nonnegative")) {
      if (any(out[[side]][, j] < 0))
        stop("Observed ", side, " violate declared nonnegative support")
    }
  }
  out$active <- vapply(seq_len(ncol(out$effects)), function(j)
    any(out$activity[, j] != 0) && any(out$effects[, j] != 0), logical(1))
  out
}

.program_least_squares <- function(B, b) {
  if (!ncol(B)) return(list(coefficients = numeric(), rank = 0L, residual = b))
  z <- svd(B)
  keep <- which(z$d > max(dim(B)) * max(z$d, 0) * .Machine$double.eps)
  coef <- if (length(keep)) drop(z$v[, keep, drop = FALSE] %*%
    (drop(crossprod(z$u[, keep, drop = FALSE], b)) / z$d[keep])) else rep(0, ncol(B))
  list(coefficients = coef, rank = length(keep), residual = b - drop(B %*% coef))
}

.program_bounded_projection <- function(B, b, nonnegative, tolerance = 1e-7) {
  ls <- .program_least_squares(B, b)
  if (!ncol(B) || !any(nonnegative)) return(c(ls, list(status = "optimal_linear_span", kkt = 0, objective_gap = 0)))
  free <- which(!nonnegative); bounded <- which(nonnegative)
  reduced <- cbind(b, B[, bounded, drop = FALSE])
  if (length(free)) {
    s <- svd(B[, free, drop = FALSE], nv = 0)
    keep <- which(s$d > max(dim(B[, free, drop = FALSE])) * max(s$d) * .Machine$double.eps)
    if (length(keep)) {
      reduced <- qr.qty(qr(s$u[, keep, drop = FALSE]), reduced)
      reduced <- reduced[setdiff(seq_len(nrow(B)), seq_along(keep)), , drop = FALSE]
    }
  }
  C <- reduced[, -1, drop = FALSE]; d <- reduced[, 1]
  norms <- vapply(seq_len(ncol(C)), function(j) exp(.program_log_norm(C[, j])), numeric(1))
  useful <- which(norms > 0)
  coef <- rep(0, ncol(B))
  finish <- function(beta, status, gap = NA_real_, kkt = NA_real_) {
    coef[bounded[useful]] <- beta / norms[useful]
    if (length(free)) coef[free] <- .program_least_squares(B[, free, drop = FALSE],
      b - drop(B[, bounded, drop = FALSE] %*% coef[bounded]))$coefficients
    residual <- b - drop(B %*% coef)
    if (any(!is.finite(coef)) || status %in% c("solver_failed", "solver_not_qualified")) residual[] <- NA_real_
    list(coefficients = coef, rank = ls$rank, residual = residual,
      status = status, kkt = kkt, objective_gap = gap)
  }
  if (!length(useful)) return(finish(numeric(), "optimal_free_subspace", 0, 0))
  C <- sweep(C[, useful, drop = FALSE], 2, norms[useful], "/")
  start <- pmax(.program_least_squares(C, d)$coefficients, 0)
  gradient <- function(c) 2 * drop(crossprod(C, drop(C %*% c) - d))
  result <- tryCatch(stats::optim(start,
    function(c) sum((d - drop(C %*% c))^2), gradient, method = "L-BFGS-B",
    lower = 0, control = list(maxit = 2000L, factr = 1, pgtol = 1e-11)), error = identity)
  if (inherits(result, "error") || any(!is.finite(result$par)))
    return(finish(rep(0, length(useful)), "solver_failed"))
  beta <- result$par
  candidate <- finish(beta, "candidate")
  primal <- sum(candidate$residual^2)
  g <- gradient(beta)
  error <- max(abs(beta - pmax(beta - g, 0)), 0)
  objective_tolerance <- tolerance / 10
  # A small gradient alone is unsafe on cancelling cone faces. Use a valid
  # primal/dual objective-gap bound; uncertain faces remain unassessed.
  if (is.finite(primal) && primal <= objective_tolerance)
    return(finish(beta, "certified_primal_zero_bound", primal, error))
  dual <- d - drop(C %*% beta)
  violations <- drop(crossprod(C, dual))
  if (any(violations > 0)) {
    inward <- rowSums(C)
    gain <- drop(crossprod(C, inward))
    if (all(gain > 0)) {
      adjustment <- max(pmax(violations, 0) / gain) +
        64 * .Machine$double.eps * max(1, sum(abs(dual))) / min(gain)
      dual <- dual - adjustment * inward
    }
  }
  feasible_dual <- all(drop(crossprod(C, dual)) <= 0)
  lower_bound <- 2 * sum(d * dual) - sum(dual^2)
  gap <- primal - lower_bound
  if (!is.finite(primal) || !feasible_dual || !is.finite(gap) ||
      gap < -objective_tolerance || gap > objective_tolerance)
    return(finish(beta, "solver_not_qualified", gap, error))
  finish(beta, "certified_primal_dual_gap", max(0, gap), error)
}

.program_evidence <- function(evidence, sources) {
  out <- list(stability = list(reference = NULL, target = NULL),
    prediction = NULL, state = "not_assessed")
  if (is.null(evidence)) return(out)
  if (!is.list(evidence) || (length(evidence) && (is.null(names(evidence)) || anyNA(names(evidence)) ||
      anyDuplicated(names(evidence)))) || any(!names(evidence) %in% c("reference", "target")))
    stop("evidence must name reference and/or target")
  for (side in names(evidence)) {
    entry <- evidence[[side]]
    if (!is.list(entry) || any(!names(entry) %in% c("stability", "holdout", "source_basis_id")))
      stop("Evidence entries accept stability, holdout and source_basis_id")
    if (!is.null(entry$stability)) {
      x <- entry$stability
      if (!inherits(x, "stability_result") ||
          is.null(sources[[side]]$manifest) ||
          !identical(x$analysis_metadata$reference_basis, sources[[side]]$manifest) ||
          !is.data.frame(x$summary) || !"factor" %in% names(x$summary) ||
          any(!x$summary$factor %in% colnames(sources[[side]]$effects)))
        stop("Stability evidence does not belong to this exact source basis")
      out$stability[[side]] <- list(summary = x$summary,
        config_hash = x$analysis_metadata$config_hash, unit_col = x$analysis_metadata$unit_col)
      out$state <- "bound_stability_attached"
    }
    if (!is.null(entry$holdout) && !inherits(entry$holdout, "holdout_result"))
      stop("holdout evidence must be a holdout_result")
    if (!is.null(entry$source_basis_id) && !identical(entry$source_basis_id, sources[[side]]$basis_id))
      stop("Evidence source_basis_id differs from supplied source")
  }
  a <- evidence$reference; b <- evidence$target
  if (!is.null(a$holdout) && !is.null(b$holdout)) {
    ma <- a$holdout$analysis_metadata; mb <- b$holdout$analysis_metadata
    bound <- identical(a$source_basis_id, sources$reference$basis_id) &&
      identical(b$source_basis_id, sources$target$basis_id)
    paired <- !is.null(ma$evaluation_digest) && !is.null(ma$input_digest) &&
      identical(ma$input_digest, mb$input_digest) && identical(ma$evaluation_digest, mb$evaluation_digest)
    masked <- identical(ma$evaluation_type, "masked") && identical(mb$evaluation_type, "masked")
    if (bound && !masked) stop("Predictive evidence requires masked held-out evaluation")
    if (bound && !paired) stop("Bound holdout artifacts do not prove paired input/folds/masks")
    if (bound && paired && masked) {
      out$prediction <- list(reference = a$holdout$metrics, target = b$holdout$metrics,
        input_digest = ma$input_digest, evaluation_digest = ma$evaluation_digest,
        reference_recipe = ma$config_hash, target_recipe = mb$config_hash,
        interpretation = "driver_bound_paired_evaluation_not_program_novelty")
      out$state <- "driver_bound_paired_evaluation_attached"
    } else out$state <- "unbound_holdout_not_assessed"
  }
  out
}

#' Compare saved feature programs and their joint coverage
#'
#' Ask whether the programs of one saved fit are representable by those of
#' another: pairwise loading cosines, bidirectional span and feasible (support-
#' constrained) coverage, residual gene profiles and optional bound validation
#' evidence. No fitting or biological novelty test is performed.
#' @param reference,target Native fits, validated factor representations, GBCD
#'   results, or declared lists with named activity/effects matrices, support
#'   (activity/effects), feature_scale and optional preprocessing, roles and offset.
#'   Matrix inputs are comparison sources, not invented native FLASH fits.
#' @param feature_map Optional one-to-one reference-to-target feature mapping.
#' @param min_shared_features Computational minimum shared features, default two.
#' @param activity_poles observed evaluates directions used by each source's
#'   scores; supported evaluates its declared signed/nonnegative score support.
#' @param evidence Optional reference/target lists with exact-basis stability
#'   artifacts and explicitly driver-bound paired holdout artifacts.
#' @return A prior_comparison with similarities, bidirectional coverage, residual
#'   gene profiles, exclusions, compact manifests and separate evidence states.
#' @details Effects are compared after L2 normalization on shared genes. Span
#'   coverage allows signed coefficients; feasible coverage constrains target
#'   coefficients according to declared activity support. Unknown support and
#'   failed solver certification remain unassessed. Residuals are in unit-loading
#'   coordinates, not expression variance or absence probabilities. Offsets are
#'   reported separately. Equal units do not establish equal fitting inputs.
#' @seealso [factor_replicability()] for one-to-one factor matching;
#'   [factor_holdout()] and [factor_stability()] produce the optional evidence.
#'   Tutorial: [compare priors](https://mdmanurung.github.io/flashier-utils/articles/compare-priors.html).
#' @export
#' @author Mikhael Manurung; R Core Team (linear algebra and stats::optim).
#'   Native EBMF: Jason Willwerscheid, Peter Carbonetto, Wei Wang, Matthew Stephens.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: NEW_ORCHESTRATION / none.
#'   Credits base R linear algebra and stats::optim; use provenance("compare_prior_fits")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   coverage describes loading representability, not program usage or biological replication.
#' @examples
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48), 24, 2), matrix(rnorm(32), 16, 2)) +
#'   matrix(rnorm(384, sd = 0.2), 24, 16)
#' dimnames(X) <- list(paste0("s", 1:24), paste0("g", 1:16))
#' fit <- fit_ebmf(X, "rows", max_factors = 2, seed = 4, verbose = 0,
#'   ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace),
#'   nullcheck = FALSE, backfit = TRUE, feature_scale = "simulated_centered_intensity")
#' comparison <- compare_prior_fits(fit, fit)
#' comparison$coverage[, c("source", "factor", "pole", "span_residual", "feasible_residual")]
compare_prior_fits <- function(reference, target, feature_map = NULL,
                               min_shared_features = 2L,
                               activity_poles = c("observed", "supported"), evidence = NULL) {
  activity_poles <- match.arg(activity_poles)
  .program_number(min_shared_features, "min_shared_features", TRUE, 1)
  sources <- list(reference = .program_source(reference), target = .program_source(target))
  r <- sources$reference; t <- sources$target
  if (!identical(r$feature_scale, t$feature_scale)) stop("Program feature units differ")
  if (is.null(feature_map)) feature_map <- stats::setNames(intersect(rownames(r$effects), rownames(t$effects)),
    intersect(rownames(r$effects), rownames(t$effects)))
  if (!length(feature_map)) stop("Insufficient shared features")
  .validate_ids(names(feature_map), "Reference feature mapping")
  .validate_ids(unname(feature_map), "Target feature mapping")
  if (!all(names(feature_map) %in% rownames(r$effects)) || !all(feature_map %in% rownames(t$effects)))
    stop("Feature mapping names unknown features")
  features <- sort(names(feature_map), method = "radix")
  if (length(features) < min_shared_features) stop("Insufficient shared features")
  matrices <- list(reference = r$effects[features, , drop = FALSE],
    target = t$effects[unname(feature_map[features]), , drop = FALSE])
  rownames(matrices$target) <- features
  units <- lapply(matrices, function(B) {
    result <- B
    for (j in seq_len(ncol(B))) result[, j] <- .program_unit(B[, j])
    result
  })
  raw <- crossprod(units$reference, units$target)
  raw[] <- pmax(-1, pmin(1, raw))
  if (nrow(raw)) raw[!r$active | colSums(abs(matrices$reference)) == 0, ] <- NA_real_
  if (ncol(raw)) raw[, !t$active | colSums(abs(matrices$target)) == 0] <- NA_real_
  flips <- outer(r$support$activity_support == "signed" & r$support$effect_support == "signed",
    t$support$activity_support == "signed" & t$support$effect_support == "signed", "|")
  matching <- raw; matching[flips] <- abs(matching[flips])
  pairs <- expand.grid(reference_factor = colnames(matrices$reference),
    target_factor = colnames(matrices$target), stringsAsFactors = FALSE)
  pairs$cosine <- as.vector(raw); pairs$matching_similarity <- as.vector(matching)
  pairs$reorientation_allowed <- as.vector(flips)
  shared_samples <- intersect(rownames(r$activity), rownames(t$activity))
  pairs$shared_samples <- rep(length(shared_samples), nrow(pairs))
  pairs$activity_correlation <- rep(NA_real_, nrow(pairs))
  if (length(shared_samples) > 1L) for (i in seq_len(nrow(pairs))) {
    a <- r$activity[shared_samples, pairs$reference_factor[i]]
    b <- t$activity[shared_samples, pairs$target_factor[i]]
    a <- a / max(abs(a), 1e-300); b <- b / max(abs(b), 1e-300)
    a <- a - mean(a); b <- b - mean(b)
    if (any(a != 0) && any(b != 0)) pairs$activity_correlation[i] <-
      sum(.program_unit(a) * .program_unit(b))
  }
  rows <- residuals <- exclusions <- list()
  for (side in names(sources)) {
    other <- if (side == "reference") "target" else "reference"
    src <- sources[[side]]; dst <- sources[[other]]
    keep <- dst$active & colSums(abs(matrices[[other]])) != 0
    B <- units[[other]][, keep, drop = FALSE]
    target_support <- dst$support$activity_support[keep]
    for (j in seq_len(ncol(src$activity))) {
      id <- colnames(src$activity)[j]; a <- src$activity[, j]
      active <- src$active[j] && any(matrices[[side]][, j] != 0)
      poles <- if (activity_poles == "observed") sort(unique(sign(a[a != 0]))) else
        switch(src$support$activity_support[j], signed = c(-1, 1), nonnegative = 1, unknown = numeric())
      if (!active || !length(poles)) {
        exclusions[[length(exclusions) + 1L]] <- data.frame(side = side, factor = id,
          reason = if (!active) "inactive_or_zero_shared_loading" else "unknown_source_support")
        next
      }
      scaled_a <- a / max(abs(a))
      log_magnitude <- .program_log_norm(a) + .program_log_norm(matrices[[side]][, j])
      for (pole in poles) {
        b <- pole * units[[side]][, j]
        span <- .program_least_squares(B, b)
        feasible <- if (any(target_support == "unknown")) list(residual = rep(NA_real_, length(b)),
          status = "not_assessed_unknown_target_support", kkt = NA_real_) else
          .program_bounded_projection(B, b, target_support == "nonnegative")
        sr <- sum(span$residual^2); fr <- if (all(is.finite(feasible$residual))) sum(feasible$residual^2) else NA_real_
        rows[[length(rows) + 1L]] <- data.frame(source = side, target = other, factor = id,
          role = src$support$role[j], pole = pole, span_residual = sr, feasible_residual = fr,
          target_active_factors = sum(keep), target_numerical_rank = span$rank,
          feasible_status = feasible$status, solver_kkt = feasible$kkt,
          solver_objective_gap = if (is.null(feasible$objective_gap)) NA_real_ else feasible$objective_gap,
          pole_score_fraction = sum(scaled_a[sign(a) == pole]^2) / sum(scaled_a^2),
          score_effective_cells = sum(scaled_a^2)^2 / sum(scaled_a^4),
          shared_component_log_norm = log_magnitude,
          shared_component_norm = if (log_magnitude > log(.Machine$double.xmax)) NA_real_ else exp(log_magnitude),
          evidence_state = "not_assessed", stringsAsFactors = FALSE)
        residuals[[length(residuals) + 1L]] <- data.frame(source = side, target = other,
          factor = id, pole = pole, feature = features, normalized_loading = b,
          span_residual = span$residual, feasible_residual = feasible$residual)
      }
    }
  }
  coverage <- if (length(rows)) do.call(rbind, rows) else data.frame(source = character(), target = character(),
    factor = character(), role = character(), pole = numeric(), span_residual = numeric(), feasible_residual = numeric(),
    target_active_factors = integer(), target_numerical_rank = integer(), feasible_status = character(),
    solver_kkt = numeric(), solver_objective_gap = numeric(), pole_score_fraction = numeric(), score_effective_cells = numeric(),
    shared_component_log_norm = numeric(), shared_component_norm = numeric(), evidence_state = character())
  residual <- if (length(residuals)) do.call(rbind, residuals) else data.frame(source = character(), target = character(),
    factor = character(), pole = numeric(), feature = character(), normalized_loading = numeric(), span_residual = numeric(), feasible_residual = numeric())
  excluded <- if (length(exclusions)) do.call(rbind, exclusions) else data.frame(side = character(), factor = character(), reason = character())
  attached <- .program_evidence(evidence, sources)
  metadata <- lapply(sources, function(s) list(basis_id = s$basis_id, fit_id = s$fit_id,
    feature_scale = s$feature_scale, preprocessing = s$preprocessing, support = s$support,
    n_samples = nrow(s$activity), n_features = nrow(s$effects), n_factors = ncol(s$effects),
    active_factors = sum(s$active), offset_nonzero = any(s$offset != 0)))
  feature_coverage <- data.frame(side = c("reference", "target"),
    total_features = c(nrow(r$effects), nrow(t$effects)), shared_features = length(features),
    fraction = length(features) / c(nrow(r$effects), nrow(t$effects)))
  result <- structure(list(schema_version = "1.0.0", similarities = pairs,
    similarity_matrix = raw, matching_matrix = matching, coverage = coverage,
    residual_features = residual, exclusions = excluded, feature_coverage = feature_coverage,
    feature_map = feature_map[features], evidence = attached,
    analysis_metadata = list(reference = metadata$reference, target = metadata$target,
      preprocessing_equal = length(r$preprocessing) > 0 && identical(r$preprocessing, t$preprocessing),
      controlled_input_comparison = "not_assessed", normalization = "shared_feature_L2",
      activity_poles = activity_poles, solver = "stats::optim(L-BFGS-B), free-subspace reduction, primal/dual gap<=1e-8",
      interpretation = "representability_not_program_usage_or_biological_novelty")), class = "prior_comparison")
  .attach_provenance(result, "compare_prior_fits", match.call(),
    list(min_shared_features = min_shared_features, activity_poles = activity_poles,
      evidence_state = attached$state))
}

#' @rdname compare_prior_fits
#' @export
print.prior_comparison <- function(x, ...) {
  cat("Saved program comparison:", x$analysis_metadata$reference$n_factors, "reference factors;",
    x$analysis_metadata$target$n_factors, "target factors\n")
  cat(length(x$feature_map), "shared features;", nrow(x$coverage), "assessed directions\n")
  cat("Coverage describes loading representability. Evidence:", x$evidence$state, "\n")
  invisible(x)
}

#' @rdname compare_prior_fits
#' @param x A prior_comparison result.
#' @param type `"similarity"` (default, cosine tiles with printed values),
#'   `"coverage"` (uncovered loading fraction per source factor and pole) or
#'   `"residual"` (per-feature residual for one source factor).
#' @param direction Source side for coverage/residual: `"reference"` (default) or `"target"`.
#' @param factor Source factor ID for a residual plot.
#' @param ... Passed through by the `print` and `plot` generics; unused.
#' @export
plot.prior_comparison <- function(x, type = c("similarity", "coverage", "residual"),
                                  direction = c("reference", "target"), factor = NULL, ...) {
  type <- match.arg(type); direction <- match.arg(direction)
  shared <- paste(length(x$feature_map), "shared features")
  if (type == "similarity") {
    d <- x$similarities
    d$label <- ifelse(is.na(d$cosine), "n/a", formatC(d$cosine, format = "f", digits = 2))
    return(ggplot2::ggplot(d, ggplot2::aes(.data$target_factor, .data$reference_factor, fill = .data$cosine)) +
      ggplot2::geom_tile() + ggplot2::geom_text(ggplot2::aes(label = .data$label)) +
      ggplot2::scale_fill_gradient2(low = "#D55E00", mid = "white", high = "#0072B2", limits = c(-1, 1), na.value = "grey85") +
      ggplot2::labs(x = "Target factor", y = "Reference factor", fill = "Loading cosine", subtitle = shared))
  }
  if (type == "coverage") {
    d <- x$coverage[x$coverage$source == direction, , drop = FALSE]
    n_unassessed <- sum(!is.finite(d$feasible_residual))
    long <- rbind(data.frame(factor = d$factor, pole = d$pole, metric = rep("span", nrow(d)), residual = d$span_residual),
      data.frame(factor = d$factor, pole = d$pole, metric = rep("feasible", nrow(d)), residual = d$feasible_residual))
    long <- long[is.finite(long$residual), , drop = FALSE]
    if (!nrow(long)) return(ggplot2::ggplot(data.frame(x = 0, y = 0, label = "No assessed loading directions"),
      ggplot2::aes(.data$x, .data$y, label = .data$label)) + ggplot2::geom_text() + ggplot2::theme_void())
    long$pole <- base::factor(long$pole, c(-1, 1), c("negative", "positive"))
    return(ggplot2::ggplot(long, ggplot2::aes(.data$factor, .data$residual, colour = .data$pole, shape = .data$pole)) +
      ggplot2::geom_point(size = 2.5, position = ggplot2::position_dodge(width = 0.5)) + ggplot2::facet_wrap(ggplot2::vars(.data$metric)) +
      .discrete_scale(long$pole) + ggplot2::scale_shape_manual(values = c(negative = 25, positive = 24), drop = FALSE) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(x = paste(direction, "factor"), y = "Uncovered squared loading fraction", colour = "Pole", shape = "Pole",
        subtitle = shared, caption = paste(n_unassessed, "feasible directions unassessed")))
  }
  d <- x$residual_features[x$residual_features$source == direction, , drop = FALSE]
  if (is.null(factor) || length(factor) != 1L || !factor %in% d$factor) stop("Choose an assessed source factor")
  d <- d[d$factor == factor, , drop = FALSE]
  d$pole <- base::factor(d$pole, c(-1, 1), c("negative", "positive"))
  ggplot2::ggplot(d, ggplot2::aes(.data$feature, .data$span_residual, fill = .data$pole)) +
    ggplot2::geom_col(position = "dodge") + ggplot2::coord_flip() + .discrete_scale(d$pole, "fill") +
    ggplot2::labs(x = "Feature", y = "Span residual in unit-loading coordinates", fill = "Pole", subtitle = shared)
}
