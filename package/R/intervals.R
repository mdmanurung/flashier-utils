# Native conditional draws, transformed using one fixed view; no Gaussian substitute.
.draw_summary <- function(values, level) {
  alpha <- (1 - level) / 2
  intervals <- t(apply(values, 2L, stats::quantile, probs = c(alpha, 1 - alpha), names = FALSE))
  data.frame(mean = colMeans(values), posterior.sd = apply(values, 2L, stats::sd),
             lower = intervals[, 1L], upper = intervals[, 2L],
             prob_positive = colMeans(values > 0), prob_zero = colMeans(values == 0))
}

#' Conditional equal-tail intervals from the native sampler
#'
#' Uses fit$sampler batches with fixed representation multipliers, conditional on
#' empirical priors. Quantiles are exact sample quantiles. Draws do not include
#' rank, optimization or biological resampling uncertainty.
#' @inheritParams factor_uncertainty
#' @param level Coverage strictly between zero and one.
#' @param nsamp At least two native draws.
#' @param seed Optional scoped seed.
#' @return Selected posterior means, SDs, quantiles and sign/zero probabilities.
#' @export
#' @details Origin: WRAPPER / call_only.
#'   Credits sampler; use provenance("factor_intervals")
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
#' factor_intervals(fit,"activity",factors="F1",ids=c("s1","s2"),nsamp=20,seed=9)
factor_intervals <- function(fit, side, representation = NULL, factors = NULL, ids = NULL,
                             level = 0.95, nsamp = 2000L, seed = NULL, sample_side = NULL) {
  if (!inherits(fit, "flash") || !is.function(fit$sampler)) stop("Native sampler is unavailable for this fit/prior")
  side <- match.arg(side, c("activity", "effects"))
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) stop("level must be between zero and one")
  if (length(nsamp) != 1L || !is.finite(nsamp) || nsamp < 2 || nsamp != as.integer(nsamp)) stop("nsamp must be an integer of at least two")
  view <- .resolve_view(fit, representation, sample_side)
  out <- .extract(view, side, factors = factors, ids = ids, format = "long")
  if (!nrow(out)) stop("No factor entries selected")
  # ponytail: 256 MiB exact-quantile storage; select fewer entries or add a file backend.
  if (8 * nsamp * nrow(out) > 256 * 1024^2) stop("Selected interval draws exceed 256 MiB; reduce nsamp or select fewer factors/IDs")
  native_side <- if ((side == "activity") == (view$manifest$sample_side == "rows")) 1L else 2L
  id_col <- if (side == "activity") "sample" else "feature"
  index <- cbind(match(out[[id_col]], rownames(view[[side]])), match(out$factor, view$manifest$factor_ids))
  multiplier <- view$manifest[[paste0(side, "_multiplier")]][index[, 2L]]
  values <- .with_seed(seed, {
    values <- matrix(NA_real_, nsamp, nrow(out))
    for (start in seq.int(1L, nsamp, by = 64L)) {
      rows <- start:min(nsamp, start + 63L)
      draws <- fit$sampler(length(rows))
      if (!is.list(draws) || length(draws) != length(rows)) stop("Unprobed native sampler layout")
      for (j in seq_along(draws)) {
        draw <- draws[[j]]
        if (!is.list(draw) || length(draw) != 2L || !identical(dim(draw[[native_side]]), dim(view[[side]]))) stop("Unprobed native sampler dimensions")
        values[rows[j], ] <- draw[[native_side]][index] * multiplier
      }
    }
    values
  })
  out <- cbind(out, .draw_summary(values, level))
  out$uncertainty_type <- rep("conditional_variational_ebmf_draws", nrow(out))
  attr(out, "analysis_metadata") <- c(view$manifest, list(level = level, nsamp = nsamp, seed = seed, interval = "equal_tail", quantile_type = 7L, prob_positive_definition = "P(theta>0)"))
  .attach_provenance(out, "factor_intervals", match.call(), .resolved_parameters("factor_intervals", environment()))
}
