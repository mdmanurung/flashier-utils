# Independently implemented ImmGenT-GP filtering; upstream workflow source not copied.
#' Filter cells by total gene-program membership
#'
#' Drop cells whose summed membership across nonnegative programs is too high,
#' following the ImmGenT rule, and return rescaled activities and effects for
#' the kept cells. Reproduce the ImmGenT post-fit filter on nonnegative sample activities.
#' Each pass scales each factor to a maximum of one, then removes samples
#' whose summed membership exceeds max_val. The retained activities and feature
#' effects are rescaled together, preserving their reconstructed expression.
#' @param fit Native flashier fit or a view with nonnegative activities.
#' @param sample_side `"rows"` or `"columns"`; needed only for a bare native fit.
#' @param max_val Maximum retained total membership, inclusive (default 10).
#' @param numiter Number of passes (default 12, matching the ImmGenT pipeline).
#' @return A `membership_filter_result` list with `$activity`, `$effects`,
#'   `$offset`, `$kept_samples`, `$removed_samples`, `$history` (one row per
#'   pass) and `$converged`. Matrices use sample
#'   by factor and feature by factor axes, regardless of native sample_side.
#' @details Zero-max activity columns remain zero; their effects stay unchanged.
#'   Negative activities are rejected. The input fit is never modified or refitted.
#'   No native PVE or posterior sampler is recalculated for the selected population.
#'   A warning reports a pass limit that leaves memberships above max_val.
#'   This implements a declared selection rule, not a general cell-quality test.
#' @author Mikhael Manurung. Conceptual workflow: Ziang Zhang / AgueroZZ,
#'   ImmGenT-GP analysis; per-script authorship is unresolved.
#' @references \url{https://aguerozz.github.io/immgenT-GP-analysis/Methods_FlashierFit.html#5_Cell_filtering_by_total_GP_membership}
#' @export
#' @examples
#' L <- rbind(broad = c(10, 10, 0), next_cell = c(4, 4, .2),
#'            a = c(5, 0, 0), b = c(0, 5, 0), c = c(0, 0, 1))
#' F <- diag(3)
#' dimnames(F) <- list(paste0("g", 1:3), paste0("F", 1:3))
#' colnames(L) <- colnames(F)
#' native <- flashier::flash_factors_init(
#'   flashier::flash_init(tcrossprod(L, F)), list(L, F),
#'   ebnm_fn = list(ebnm::ebnm_point_exponential, ebnm::ebnm_normal))
#' filtered <- filter_cells_by_total_membership(native, "rows", max_val = 1.5)
#' filtered$kept_samples
#' filtered$history
filter_cells_by_total_membership <- function(fit, sample_side = NULL,
                                            max_val = 10, numiter = 12L) {
  .program_number(max_val, "max_val")
  .program_number(numiter, "numiter", integer = TRUE, minimum = 1)
  if (!is.null(sample_side)) sample_side <- match.arg(sample_side, c("rows", "columns"))
  view <- .resolve_view(fit, sample_side = sample_side)
  if (!is.null(sample_side) && sample_side != view$manifest$sample_side) {
    stop("sample_side conflicts with the representation")
  }
  A <- view$activity
  .require_biological_factors(colnames(A))
  if (!nrow(A)) stop("Membership filtering requires at least one sample")
  if (any(A < 0)) stop("Membership filtering requires nonnegative activities")
  kept <- seq_len(nrow(A))
  scaled <- A
  history <- vector("list", numiter)
  for (iteration in seq_len(numiter)) {
    maxima <- apply(scaled, 2L, max)
    maxima[maxima == 0] <- 1
    scaled <- sweep(scaled, 2L, 1 / maxima, "*")
    if (any(!is.finite(scaled))) stop("Membership rescaling produced nonfinite coordinates")
    retain <- rowSums(scaled) <= max_val
    if (!any(retain)) stop("Membership filtering would remove all samples")
    history[[iteration]] <- data.frame(iteration = iteration, before = length(kept),
      removed = sum(!retain), remaining = sum(retain))
    kept <- kept[retain]
    scaled <- scaled[retain, , drop = FALSE]
  }
  activity <- A[kept, , drop = FALSE]
  maxima <- apply(activity, 2L, max)
  zero_max <- names(maxima)[maxima == 0]
  maxima[maxima == 0] <- 1
  activity <- sweep(activity, 2L, 1 / maxima, "*")
  effects <- sweep(view$effects, 2L, maxima, "*")
  if (any(!is.finite(activity)) || any(!is.finite(effects))) {
    stop("Membership rescaling produced nonfinite coordinates")
  }
  converged <- all(rowSums(activity) <= max_val)
  if (!converged) warning("Membership filter reached numiter with memberships above max_val")
  out <- structure(list(activity = activity, effects = effects, offset = view$offset,
    kept_samples = rownames(A)[kept], removed_samples = rownames(A)[-kept],
    history = do.call(rbind, history), converged = converged), class = "membership_filter_result")
  attr(out, "analysis_metadata") <- list(source_coordinates = view$manifest,
    source_basis_id = view$manifest$basis_id, sample_side = view$manifest$sample_side,
    max_val = max_val, numiter = numiter, threshold_operator = ">",
    normalization = "per_factor_maximum", factor_scales = maxima,
    zero_max_factors = zero_max, native_fit_refitted = FALSE,
    posterior_uncertainty = "not_recomputed_for_selected_population")
  .attach_provenance(out, "filter_cells_by_total_membership", match.call(),
    .resolved_parameters("filter_cells_by_total_membership", environment()))
}
