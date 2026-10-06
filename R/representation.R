# Matched coordinates around public flashier::ldf; no native source copied.
.resolve_sample_side <- function(fit, sample_side = NULL) {
  recorded <- attr(fit, "flashbridge_metadata")$sample_side
  if (is.null(sample_side)) sample_side <- recorded
  if (is.null(sample_side)) stop("Supply sample_side explicitly for a bare native fit")
  sample_side <- match.arg(sample_side, c("rows", "columns"))
  if (!is.null(recorded) && sample_side != recorded) stop("sample_side conflicts with recorded orientation")
  sample_side
}

.raw_pair <- function(fit) {
  if (!inherits(fit, "flash")) stop("Expected a native flash object")
  native <- flashier::flash_fit(fit)
  if (fit$n_factors == 0L) {
    metadata <- attr(fit, "flashbridge_metadata")
    ids <- list(metadata$row_ids, metadata$column_ids)
    # Public flash_fit returns the retained input; inspect only dimnames for K=0.
    if (any(vapply(ids, is.null, logical(1)))) ids <- dimnames(native$Y)
    if (length(ids) != 2L) stop("Empty fit has no retained IDs; use fit_ebmf with explicit IDs")
    for (i in 1:2) .validate_ids(ids[[i]], "Empty native fit")
    return(lapply(ids, function(id) matrix(numeric(), length(id), 0L,
                                         dimnames = list(id, character()))))
  }
  pair <- list(fit$L_pm, fit$F_pm)
  if (!all(vapply(pair, is.matrix, logical(1)))) stop("Native posterior means are unavailable; a named two-dimensional fit is required")
  for (i in 1:2) .validate_ids(rownames(pair[[i]]), if (i == 1L) "Row" else "Column")
  pair
}

.make_representation <- function(A, B, manifest, offset = stats::setNames(rep(0, nrow(B)), rownames(B))) {
  if (!is.matrix(A) || !is.matrix(B) || ncol(A) != ncol(B) || !identical(colnames(A), colnames(B))) stop("Activity/effects factor axes differ")
  for (ids in list(rownames(A), rownames(B), if (ncol(A)) colnames(A) else character())) .validate_ids(ids, "Representation")
  if (any(!is.finite(A)) || any(!is.finite(B))) stop("Representation contains nonfinite values")
  if (length(offset) != nrow(B) || !identical(names(offset), rownames(B)) || any(!is.finite(offset))) stop("Offset must be finite and named by feature")
  manifest$sample_ids <- rownames(A)
  manifest$feature_ids <- rownames(B)
  manifest$factor_ids <- if (ncol(A)) colnames(A) else character()
  manifest$basis_id <- NULL
  manifest$basis_id <- digest::digest(list(manifest = manifest, activity = A, effects = B, offset = offset), algo = "sha256")
  structure(list(activity = A, effects = B, offset = offset,
                 factors = data.frame(factor = if (ncol(A)) colnames(A) else character(), native_index = manifest$native_indices,
                                      constraints = rep(paste(manifest$constraints, collapse = "/"), ncol(A)),
                                      zeroed = colSums(abs(A)) == 0 | colSums(abs(B)) == 0),
                 manifest = manifest), class = "factor_representation")
}

.validate_view <- function(view) {
  if (!inherits(view, "factor_representation")) stop("Expected a factor_representation")
  rebuilt <- .make_representation(view$activity, view$effects, view$manifest, view$offset)
  if (!identical(view$manifest$basis_id, rebuilt$manifest$basis_id)) stop("Representation basis was modified; construct a new declared basis")
  view
}

.resolve_view <- function(fit, representation = NULL, sample_side = NULL) {
  if (inherits(fit, "factor_activity_percentile")) stop("Display-only percentiles cannot substitute for inferential activity")
  if (inherits(fit, "factor_representation")) {
    if (!is.null(representation)) stop("A view already defines its representation")
    return(.validate_view(fit))
  }
  if (inherits(representation, "factor_representation")) {
    view <- .validate_view(representation)
    raw <- standardize_factors(fit, sample_side = sample_side, scaling = "raw")
    if (!identical(view$manifest$fit_id, raw$manifest$fit_id)) stop("Representation belongs to a different fit")
    return(view)
  }
  if (!is.null(representation) && !is.list(representation)) stop("representation must be a view or named settings")
  do.call(standardize_factors, c(list(fit = fit, sample_side = sample_side), representation))
}

#' Create matched activity and feature coordinates
#'
#' Uses public flashier LDF normalization, placing D on exactly one side.
#' Raw and LDF views reconstruct the same posterior mean. Anchor orientation
#' requires declared signed support on both sides. No model is refitted.
#' @param fit A native flashier fit.
#' @param sample_side Explicit rows or columns for bare native fits.
#' @param scaling Raw or ldf coordinates.
#' @param type Native LDF norm: f, o, i (2, 1, m aliases).
#' @param d_location Absorb native D into effects or activity.
#' @param orientation as_fit or anchor_feature.
#' @param prior_support Named activity/effects support: signed, nonnegative, unknown.
#' @return A serializable factor_representation with basis manifest.
#' @export
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: ADAPTED / call_only.
#'   Credits ldf; use provenance("standardize_factors")
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
#' raw <- standardize_factors(fit,scaling="raw")
#' max(abs(tcrossprod(raw$activity,raw$effects)-tcrossprod(view$activity,view$effects)))
standardize_factors <- function(fit, sample_side = NULL, scaling = "ldf", type = "f",
                                d_location = "effects", orientation = "as_fit", prior_support = NULL) {
  sample_side <- .resolve_sample_side(fit, sample_side)
  scaling <- match.arg(scaling, c("ldf", "raw"))
  d_location <- match.arg(d_location, c("effects", "activity"))
  orientation <- match.arg(orientation, c("as_fit", "anchor_feature"))
  aliases <- c("2" = "f", "1" = "o", "m" = "i")
  type <- tolower(type)
  if (type %in% names(aliases)) type <- unname(aliases[type])
  type <- match.arg(type, c("f", "o", "i"))
  raw <- .raw_pair(fit)
  order <- if (sample_side == "rows") c(1L, 2L) else c(2L, 1L)
  raw <- raw[order]
  metadata <- attr(fit, "flashbridge_metadata")
  K <- ncol(raw[[1]])
  factors <- if (is.null(metadata$factor_ids)) paste0("F", seq_len(K)) else metadata$factor_ids
  if (K == 0L) factors <- character()
  pair <- raw
  multipliers <- list(rep(1, K), rep(1, K))
  if (scaling == "ldf" && K > 0L) {
    native <- flashier::ldf(fit, type = type)
    pair <- list(native$L, native$F)[order]
    target <- if (d_location == "effects") 2L else 1L
    pair[[target]] <- sweep(pair[[target]], 2L, native$D, "*")
    # The ratio uses a nonzero mean entry to recover each fixed native multiplier.
    # For an entirely zero native side use zero, matching native zero normalization.
    multipliers <- lapply(1:2, function(side) vapply(seq_len(K), function(k) {
      nonzero <- which(raw[[side]][, k] != 0)
      if (!length(nonzero)) 0 else pair[[side]][nonzero[1], k] / raw[[side]][nonzero[1], k]
    }, numeric(1)))
  }
  for (i in 1:2) colnames(pair[[i]]) <- factors
  support <- if (is.null(prior_support)) metadata$prior_support else prior_support
  if (is.null(support)) support <- c(activity = "unknown", effects = "unknown")
  if (!identical(sort(names(support)), c("activity", "effects")) || any(!support %in% c("signed", "nonnegative", "unknown"))) stop("prior_support must declare activity and effects support")
  sign <- rep(1, K)
  anchors <- rep(NA_character_, K)
  if (orientation == "anchor_feature") {
    if (any(support != "signed")) stop("Anchor flips require signed prior support on both sides")
    for (k in seq_len(K)) {
      maximum <- max(abs(pair[[2]][, k]))
      candidates <- rownames(pair[[2]])[abs(pair[[2]][, k]) == maximum]
      anchors[k] <- sort(candidates)[1]
      value <- pair[[2]][anchors[k], k]
      if (value < 0) sign[k] <- -1
    }
    for (side in 1:2) {
      pair[[side]] <- sweep(pair[[side]], 2L, sign, "*")
      multipliers[[side]] <- multipliers[[side]] * sign
    }
  }
  fit_id <- metadata$fit_id
  if (is.null(fit_id)) fit_id <- digest::digest(list(raw, sample_side), algo = "sha256")
  runtime <- .engine_metadata()
  manifest <- c(list(schema_version = "1.0.0", representation_version = "1.0.0", fit_id = fit_id,
    sample_side = sample_side, native_indices = seq_len(K), scaling = scaling,
    ldf_type = if (scaling == "ldf") type else NULL,
    d_location = if (scaling == "ldf") d_location else NULL,
    activity_multiplier = multipliers[[1]], effects_multiplier = multipliers[[2]],
    orientation = orientation, orientation_sign = sign, orientation_anchor = anchors,
    activity_center = rep(0, K), activity_scale = rep(1, K), offset_definition = "zero",
    feature_scale = if (is.null(metadata$feature_scale)) "unknown" else metadata$feature_scale,
    preprocessing = if (is.null(metadata$preprocessing)) list() else metadata$preprocessing,
    independent_unit = metadata$independent_unit, constraints = support, native_pve = fit$pve), runtime)
  .attach_provenance(.make_representation(pair[[1]], pair[[2]], manifest), "standardize_factors", match.call(), .resolved_parameters("standardize_factors", environment()))
}
