# Call-only adapter; GBCD algorithms remain in the upstream package.
.gbcd_pin <- "670e3afb50236303ddf99db33a3394abd4a711bd"

.gbcd_preflight <- function(Y, kmax, budget) {
  positive_fraction <- sum(Y > 0) / (as.double(nrow(Y)) * ncol(Y))
  lowrank <- 2 * ncol(Y) * positive_fraction < nrow(Y)
  # Allow four covariance copies, expression workspace and split-factor workspace.
  bytes <- 8 * (2 * as.double(nrow(Y)) * ncol(Y) +
    16 * (as.double(nrow(Y)) + ncol(Y)) * kmax + if (lowrank) 0 else 4 * as.double(nrow(Y))^2)
  if (!is.finite(bytes) || bytes > budget) {
    stop("GBCD predicted workspace exceeds dense_budget_bytes (", format(bytes), " bytes)")
  }
  list(strategy = if (lowrank) "lowrank_plus_sparse" else "dense_covariance",
    positive_fraction = positive_fraction, estimated_workspace_bytes = bytes,
    dense_budget_bytes = budget, memory_estimate = "preflight estimate; not a hard peak-RSS bound")
}

.gbcd_output <- function(res, Y, side, scale, preprocessing, diagnostics, keep_native) {
  if (!is.list(res)) stop("Malformed GBCD output: expected a list")
  A <- res$L
  stats <- res$F
  if (!is.matrix(A) || !is.numeric(A) || nrow(A) != nrow(Y) || ncol(A) < 2L ||
      any(!is.finite(A)) || !is.list(stats)) stop("Malformed GBCD output: need at least two finite membership columns")
  .validate_ids(colnames(A), "GBCD component")
  if (is.null(rownames(A))) rownames(A) <- rownames(Y)
  if (!identical(rownames(A), rownames(Y))) stop("GBCD sample IDs differ from input")
  for (name in c("lfc", "z_score", "lfsr")) {
    value <- stats[[name]]
    if (!is.matrix(value) || !is.numeric(value) ||
        !identical(dim(value), c(ncol(Y), ncol(A))) ||
        !identical(colnames(value), colnames(A))) stop("Malformed GBCD signature axes: ", name)
    if (is.null(rownames(value))) rownames(value) <- colnames(Y)
    if (!identical(rownames(value), colnames(Y))) stop("GBCD feature IDs differ from input")
    stats[[name]] <- value
  }
  if (any(!is.finite(stats$lfc))) stop("GBCD point effects must be finite")
  if (any(is.finite(stats$lfsr) & (stats$lfsr < 0 | stats$lfsr > 1))) stop("GBCD lfsr must lie in [0, 1]")
  available <- is.finite(stats$z_score) & is.finite(stats$lfsr)
  stats$z_score[!is.finite(stats$z_score)] <- NA_real_
  stats$lfsr[!is.finite(stats$lfsr)] <- NA_real_
  B <- log(2) * stats$lfc
  ids <- colnames(A)
  manifest <- list(schema_version = "1.1.0", representation_version = "1.1.0",
    fit_id = digest::digest(list(A, B, scale, preprocessing), algo = "sha256"),
    engine = "gbcd", mapping_scope = "gbcd_components",
    engine_component_indices = seq_along(ids), native_indices = rep(NA_integer_, length(ids)),
    component_labels = ids, sample_side = side, scaling = "raw", orientation = "as_fit",
    constraints = c(activity = "unknown", effects = "signed"),
    feature_scale = scale, preprocessing = preprocessing, native_pve = NULL,
    support_version = "1.0.0", refit_policy = "callback_required",
    factor_support = .support_records(ids, rep(NA_integer_, length(ids)),
      c(activity = "unknown", effects = "signed")),
    effect_source = "gbcd_eb_snmf", effect_units = "natural_log",
    log2_fc_source = "gbcd_eb_snmf", z_score_definition = "ashr_posterior_mean / posterior_sd",
    lfsr_source = "conditional_linear_regression_then_ashr",
    uncertainty_conditioning = "conditional_on_estimated_memberships",
    baseline_definition = "upstream label only; prior lineage unavailable; retained as factor")
  view <- .make_representation(A, B, manifest)
  signatures <- data.frame(feature = rep(rownames(B), times = ncol(B)),
    factor = rep(ids, each = nrow(B)), effect = as.vector(B),
    log2_fc = as.vector(stats$lfc), z_score = as.vector(stats$z_score),
    lfsr = as.vector(stats$lfsr), statistics_available = as.vector(available),
    basis_id = view$manifest$basis_id)
  diagnostics$unavailable_statistics <- sum(!available)
  diagnostics$n_components <- ncol(A)
  list(representation = view, signatures = signatures, diagnostics = diagnostics,
    native = if (keep_native) res else NULL, analysis_metadata = view$manifest)
}

#' Fit generalized binary covariance gene programs
#'
#' Delegates fitting to the optional, pinned gbcd engine. No normalization runs.
#' @param X Named finite numeric matrix or double sparse Matrix.
#' @param sample_side Required sample axis: rows or columns.
#' @param kmax Integer at least two; upstream initialization rank parameter.
#'   Final rank can exceed kmax; this wrapper imposes no claimed final-rank bound.
#' @param feature_scale Non-unknown units for caller-prepared normalized,
#'   natural-log expression, including log1p when declared in preprocessing.
#' @param preprocessing List with nonempty character scope and description,
#'   and \code{log_base = "e"}.
#' @param prior Recognized nonnegative ebnm prior, or a custom function carrying
#'   \code{attr(prior, "support") = "nonnegative"} as an explicit caller declaration.
#' @param dense_budget_bytes Positive workspace preflight budget; default 2 GiB.
#'   This estimates workspaces; it does not bound actual peak memory.
#' @param keep_native Retain the upstream output for inspection.
#' @param ... Native maxiter1, maxiter2, maxiter3, control and verbose arguments.
#' @return A gbcd_program_result with representation, signatures, diagnostics,
#'   optional native output and analysis_metadata. Pass the representation to
#'   factor extraction, feature ranking, summaries and plotting functions.
#' @details Effects reconstruct the natural-log input scale. Baseline is retained
#'   as a varying factor, with zero offset. All activity supports remain unknown:
#'   upstream output labels do not establish prior lineage after sorting.
#'   Signature effect/log2_fc are EB-SNMF point estimates. z_score and lfsr come
#'   from a separate conditional linear-regression/ashr stage. Unavailable
#'   conditional statistics remain NA. Native FLASH PVE and posterior uncertainty
#'   are unavailable. Automatic refitting/projection priors cannot be inferred.
#'   Install \code{stephenslab/gbcd@670e3afb50236303ddf99db33a3394abd4a711bd}.
#'   A different or unknown installed source pin is refused before fitting.
#' @author Mikhael Manurung (adapter). GBCD: Yusha Liu, Peter Carbonetto,
#'   Jason Willwerscheid and Matthew Stephens. No engine source is copied.
#' @references \url{https://github.com/stephenslab/gbcd}
#' @export
#' @examples
#' \dontrun{
#' programs <- fit_gbcd_programs(X, "rows", kmax = 5,
#'   feature_scale = "natural_log",
#'   preprocessing = list(scope = "within_dataset", description = "log1p normalized counts", log_base = "e"))
#' factor_features(programs$representation, "GEP1")
#' programs$signatures
#' }
fit_gbcd_programs <- function(X, sample_side, kmax, feature_scale, preprocessing,
                              prior = ebnm::ebnm_generalized_binary,
                              dense_budget_bytes = 2 * 1024^3, keep_native = FALSE, ...) {
  sample_side <- match.arg(sample_side, c("rows", "columns"))
  X <- .validate_matrix(X, allow_missing = FALSE)
  .program_number(kmax, "kmax", integer = TRUE, minimum = 2)
  .program_number(dense_budget_bytes, "dense_budget_bytes")
  if (dense_budget_bytes == 0) stop("dense_budget_bytes must be positive")
  .measurement_scale(feature_scale)
  if (!is.list(preprocessing) || !all(c("scope", "description") %in% names(preprocessing)) ||
      any(!vapply(preprocessing[c("scope", "description")], function(x)
        is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x)), logical(1)))) {
    stop("preprocessing requires nonempty scope and description")
  }
  if (!identical(preprocessing$log_base, "e")) stop("preprocessing requires log_base = 'e' (natural-log expression)")
  if (!is.function(prior) || (.prior_support(prior) != "nonnegative" &&
      !identical(attr(prior, "support"), "nonnegative"))) stop("prior requires declared nonnegative support")
  if (!is.logical(keep_native) || length(keep_native) != 1L || is.na(keep_native)) stop("keep_native must be TRUE or FALSE")
  dots <- list(...)
  allowed <- c("maxiter1", "maxiter2", "maxiter3", "control", "verbose")
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) ||
      anyDuplicated(names(dots)) || any(!names(dots) %in% allowed))) stop("Unknown or protected GBCD argument")
  for (name in intersect(names(dots), c("maxiter1", "maxiter2", "maxiter3"))) {
    .program_number(dots[[name]], name, integer = TRUE, minimum = 1)
  }
  if ("verbose" %in% names(dots)) .program_number(dots$verbose, "verbose", integer = TRUE)
  if ("control" %in% names(dots)) {
    ctl <- dots$control
    if (!is.list(ctl) || (length(ctl) && (is.null(names(ctl)) || anyDuplicated(names(ctl)) ||
        any(!names(ctl) %in% c("warmstart", "extrapolate", "corr_thres"))))) stop("Invalid GBCD control")
    for (name in intersect(names(ctl), c("warmstart", "extrapolate"))) {
      if (!is.logical(ctl[[name]]) || length(ctl[[name]]) != 1L || is.na(ctl[[name]])) stop("Invalid GBCD control: ", name)
    }
    if ("corr_thres" %in% names(ctl)) {
      .program_number(ctl$corr_thres, "corr_thres")
      if (ctl$corr_thres > 1) stop("corr_thres must lie in [0, 1]")
    }
  }
  Y <- if (sample_side == "rows") X else Matrix::t(X)
  diagnostics <- .gbcd_preflight(Y, kmax, dense_budget_bytes)
  if (!requireNamespace("gbcd", quietly = TRUE)) stop("Install optional gbcd at stephenslab/gbcd@", .gbcd_pin)
  info <- utils::packageDescription("gbcd")
  if (!identical(info$RemoteSha, .gbcd_pin)) stop("GBCD requires audited source pin ", .gbcd_pin, "; installed pin is unknown or different")
  warnings <- character()
  res <- tryCatch(withCallingHandlers(do.call(gbcd::fit_gbcd,
    c(list(Y = Y, Kmax = kmax, prior = prior), dots)),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w))),
    error = function(e) stop("GBCD engine failed: ", conditionMessage(e), call. = FALSE))
  diagnostics$warnings <- warnings
  diagnostics$source_commit <- info$RemoteSha
  diagnostics$version <- info$Version
  diagnostics$settings <- utils::modifyList(list(maxiter1 = 500, maxiter2 = 200,
    maxiter3 = 500, control = list(warmstart = TRUE, extrapolate = FALSE, corr_thres = .8),
    verbose = 1), dots)
  out <- .gbcd_output(res, Y, sample_side, feature_scale, preprocessing, diagnostics, keep_native)
  packages <- c("gbcd", "ashr", "flashier", "ebnm")
  out$analysis_metadata$engine_versions <- stats::setNames(lapply(packages, function(p)
    as.character(utils::packageVersion(p))), packages)
  out$analysis_metadata$source_commit_if_known <- stats::setNames(lapply(packages, function(p) {
    sha <- utils::packageDescription(p)$RemoteSha
    if (is.null(sha)) "unknown" else sha
  }), packages)
  class(out) <- "gbcd_program_result"
  .attach_provenance(out, "fit_gbcd_programs", match.call(),
    .resolved_parameters("fit_gbcd_programs", environment()))
}
