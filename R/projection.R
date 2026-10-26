# Native EBMF transfer using public flashier calls. Conceptual reference-program
# transfer credit: STARCAT (Dylan Kotliar, Michelle Curtis). No STARCAT code,
# query-specific scaling or usage normalization is used.
.measurement_scale <- function(scale) {
  if (length(scale) != 1 || !is.character(scale) || is.na(scale) || !nzchar(trimws(scale)) || scale == "unknown") stop("Declare non-unknown measurement units explicitly")
}
.projection_priors <- function(view, side, priors) {
  if (is.null(priors)) {
    ids <- view$manifest$prior_identifiers
    if (is.null(ids) || length(ids) != 2 || any(ids == "custom_callback_required")) stop("Supply explicit ebnm_fn prior functions; recorded native priors are unavailable")
    priors <- lapply(ids, function(id) getExportedValue("ebnm", strsplit(id, ";", fixed = TRUE)[[1]][1]))
    if (!identical(view$manifest$prior_sample_side, side)) priors <- rev(priors)
  }
  if (is.function(priors)) priors <- rep(list(priors), 2)
  if (!is.list(priors) || length(priors) != 2 || !all(vapply(priors, is.function, logical(1)))) stop("ebnm_fn requires two native row/column prior functions")
  priors
}
.projection_prior_ids <- function(priors) {
  exports <- getNamespaceExports("ebnm")
  vapply(priors, function(fn) {
    matches <- exports[vapply(exports, function(name) identical(fn, getExportedValue("ebnm", name)), logical(1))]
    if (length(matches)) paste(sort(matches), collapse = ";") else "custom_callback_required"
  }, character(1))
}
.projection_budget <- function(bytes, budget) {
  .program_number(budget, "dense_budget_bytes", minimum = 1)
  if (bytes > budget) stop("Required dense allocations exceed dense_budget_bytes; reduce query size or increase the declared budget")
}
.projection_query <- function(query, side, S) {
  query <- .validate_matrix(query, allow_missing = FALSE)
  if (inherits(query, "sparseMatrix") && !is.null(S)) stop("Sparse observations with explicit S are unsupported by native flashier")
  if (side == "rows") query else Matrix::t(query)
}
.projection_native <- function(code) {
  warnings <- character()
  value <- withCallingHandlers(force(code), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    # Retain and emit native warnings.
  })
  list(value = value, warnings = warnings)
}
.projection_metrics <- function(Y, A, B, offset, budget) {
  # ponytail: stream reconstruction by rows; no full dense query copy is required.
  width <- ncol(Y)
  block <- max(1L, min(nrow(Y), floor(budget / (8 * max(1, width) * 4))))
  n <- so <- sp <- soo <- spp <- sop <- sse <- sae <- 0
  for (start in seq.int(1L, nrow(Y), by = block)) {
    rows <- seq.int(start, min(nrow(Y), start + block - 1L))
    observed <- as.matrix(Y[rows, , drop = FALSE])
    predicted <- sweep(tcrossprod(A[rows, , drop = FALSE], B), 2, offset, "+")
    n <- n + length(observed)
    so <- so + sum(observed)
    sp <- sp + sum(predicted)
    soo <- soo + sum(observed^2)
    spp <- spp + sum(predicted^2)
    sop <- sop + sum(observed * predicted)
    sse <- sse + sum((observed - predicted)^2)
    sae <- sae + sum(abs(observed - predicted))
  }
  sst <- max(0, soo - so^2 / n)
  pst <- max(0, spp - sp^2 / n)
  c(rmse = sqrt(sse / n), mae = sae / n, fit_r2 = if (sst > 0) 1 - sse / sst else NA_real_, correlation = if (sst > 0 && pst > 0) (sop - so * sp / n) / sqrt(sst * pst) else NA_real_, n_scored_entries = n)
}
.projection_source_subset <- function(source, factors) {
  index <- match(factors, source$manifest$factor_ids)
  source$manifest$factor_support <- .factor_support(source$manifest)[index, , drop = FALSE]
  rownames(source$manifest$factor_support) <- NULL
  source$manifest$support_version <- "1.0.0"
  for (key in c("native_indices", "activity_multiplier", "effects_multiplier", "orientation_sign", "orientation_anchor", "activity_center", "activity_scale", "native_pve")) {
    if (!is.null(source$manifest[[key]])) source$manifest[[key]] <- source$manifest[[key]][index]
  }
  source
}
.projection_S <- function(S, query, side) {
  .check_S(S)
  if (is.null(S) || length(S) == 1) {
    return(S)
  }
  native <- if (side == "rows") query else Matrix::t(query)
  if (is.matrix(S)) {
    .validate_matrix(S, allow_missing = FALSE)
    if (any(S <= 0)) stop("S must be finite and positive")
    if (!setequal(rownames(S), rownames(native)) || !setequal(colnames(S), colnames(native))) stop("S axes differ from query")
    return(S[rownames(native), colnames(native), drop = FALSE])
  }
  if (!is.numeric(S) || any(!is.finite(S)) || any(S <= 0)) stop("S must be finite and positive")
  .validate_ids(names(S), "S axis")
  row <- setequal(names(S), rownames(native))
  col <- setequal(names(S), colnames(native))
  if (row == col) stop("S vector has missing or ambiguous query axis IDs")
  S[if (row) rownames(native) else colnames(native)]
}
.projection_output <- function(native, source, A, B, offset, query, side, scale, preprocessing, settings, name, call, budget) {
  manifest <- source$manifest
  if (identical(source$manifest$engine, "gbcd")) {
    manifest$source_engine <- "gbcd"
    manifest$source_component_mapping <- source$manifest[c("mapping_scope", "engine_component_indices", "component_labels")]
    manifest$engine <- "flashier"
    manifest$mapping_scope <- "flashier"
    manifest$engine_component_indices <- NULL
    manifest$component_labels <- NULL
    fields <- c(
      "effect_source", "effect_units", "log2_fc_source", "z_score_definition",
      "lfsr_source", "uncertainty_conditioning", "baseline_definition"
    )
    manifest$source_gbcd_statistics <- manifest[fields]
    manifest[fields] <- NULL
  }
  manifest$source_basis_id <- source$manifest$basis_id
  manifest$source_coordinates <- source$manifest[c("native_indices", "scaling", "activity_multiplier", "effects_multiplier", "native_pve", "intercept_index", "offset_model", "native_factor_ids", "factor_roles", "support_version", "factor_support")]
  intercept_index <- native$value[["intercept_index"]]
  manifest$intercept_index <- intercept_index
  object <- native$value$fit
  K <- ncol(A)
  manifest$native_factor_ids <- c(if (!is.null(intercept_index)) "Intercept", if (K) colnames(A) else character())
  manifest$factor_roles <- stats::setNames(c(if (!is.null(intercept_index)) "offset", rep("factor", K)), manifest$native_factor_ids)
  manifest$native_indices <- if (is.null(intercept_index)) seq_len(K) else setdiff(seq_len(object$n_factors), intercept_index)
  records <- manifest$factor_support
  records$native_index <- manifest$native_indices
  records$role <- rep("factor", K)
  estimated_side <- if (settings$fixed == "activity") "effects" else "activity"
  records[[if (estimated_side == "effects") "effect_support" else "activity_support"]] <- rep(unname(manifest$constraints[estimated_side]), K)
  manifest$factor_support <- records
  manifest$constraints <- .support_summary(records, manifest$constraints)
  manifest$schema_version <- manifest$representation_version <- "1.1.0"
  manifest$native_pve <- object$pve[manifest$native_indices]
  manifest$scaling <- "raw"
  manifest$ldf_type <- NULL
  manifest$d_location <- NULL
  manifest$activity_multiplier <- manifest$effects_multiplier <- rep(1, K)
  manifest$orientation <- "as_fit"
  manifest$orientation_sign <- rep(1, K)
  manifest$orientation_anchor <- rep(NA_character_, K)
  manifest$activity_center <- rep(0, K)
  manifest$activity_scale <- rep(1, K)
  manifest$offset_definition <- if (!is.null(intercept_index)) "native_feature_intercept" else if (any(offset != 0)) "external_projection_offset" else "zero"
  if (!is.null(intercept_index)) manifest$offset_model <- "feature"
  manifest$sample_side <- side
  manifest$feature_scale <- scale
  manifest$preprocessing <- preprocessing
  manifest$projection <- settings
  manifest$prior_sample_side <- side
  manifest$fit_id <- digest::digest(list(source$manifest$basis_id, A, B, offset, settings, Sys.time()), algo = "sha256")
  view <- .make_representation(A, B, manifest, offset)
  attr(object, "flashbridge_metadata") <- list(sample_side = side, factor_ids = if (K) colnames(A) else character(), fit_id = view$manifest$fit_id, feature_scale = scale, preprocessing = preprocessing, prior_support = view$manifest$constraints, prior_identifiers = view$manifest$prior_identifiers, offset = if (is.null(intercept_index)) offset else NULL, offset_model = manifest$offset_model, intercept_index = intercept_index, native_indices = manifest$native_indices, source_basis_id = source$manifest$basis_id, projection = settings, recipe_version = "1.0.0", refit_policy = "callback_required")
  metadata <- attr(object, "flashbridge_metadata")
  metadata$native_factor_ids <- manifest$native_factor_ids
  metadata$factor_roles <- manifest$factor_roles
  metadata$support_version <- manifest$support_version
  metadata$factor_support <- manifest$factor_support
  attr(object, "flashbridge_metadata") <- metadata
  native_fit <- flashier::flash_fit(object)
  iteration_fields <- native_fit[grep("iter|converg|backfit", names(native_fit))]
  iteration_fields$iteration_limit <- settings$maxiter
  iteration_fields$iteration_limit_reached <- any(grepl("Maximum number of iterations reached", native$warnings, fixed = TRUE))
  metrics <- .projection_metrics(query, A, B[colnames(query), , drop = FALSE], offset[colnames(query)], budget)
  result <- list(
    fit = object, representation = view, reconstruction = metrics, warnings = native$warnings, native_iteration_diagnostics = iteration_fields,
    analysis_metadata = utils::modifyList(view$manifest, c(settings, list(engine = "flashier", source_basis_id = source$manifest$basis_id, native_warnings = native$warnings, prior_noise_policy = if (identical(settings$prior_policy, "fixed_reference")) "activity_prior_fixed_to_reference; noise_recalibrated_on_query" else "recalibrate_on_query_observations")))
  )
  .attach_provenance(result, name, call, settings)
}

#' Estimate effects for new measurements with fixed existing activities
#'
#' Given new features measured on the same samples (for example a second
#' assay), keep the fitted sample activities fixed and estimate effects for the
#' new features. For new samples use [project_factor_activity()] instead.
#' @param fit Native reference fit or a view from [standardize_factors()].
#' @param query Named matrix of the new measurements on the same samples;
#'   `sample_side` says which side holds samples.
#' @param sample_side `"rows"` or `"columns"` of `query`. Required.
#' @param feature_scale Free-text label for what the query values measure. Required.
#' @param factors Optional character factor IDs to keep fixed, in source order. Default all.
#' @param ebnm_fn Optional native prior for the new features. The recorded
#'   activity prior is reused; the default new-feature prior is signed point-Laplace.
#' @param maxiter Native iteration limit (default 200).
#' @param S Optional native standard errors; sparse plus S is refused.
#' @param dense_budget_bytes Budget for dense allocations in bytes (default 256 MiB).
#' @return List with `$fit` (native), `$representation` (view with the new
#'   features), `$reconstruction` (in-sample rmse, mae, fit_r2, correlation; not
#'   prediction) and the source basis. Cells are matched by ID; QR initialization
#'   rejects rank deficiency.
#' @seealso [project_factor_activity()]
#' @details Feature-specific residual variance is estimated. Native noise/prior
#'   hyperparameters may recalibrate. Fixed activities retain their exact values.
#'   Sources with an intercept or external offset use an augmented ones/activity
#'   basis to estimate new signed point-normal feature offsets. Offsets are kept
#'   separate from semantic factors; dependent augmented columns are refused.
#' @author Mikhael Manurung; native engine Jason Willwerscheid, Peter Carbonetto,
#'   Wei Wang and Matthew Stephens. Conceptual transfer: STARCAT authors.
#' @export
project_factor_features <- function(fit, query, sample_side, feature_scale, factors = NULL, ebnm_fn = NULL, maxiter = 200L, S = NULL, dense_budget_bytes = 256 * 1024^2) {
  side <- match.arg(sample_side, c("rows", "columns"))
  .measurement_scale(feature_scale)
  .program_number(maxiter, "maxiter", TRUE, 1)
  source <- .resolve_view(fit, sample_side = if (inherits(fit, "flash")) .resolve_sample_side(fit, if (is.null(attr(fit, "flashbridge_metadata")$sample_side)) side else NULL) else NULL)
  Y <- .projection_query(query, side, S)
  requested <- .select_ids(factors, colnames(source$activity), "factor")
  index <- which(seq_len(ncol(source$activity)) %in% requested)
  affine <- any(source$offset != 0) || identical(source$manifest[["offset_model"]], "feature") ||
    isTRUE(source$manifest[["offset_definition"]] %in% c("native_feature_intercept", "external_projection_offset"))
  if (!length(index) && !isTRUE(affine)) stop("Select at least one factor")
  if (!all(rownames(Y) %in% rownames(source$activity))) stop("Query contains cells absent from reference activities")
  rank <- length(index) + as.integer(isTRUE(affine))
  .projection_budget(8 * (length(Y) * 4 + nrow(Y) * rank * 6 + ncol(Y) * rank * 4), dense_budget_bytes)
  A <- source$activity[rownames(Y), index, drop = FALSE]
  S <- .projection_S(S, Y, side)
  design <- if (isTRUE(affine)) cbind("(Intercept)" = rep(1, nrow(A)), A) else A
  qrA <- qr(design)
  if (qrA$rank < ncol(design)) stop("Rank-deficient activity/intercept initialization; select fewer independent factors or more measured cells")
  B0 <- t(qr.coef(qrA, as.matrix(Y)))
  colnames(B0) <- colnames(design)
  if (is.null(ebnm_fn)) {
    priors <- .projection_priors(source, side, NULL)
    priors[[if (side == "rows") 2 else 1]] <- ebnm::ebnm_point_laplace
  } else {
    priors <- .projection_priors(source, side, ebnm_fn)
  }
  native <- .projection_native({
    native_data <- if (side == "rows") Y else Matrix::t(Y)
    S_dim <- if (!is.null(S) && is.null(dim(S)) && length(S) > 1) {
      if (setequal(names(S), rownames(native_data))) 1L else 2L
    } else {
      NULL
    }
    object <- flashier::flash_init(native_data, S = S, S_dim = S_dim, var_type = if (side == "rows") 2L else 1L)
    if (isTRUE(affine)) {
      offset_priors <- priors
      offset_priors[[if (side == "rows") 2 else 1]] <- ebnm::ebnm_point_normal
      pair <- list(design[, 1, drop = FALSE], B0[, 1, drop = FALSE])
      object <- flashier::flash_factors_init(object, if (side == "rows") pair else rev(pair), ebnm_fn = offset_priors)
      if (ncol(A)) {
        pair <- list(A, B0[, -1, drop = FALSE])
        object <- flashier::flash_factors_init(object, if (side == "rows") pair else rev(pair), ebnm_fn = priors)
      }
    } else {
      pair <- if (side == "rows") list(A, B0) else list(B0, A)
      object <- flashier::flash_factors_init(object, pair, ebnm_fn = priors)
    }
    object <- flashier::flash_factors_fix(object, seq_len(ncol(design)), which_dim = if (side == "rows") "loadings" else "factors")
    object <- flashier::flash_backfit(object, maxiter = maxiter, extrapolate = FALSE, verbose = 0)
    fixed <- if (side == "rows") object$L_pm else object$F_pm
    if (!isTRUE(all.equal(unname(fixed), unname(design), tolerance = 1e-12))) stop("Native projection changed fixed activities")
    B <- if (side == "rows") object$F_pm else object$L_pm
    offset <- stats::setNames(rep(0, nrow(B)), rownames(B))
    if (isTRUE(affine)) {
      offset <- B[, 1]
      B <- B[, -1, drop = FALSE]
    }
    colnames(B) <- colnames(A)
    list(fit = object, effects = B, offset = offset, intercept_index = if (isTRUE(affine)) 1L else NULL)
  })
  B <- native$value$effects
  source <- .projection_source_subset(source, colnames(A))
  source$manifest$prior_identifiers <- .projection_prior_ids(priors)
  source$manifest$constraints["effects"] <- .prior_support(priors[[if (side == "rows") 2 else 1]])
  .projection_output(native, source, A, B, native$value$offset, Y, side, feature_scale, list(scope = "new_measurements_user_declared", reference_preprocessing = source$manifest$preprocessing), list(fixed = "activity", maxiter = maxiter, dense_budget_bytes = dense_budget_bytes, var_type = if (side == "rows") 2L else 1L), "project_factor_features", match.call(), dense_budget_bytes)
}

#' Transfer fixed reference programs to a new dataset
#'
#' Keep the reference feature effects fixed and estimate activities for new
#' samples, so programs learned on one dataset can be scored on another. The
#' query must be preprocessed the same way as the reference; you declare that.
#' @inheritParams project_factor_features
#' @param query Named matrix of new samples measured on the reference features;
#'   `sample_side` says which side holds samples. Needs at least two samples.
#' @param preprocessing Named list with character `scope` and `description` of
#'   the query preprocessing, e.g. `list(scope = "within_dataset", description =
#'   "log1p library-size normalized counts")`. Required; the caller declares
#'   compatibility with the reference.
#' @param feature_match `"exact"` (default, all reference features must be
#'   present) or `"intersection"` (use the shared features, subject to `min_coverage`).
#' @param gene_mapping Optional character vector mapping query feature names
#'   (names) to reference feature names (values), one to one.
#' @param min_coverage Minimum fraction of reference features present (default 0.5).
#' @param prior_policy "query" (default) re-estimates activity priors on the
#'   query; "fixed_reference" reuses the reference fit's estimated priors.
#' @param var_type Native query residual variance setting. NULL (default) inherits
#'   the reference fit's recorded `var_type`; views without one use 0L.
#' @return List with `$fit` (native), `$representation` (view with the query
#'   samples and the exact source effects and offsets), `$reconstruction`
#'   (in-sample fit_r2, not prediction), coverage and exclusions in the
#'   settings, and `$warnings`.
#' @seealso [factor_holdout()] for held-out evaluation, [project_factor_features()]
#' @details Uses the same native fixed-program projection as holdout evaluation.
#'   No automatic preprocessing or normalized usages. Sparse projection is retained
#'   when offsets are zero; necessary dense copies obey the declared budget.
#'
#'   Query-set dependence: with `prior_policy = "query"` activity priors are
#'   re-estimated by empirical Bayes over the query rows. `"fixed_reference"`
#'   instead fixes each activity prior to the reference fit's estimate
#'   (`L_ghat`/`F_ghat`, raw native coordinates); it needs the native reference
#'   fit and refuses views. Residual variance is still estimated on the query:
#'   `var_type = 0` (and the feature-wise setting) pools it over samples, so the
#'   same cell projected within two different query sets can still receive
#'   different activities. Only `"fixed_reference"` with sample-wise variance
#'   (`var_type = 1` for row samples, `2` for column samples) makes each sample's
#'   activity independent of the other query samples. Row-chunked projection is therefore an approximation,
#'   not an exact decomposition, and is not offered; oversize queries are refused
#'   pre-flight against `dense_budget_bytes` (see
#'   `reports/FU3B_CP6_CHUNKING_DECISION.md`).
#' @author Mikhael Manurung; native flashier authors Jason Willwerscheid,
#'   Peter Carbonetto, Wei Wang and Matthew Stephens. Conceptual STARCAT credit:
#'   Dylan Kotliar and Michelle Curtis; no source reuse.
#' @export
project_factor_activity <- function(fit, query, sample_side, feature_scale, preprocessing, factors = NULL, ebnm_fn = NULL, prior_policy = c("query", "fixed_reference"), feature_match = "exact", gene_mapping = NULL, min_coverage = 0.5, maxiter = 200L, S = NULL, var_type = NULL, dense_budget_bytes = 256 * 1024^2) {
  side <- match.arg(sample_side, c("rows", "columns"))
  .measurement_scale(feature_scale)
  prior_policy <- match.arg(prior_policy)
  if (prior_policy == "fixed_reference" && !inherits(fit, "flash")) stop("fixed_reference requires the native reference fit")
  if (is.null(var_type)) var_type <- attr(fit, "flashbridge_metadata")$configuration$var_type
  if (is.null(var_type)) var_type <- 0L
  source <- .resolve_view(fit, sample_side = if (inherits(fit, "flash")) .resolve_sample_side(fit, if (is.null(attr(fit, "flashbridge_metadata")$sample_side)) side else NULL) else NULL)
  if (prior_policy == "fixed_reference") source <- standardize_factors(fit, sample_side = source$manifest$sample_side, scaling = "raw")
  if (!identical(feature_scale, source$manifest$feature_scale)) stop("Query and reference feature units are incompatible")
  if (!is.list(preprocessing) || is.null(preprocessing$scope) || !is.character(preprocessing$scope) || length(preprocessing$scope) != 1 || !nzchar(preprocessing$scope)) stop("Supply a preprocessing manifest with an explicit scope")
  feature_match <- match.arg(feature_match, c("exact", "intersection"))
  .program_number(maxiter, "maxiter", TRUE, 1)
  .program_number(min_coverage, "min_coverage")
  if (min_coverage > 1) stop("min_coverage must be <= 1")
  Y <- .projection_query(query, side, S)
  if (nrow(Y) < 2L) stop("Projection needs at least two query samples; native flashier drops a single-sample dimension")
  if (!is.null(gene_mapping)) {
    .validate_ids(names(gene_mapping), "Query mapping")
    .validate_ids(unname(gene_mapping), "Reference mapping")
    if (!all(colnames(Y) %in% names(gene_mapping))) stop("Gene mapping must cover every query gene")
    colnames(Y) <- unname(gene_mapping[colnames(Y)])
    if (!is.null(S)) stop("Explicit S with gene mapping requires pre-mapped query and S axes")
  }
  query_genes <- colnames(Y)
  genes <- intersect(rownames(source$effects), query_genes)
  coverage <- length(genes) / nrow(source$effects)
  if (feature_match == "exact" && !setequal(colnames(Y), rownames(source$effects))) stop("Exact reference-gene matching failed; explicitly choose intersection if intended")
  requested <- .select_ids(factors, colnames(source$effects), "factor")
  index <- which(seq_len(ncol(source$effects)) %in% requested)
  if (!length(index)) stop("Select at least one factor")
  .projection_budget(8 * (nrow(source$effects) * length(index) * 2 + length(genes) * length(index) * 3 + nrow(Y) * length(index) * 4 + ncol(Y) * 4), dense_budget_bytes)
  B <- source$effects[, index, drop = FALSE]
  shared <- B[genes, , drop = FALSE]
  if (coverage < min_coverage || nrow(shared) < ncol(shared) || qr(shared)$rank < ncol(shared)) stop("Insufficient or rank-deficient shared gene support; select fewer independent factors or increase gene coverage")
  if (!is.null(S) && feature_match == "intersection") stop("Explicit S with intersection requires pre-subset query and S axes and exact matching")
  .projection_budget(8 * (nrow(Y) * ncol(B) * 4 + length(shared) * 3 + 4 * ncol(Y)), dense_budget_bytes)
  Y <- Y[, genes, drop = FALSE]
  inference <- Y
  if (any(source$offset[genes] != 0)) {
    .projection_budget(8 * (length(Y) * 4 + nrow(Y) * ncol(B) * 4 + length(shared) * 3), dense_budget_bytes)
    inference <- sweep(as.matrix(Y), 2, source$offset[genes], "-")
  } else if (!inherits(Y, "sparseMatrix")) .projection_budget(8 * (length(Y) * 4 + nrow(Y) * ncol(B) * 4 + length(shared) * 3), dense_budget_bytes)
  S <- .projection_S(S, Y, side)
  priors <- .projection_priors(source, side, ebnm_fn)
  if (prior_policy == "fixed_reference" && !is.null(ebnm_fn)) {
    recorded <- tryCatch(.projection_priors(source, side, NULL), error = function(e) NULL)
    activity <- if (side == "rows") 1 else 2
    if (!is.null(recorded) && !identical(priors[[activity]], recorded[[activity]])) stop("fixed_reference activity prior must be the reference fit's ebnm family")
  }
  g <- if (prior_policy == "fixed_reference") .reference_activity_g(fit, .projection_source_subset(source, colnames(B)))
  native <- .projection_native(.native_projection(inference, list(effects = shared), side, priors, maxiter, S = S, var_type = var_type, return_fit = TRUE, g = g))
  A <- native$value$activity
  source <- .projection_source_subset(source, colnames(B))
  source$manifest$prior_identifiers <- .projection_prior_ids(priors)
  source$manifest$constraints["activity"] <- .prior_support(priors[[if (side == "rows") 1 else 2]])
  settings <- list(fixed = "effects_and_offsets", maxiter = maxiter, feature_match = feature_match, shared_genes = genes, coverage = coverage, missing_reference_genes = setdiff(rownames(B), genes), excluded_query_genes = setdiff(query_genes, genes), gene_mapping = gene_mapping, source_preprocessing = source$manifest$preprocessing, dense_budget_bytes = dense_budget_bytes, var_type = var_type, prior_policy = prior_policy)
  .projection_output(native, source, A, B, source$offset, Y, side, feature_scale, preprocessing, settings, "project_factor_activity", match.call(), dense_budget_bytes)
}
