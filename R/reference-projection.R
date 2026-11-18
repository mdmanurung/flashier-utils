# Native EBMF projection of query data onto a consensus reference (flashier-starcat M6).
# Reuses .native_projection() (R/holdout.R): public flashier calls only.
# Conceptual reference-program transfer: STARCAT (Kotliar, Curtis et al.); no source reuse.

#' Project query data onto a consensus reference
#'
#' Keep the consensus programs of a [build_factor_reference()] result fixed and
#' estimate each query sample's activity on them with native flashier
#' (`flash_factors_init()`, `flash_factors_fix()`, `flash_backfit()`). The query
#' must be in the measurement basis declared by the reference's
#' `projection_manifest`; you declare that by matching `feature_scale`.
#' @param reference A projection-ready `factor_reference`.
#' @param query Named matrix of query samples on reference features;
#'   `sample_side` says which side holds samples. Needs at least two samples.
#' @param sample_side `"rows"` or `"columns"` of `query`. Required.
#' @param feature_scale Free-text label for what the query values measure; must
#'   equal the reference `projection_manifest$feature_scale`. Required.
#' @param preprocessing Named list with a character `scope` describing the query
#'   preprocessing. Required.
#' @param programs Optional character program IDs. Default all.
#' @param ebnm_fn Activity prior: one ebnm prior function. Default: the
#'   `projection_manifest$activity_prior`. One of the two is required; its
#'   support must match the programs' activity support when recognizable.
#'   Hyperparameters are re-estimated on the query.
#' @param feature_match `"exact"` (default, every reference feature present) or
#'   `"intersection"` (shared features, subject to `min_coverage`).
#' @param gene_mapping Optional character vector mapping query feature names
#'   (names) to reference feature names (values), one to one.
#' @param min_coverage Minimum fraction of reference features present (default 0.5).
#' @param maxiter Native iteration limit (default 200).
#' @param S Optional native standard errors.
#' @param var_type Native residual variance setting (default 0).
#' @param normalize_activity `"none"` (default) or `"l1"`: also return activities
#'   divided by their row sum; allowed only for non-negative activity priors.
#'   Raw activity is always returned.
#' @param dense_budget_bytes Budget for dense allocations in bytes (default 256 MiB).
#' @return A `reference_projection` list: `$activity` (samples by programs),
#'   `$activity_normalized` and `$zero_sum_activity` (L1 only),
#'   `$representation` (a view whose effects are the fixed programs),
#'   `$reconstruction` (in-sample conditional fit diagnostics, not prediction),
#'   `$coverage`, `$missing_reference_features`, `$excluded_query_features`,
#'   `$warnings`, `$native_iteration_diagnostics`, `$reference_id` and `$fit`.
#' @details Origin: NEW_ORCHESTRATION / call_only. Uses the same fixed-program
#'   path as [project_factor_activity()] and checks that the programs did not
#'   change. Activities depend on the other query samples through the
#'   re-estimated prior and pooled noise; see [project_factor_activity()].
#' @seealso [build_factor_reference()], [project_factor_activity()]
#' @author Mikhael Manurung; native flashier authors Jason Willwerscheid, Peter
#'   Carbonetto, Wei Wang and Matthew Stephens. Conceptual STARCAT credit: Dylan
#'   Kotliar and Michelle Curtis; no source reuse.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @examples
#' set.seed(1)
#' programs <- matrix(rnorm(60), 30, 2)
#' simulate <- function(d, n = 40) {
#'   X <- tcrossprod(matrix(rnorm(2 * n), n, 2), programs) + matrix(rnorm(30 * n, sd = 0.3), n, 30)
#'   dimnames(X) <- list(paste0(d, "_s", seq_len(n)), paste0("g", 1:30))
#'   X
#' }
#' fits <- lapply(c("A", "B"), function(d) {
#'   fit_ebmf(simulate(d), "rows",
#'     max_factors = 2, seed = 1, verbose = 0, nullcheck = FALSE, backfit = TRUE,
#'     ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace),
#'     feature_scale = "simulated_intensity"
#'   )
#' })
#' names(fits) <- c("A", "B")
#' ref <- build_factor_reference(fits,
#'   comparison_manifest = list(
#'     mode = "generic", feature_scale = "simulated_intensity",
#'     preprocessing_compatibility = "user_declared"
#'   ),
#'   min_similarity = 0.8,
#'   projection_manifest = list(
#'     feature_scale = "simulated_intensity", scope = "external_reference",
#'     transformation = "identity", compatibility = "exact_declared_contract",
#'     activity_prior = "ebnm_normal"
#'   )
#' )
#' projected <- project_reference(ref, simulate("Q", 20), "rows", "simulated_intensity",
#'   preprocessing = list(scope = "identity"), maxiter = 50
#' )
#' head(projected$activity)
#' projected$reconstruction
#' @export
project_reference <- function(
  reference,
  query,
  sample_side,
  feature_scale,
  preprocessing,
  programs = NULL,
  ebnm_fn = NULL,
  feature_match = c("exact", "intersection"),
  gene_mapping = NULL,
  min_coverage = 0.5,
  maxiter = 200L,
  S = NULL,
  var_type = 0L,
  normalize_activity = c("none", "l1"),
  dense_budget_bytes = 256 * 1024^2
) {
  reference <- .validate_factor_reference(reference)
  if (!isTRUE(reference$projection_ready)) {
    stop(
      "This reference is not projection-ready: rebuild it with a projection_manifest ",
      "and programs that share one activity support"
    )
  }
  side <- match.arg(sample_side, c("rows", "columns"))
  .measurement_scale(feature_scale)
  if (!identical(feature_scale, reference$projection_manifest$feature_scale)) {
    stop("Query feature_scale differs from the reference projection_manifest$feature_scale")
  }
  if (
    !is.list(preprocessing) ||
      !is.character(preprocessing$scope) ||
      length(preprocessing$scope) != 1L ||
      is.na(preprocessing$scope) ||
      !nzchar(preprocessing$scope)
  ) {
    stop("Supply a preprocessing manifest with an explicit scope")
  }
  feature_match <- match.arg(feature_match)
  normalize_activity <- match.arg(normalize_activity)
  .program_number(maxiter, "maxiter", TRUE, 1)
  .program_number(min_coverage, "min_coverage")
  if (min_coverage > 1) {
    stop("min_coverage must be <= 1")
  }
  .check_var_type(var_type, S)
  keep <- .select_ids(programs, colnames(reference$programs), "program")
  B <- reference$programs[, seq_len(ncol(reference$programs)) %in% keep, drop = FALSE]
  if (!ncol(B)) {
    stop("Select at least one program")
  }
  support <- reference$program_constraints
  support <- support[match(colnames(B), support$program_id), , drop = FALSE]
  activity_support <- unique(support$activity_support)
  if (is.null(ebnm_fn)) {
    prior_name <- reference$projection_manifest$activity_prior
    if (is.null(prior_name)) {
      stop("Supply ebnm_fn or record projection_manifest$activity_prior; no source prior is inherited")
    }
    ebnm_fn <- getExportedValue("ebnm", prior_name)
  }
  if (!is.function(ebnm_fn)) {
    stop("ebnm_fn must be one ebnm prior function for the query activities")
  }
  prior_support <- .prior_support(ebnm_fn)
  if (prior_support != "unknown" && activity_support != "unknown" && prior_support != activity_support) {
    stop("ebnm_fn support (", prior_support, ") differs from the programs' activity support (", activity_support, ")")
  }
  if (normalize_activity == "l1" && prior_support != "nonnegative") {
    stop("normalize_activity = \"l1\" needs a recognized nonnegative activity prior")
  }
  Y <- .projection_query(query, side, S)
  if (nrow(Y) < 2L) {
    stop("Projection needs at least two query samples; native flashier drops a single-sample dimension")
  }
  if (!is.null(gene_mapping)) {
    if (!is.null(S)) {
      stop("Explicit S with gene mapping requires pre-mapped query and S axes")
    }
    .validate_ids(names(gene_mapping), "Query mapping")
    .validate_ids(unname(gene_mapping), "Reference mapping")
    if (!all(colnames(Y) %in% names(gene_mapping))) {
      stop("Gene mapping must cover every query gene")
    }
    colnames(Y) <- unname(gene_mapping[colnames(Y)])
  }
  query_features <- colnames(Y)
  features <- intersect(rownames(B), query_features)
  coverage <- length(features) / nrow(B)
  if (feature_match == "exact" && !setequal(query_features, rownames(B))) {
    stop("Exact reference-feature matching failed; explicitly choose intersection if intended")
  }
  if (!is.null(S) && feature_match == "intersection") {
    stop("Explicit S with intersection requires pre-subset query and S axes and exact matching")
  }
  .projection_budget(8 * (nrow(Y) * ncol(B) * 4 + length(features) * ncol(B) * 3 + ncol(Y) * 4), dense_budget_bytes)
  shared <- B[features, , drop = FALSE]
  if (coverage < min_coverage || nrow(shared) < ncol(shared) || qr(shared)$rank < ncol(shared)) {
    stop("Insufficient or rank-deficient reference feature support; select fewer programs or increase coverage")
  }
  Y <- Y[, features, drop = FALSE]
  S <- .projection_S(S, Y, side)
  native <- .projection_native(.native_projection(
    Y,
    list(effects = shared),
    side,
    list(ebnm_fn, ebnm_fn),
    maxiter,
    S = S,
    var_type = var_type,
    return_fit = TRUE
  ))
  A <- native$value$activity
  normalized <- NULL
  zero_sum <- NULL
  if (normalize_activity == "l1") {
    totals <- rowSums(A)
    zero_sum <- totals == 0
    normalized <- A / ifelse(zero_sum, 1, totals)
  }
  offset <- stats::setNames(rep(0, nrow(B)), rownames(B))
  constraints <- c(activity = prior_support, effects = unique(support$effect_support))
  if (length(constraints) != 2L) {
    constraints <- c(activity = prior_support, effects = "unknown")
  }
  view <- .make_representation(
    A,
    B,
    list(
      schema_version = "1.0.0",
      representation_version = "1.0.0",
      fit_id = digest::digest(list(reference$reference_id, A, features, var_type), algo = "sha256"),
      sample_side = side,
      native_indices = seq_len(ncol(B)),
      scaling = "raw",
      orientation = "as_fit",
      offset_definition = "zero",
      feature_scale = feature_scale,
      preprocessing = preprocessing,
      constraints = constraints,
      prior_identifiers = .prior_identifiers(list(ebnm_fn, ebnm_fn)),
      prior_sample_side = side,
      reference_id = reference$reference_id
    ),
    offset
  )
  object <- native$value$fit
  iteration <- list(
    iteration_limit = maxiter,
    iteration_limit_reached = any(grepl("Maximum number of iterations reached", native$warnings, fixed = TRUE))
  )
  settings <- list(
    fixed = "reference_programs",
    programs = colnames(B),
    feature_match = feature_match,
    coverage = coverage,
    min_coverage = min_coverage,
    var_type = var_type,
    maxiter = maxiter,
    normalize_activity = normalize_activity,
    prior_policy = "query_reestimated_declared_family",
    dense_budget_bytes = dense_budget_bytes
  )
  result <- structure(
    list(
      activity = A,
      activity_normalized = normalized,
      zero_sum_activity = zero_sum,
      representation = view,
      reconstruction = .projection_metrics(Y, A, shared, offset[features], dense_budget_bytes),
      coverage = coverage,
      missing_reference_features = setdiff(rownames(B), features),
      excluded_query_features = setdiff(query_features, features),
      warnings = native$warnings,
      native_iteration_diagnostics = iteration,
      reference_id = reference$reference_id,
      fit = object,
      analysis_metadata = c(
        settings,
        list(
          engine = "flashier",
          reference_id = reference$reference_id,
          reconstruction_interpretation = "in_sample_conditional_fit_not_prediction",
          feature_scale = feature_scale,
          sample_side = side
        )
      )
    ),
    class = c("reference_projection", "list")
  )
  result <- .attach_provenance(result, "project_reference", match.call(), settings)
  result$provenance <- result$analysis_metadata[c("package_version", "function", "engine", "engine_version", "date")]
  result
}
