# factor_reference class, deterministic reference_id, builder and accessors
# (flashier-starcat M4-M5). Conceptual method credit: STARCAT (Kotliar, Curtis
# et al. 2025); independent R implementation.

.reference_required_fields <- c(
  "programs",
  "program_constraints",
  "programs_table",
  "members_table",
  "feature_table",
  "clusters",
  "similarity",
  "graph",
  "excluded_members",
  "feature_sets",
  "comparison_manifest",
  "projection_manifest",
  "projection_ready",
  "source_table",
  "build_settings",
  "reference_id",
  "schema_version",
  "provenance"
)
.reference_projection_keys <- c("feature_scale", "scope", "transformation", "compatibility", "activity_prior")

# Spec section 18: the declared query measurement basis.
.reference_projection_manifest <- function(manifest) {
  if (is.null(manifest)) {
    return(NULL)
  }
  if (!is.list(manifest) || is.null(names(manifest)) || anyDuplicated(names(manifest))) {
    stop("projection_manifest must be a named list")
  }
  unknown <- setdiff(names(manifest), .reference_projection_keys)
  if (length(unknown)) {
    stop("Unknown projection_manifest fields: ", paste(unknown, collapse = ", "))
  }
  .measurement_scale(manifest$feature_scale)
  for (key in c("scope", "transformation", "compatibility")) {
    value <- manifest[[key]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
      stop("projection_manifest$", key, " must be one nonempty character string")
    }
  }
  prior <- manifest$activity_prior
  if (!is.null(prior)) {
    if (
      !is.character(prior) ||
        length(prior) != 1L ||
        !prior %in% getNamespaceExports("ebnm") ||
        !startsWith(prior, "ebnm_")
    ) {
      stop(
        "projection_manifest$activity_prior must name one exported ebnm prior function, e.g. \"ebnm_point_exponential\""
      )
    }
  }
  manifest[intersect(.reference_projection_keys, names(manifest))]
}

# Permutation-invariant digest of one source: effects with sorted rows and columns.
.reference_source_digest <- function(view) {
  B <- view$effects
  B <- B[sort(rownames(B), method = "radix"), sort(colnames(B), method = "radix"), drop = FALSE]
  support <- .factor_support(view$manifest)
  support <- support[order(support$factor, method = "radix"), c("factor", "activity_support", "effect_support")]
  rownames(support) <- NULL
  digest::digest(
    list(effects = signif(B, 12), support = support, feature_scale = view$manifest$feature_scale),
    algo = "sha256"
  )
}

# Spec section 16, adapted by D024: stable scientific content only.
.reference_digest <- function(x) {
  digest::digest(
    list(
      schema_version = x$schema_version,
      programs = signif(x$programs, 12),
      program_constraints = x$program_constraints,
      membership = x$clusters[c("cluster_id", "member_id", "dataset_id")],
      reference_features = rownames(x$programs),
      source_digests = x$source_table[c("dataset_id", "source_digest")],
      comparison_manifest = x$comparison_manifest,
      build_settings = x$build_settings,
      projection_manifest = x$projection_manifest
    ),
    algo = "sha256"
  )
}

.validate_factor_reference <- function(x) {
  if (!inherits(x, "factor_reference")) {
    stop("Expected a factor_reference from build_factor_reference()")
  }
  missing <- setdiff(.reference_required_fields, names(x))
  if (length(missing)) {
    stop("factor_reference lacks fields: ", paste(missing, collapse = ", "))
  }
  if (!identical(x$schema_version, "1.0.0")) {
    stop("Unsupported factor_reference schema_version")
  }
  if (
    !is.matrix(x$programs) ||
      !identical(as.character(colnames(x$programs)), as.character(x$program_constraints$program_id)) ||
      any(!is.finite(x$programs))
  ) {
    stop("factor_reference programs and program_constraints disagree")
  }
  if (!identical(.reference_digest(x), x$reference_id)) {
    stop("factor_reference content changed after construction; rebuild it")
  }
  x
}

#' Build a consensus program reference from several datasets
#'
#' Combine factors fitted independently in several source datasets into one
#' catalog of consensus programs. Each factor's feature effects become a
#' comparison profile. All cross-dataset pairs are compared, a nearest-neighbour
#' graph links similar factors, and a greedy, density-checked merge groups them
#' with at most one factor per dataset. A program is the mean member profile on
#' a common feature basis. Source fits are not modified.
#'
#' Experimental: arguments and outputs may change in a later release.
#'
#' Two modes exist. `"generic"` suits ordinary EBMF factors: congruence on the
#' features two datasets share, L1 profiles for non-negative effects and L2
#' profiles otherwise, sign flips only between signed/signed factors, and no
#' default similarity cut-off. `"starcat2025"` follows the STARCAT catalog
#' construction (Kotliar et al. 2025): spectra summed to 1e6 on shared
#' features and divided by each dataset's feature SD, Pearson correlation on
#' the union of the two datasets' feature sets, reciprocal top-ranked
#' neighbours above 0.5, and starCAT's cross-cluster density rule (> 0.666).
#' @param fits List of native fits from [fit_ebmf()] or views from
#'   [standardize_factors()], one per source dataset, each with a declared
#'   `feature_scale` and zero offsets.
#' @param dataset_id Character dataset IDs, one per element of `fits`. Defaults
#'   to `names(fits)`.
#' @param comparison_manifest Named list declaring the comparison contract:
#'   `mode` (`"generic"` or `"starcat2025"`), `feature_scale` (must match every
#'   source) and `preprocessing_compatibility = "user_declared"`. Optional
#'   `profile_normalization` (`"auto"`) and `feature_standardization`
#'   (`"none"` in generic mode, `"source_sd"` in starcat2025 mode).
#' @param source_feature_sets,source_feature_sd For `"starcat2025"` only: lists
#'   named by dataset ID with each dataset's discovery (overdispersed) feature
#'   IDs and its named per-feature SD in the declared linear scale.
#' @param reference_features Optional feature IDs measured in every dataset to
#'   use as the reference basis. Default: all shared features (generic) or the
#'   union of source feature sets within the shared features (starcat2025).
#' @param metric Generic mode: `"congruence"` (default), `"cosine"` (the same
#'   statistic) or `"pearson"`. starcat2025 always uses Pearson.
#' @param min_shared_features Minimum features a pair must share (default 2).
#' @param min_similarity Similarity a candidate edge must strictly exceed.
#'   Required in generic mode, where no threshold is universal; 0.5 in
#'   starcat2025 mode.
#' @param rank_max Keep only each factor's `rank_max` most similar factors from
#'   other datasets. `NULL`: no rank limit (generic) or the number of datasets
#'   (starcat2025; STARCAT used 7 for its T cell catalog).
#' @param reciprocal Keep only edges chosen from both ends. Default `FALSE`
#'   (generic); starcat2025 requires `TRUE`.
#' @param min_cluster_density Merge threshold. Generic mode: accepted edges
#'   inside the merged set divided by its pair count must be at least this
#'   (default 2/3). starcat2025: directed candidate links between the two sets
#'   divided by `2 |X| |Y|` must exceed this (default 0.666).
#' @param min_datasets Datasets a cluster needs to become a program (default 2).
#'   Singletons stay in the diagnostics.
#' @param projection_manifest Optional named list declaring the query basis:
#'   `feature_scale`, `scope`, `transformation`, `compatibility`, and an
#'   optional `activity_prior` (an ebnm function name). Without it the
#'   reference is for inspection only and [project_reference()] refuses it.
#' @param keep_source_objects `TRUE` stores the source views for debugging.
#'   Default `FALSE` keeps the object small.
#' @return A `factor_reference`: `$programs` (features by programs),
#'   `$program_constraints`, `$programs_table`, `$members_table`,
#'   `$similarity` (all pairs, with reasons), `$graph`, `$clusters`,
#'   `$excluded_members`, `$source_table`, `$build_settings`,
#'   `$projection_ready` and a content hash `$reference_id`.
#' @references Kotliar D, Curtis M, Agnew R, et al. (2025). Reproducible
#'   single-cell annotation of programs underlying T cell subsets, activation
#'   states and functions. Nature Methods 22, 1964-1980.
#'   \doi{10.1038/s41592-025-02793-1}
#' @seealso [project_reference()], [reference_programs()], [factor_replicability()]
#' @author Mikhael Manurung. Conceptual method: STARCAT, Dylan Kotliar and
#'   Michelle Curtis; no source reuse.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: NEW_ORCHESTRATION / none.
#'   Use factor_provenance("build_factor_reference") for audited sources.
#'   Cross-dataset recurrence is descriptive; it does not by itself establish a
#'   biological program.
#' @examples
#' set.seed(1)
#' programs <- matrix(rnorm(60), 30, 2)
#' fits <- lapply(1:3, function(d) {
#'   X <- tcrossprod(matrix(rnorm(80), 40, 2), programs) + matrix(rnorm(1200, sd = 0.3), 40, 30)
#'   dimnames(X) <- list(paste0("d", d, "_s", 1:40), paste0("g", 1:30))
#'   fit_ebmf(X, "rows",
#'     max_factors = 2, seed = 1, verbose = 0, nullcheck = FALSE, backfit = TRUE,
#'     ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace),
#'     feature_scale = "simulated_intensity"
#'   )
#' })
#' names(fits) <- c("A", "B", "C")
#' ref <- build_factor_reference(fits,
#'   comparison_manifest = list(
#'     mode = "generic", feature_scale = "simulated_intensity",
#'     preprocessing_compatibility = "user_declared"
#'   ),
#'   min_similarity = 0.8
#' )
#' ref
#' reference_members(ref)
#' @export
build_factor_reference <- function(
  fits,
  dataset_id = names(fits),
  comparison_manifest,
  source_feature_sets = NULL,
  source_feature_sd = NULL,
  reference_features = NULL,
  metric = NULL,
  min_shared_features = 2L,
  min_similarity = NULL,
  rank_max = NULL,
  reciprocal = NULL,
  min_cluster_density = NULL,
  min_datasets = 2L,
  projection_manifest = NULL,
  keep_source_objects = FALSE
) {
  if (missing(comparison_manifest)) {
    stop("Supply comparison_manifest; the package does not assume datasets are comparable")
  }
  manifest <- .reference_manifest(comparison_manifest)
  starcat_mode <- manifest$mode == "starcat2025"
  projection <- .reference_projection_manifest(projection_manifest)
  collected <- .reference_collect_sources(fits, dataset_id, manifest)
  n_datasets <- length(collected$views)
  starcat <- NULL
  if (starcat_mode) {
    starcat <- .reference_starcat_inputs(collected$views, source_feature_sets, source_feature_sd)
  } else if (!is.null(source_feature_sets) || !is.null(source_feature_sd)) {
    stop("source_feature_sets and source_feature_sd apply only to starcat2025 mode")
  }
  if (is.null(metric)) {
    metric <- if (starcat_mode) "pearson" else "congruence"
  }
  metric <- match.arg(metric, c("congruence", "cosine", "pearson"))
  if (starcat_mode && metric != "pearson") {
    stop("starcat2025 mode uses Pearson correlation")
  }
  .program_number(min_shared_features, "min_shared_features", TRUE, 1)
  if (is.null(min_similarity)) {
    if (!starcat_mode) {
      stop("Generic mode needs an explicit min_similarity; no factor-similarity threshold is universal")
    }
    min_similarity <- 0.5
  }
  if (
    length(min_similarity) != 1L ||
      !is.numeric(min_similarity) ||
      !is.finite(min_similarity) ||
      min_similarity < -1 ||
      min_similarity >= 1
  ) {
    stop("min_similarity must be one number in [-1, 1)")
  }
  if (is.null(rank_max) && starcat_mode) {
    rank_max <- n_datasets
  }
  if (!is.null(rank_max)) {
    .program_number(rank_max, "rank_max", TRUE, 1)
  }
  if (is.null(reciprocal)) {
    reciprocal <- starcat_mode
  }
  if (!.is_flag(reciprocal)) {
    stop("reciprocal must be TRUE or FALSE")
  }
  if (starcat_mode && !reciprocal) {
    stop("starcat2025 mode uses reciprocal edges")
  }
  if (is.null(min_cluster_density)) {
    min_cluster_density <- if (starcat_mode) 0.666 else 2 / 3
  }
  if (
    length(min_cluster_density) != 1L ||
      !is.numeric(min_cluster_density) ||
      !is.finite(min_cluster_density) ||
      min_cluster_density <= 0 ||
      min_cluster_density > 1
  ) {
    stop("min_cluster_density must be one number in (0, 1]")
  }
  .program_number(min_datasets, "min_datasets", TRUE, 1)
  if (min_datasets > n_datasets) {
    stop("min_datasets exceeds the number of source datasets")
  }
  if (!.is_flag(keep_source_objects)) {
    stop("keep_source_objects must be TRUE or FALSE")
  }
  features <- .reference_feature_universe(collected$views, reference_features, starcat)
  similarity <- .reference_similarity_table(collected, metric, min_shared_features, starcat)
  graph <- .reference_build_graph(
    similarity,
    collected$members,
    min_similarity,
    rank_max,
    reciprocal,
    best_per_dataset = starcat_mode
  )
  density_rule <- if (starcat_mode) "starcat_cross" else "merged_set"
  clusters <- .reference_cluster_graph(graph, collected$members, density_rule, min_cluster_density)
  consensus <- .reference_estimate_consensus(collected, clusters, similarity, features, starcat, min_datasets)
  table <- consensus$programs_table
  is_program <- table$status == "program"
  constraints <- data.frame(
    program_id = table$program_id[is_program],
    activity_support = table$activity_support[is_program],
    effect_support = table$effect_support[is_program]
  )
  members <- collected$members
  members$cluster_id <- clusters$membership$cluster_id[match(members$member_id, clusters$membership$member_id)]
  members$program_id <- ifelse(members$cluster_id %in% constraints$program_id, members$cluster_id, NA_character_)
  sources <- collected$sources
  sources$source_digest <- vapply(collected$views[sources$dataset_id], .reference_source_digest, character(1))
  activity <- unique(constraints$activity_support)
  ready <- !is.null(projection) && nrow(constraints) > 0L && length(activity) == 1L
  result <- list(
    programs = consensus$programs,
    program_constraints = constraints,
    programs_table = table,
    members_table = members,
    feature_table = data.frame(
      feature_id = features,
      n_source_feature_sets = if (is.null(starcat)) {
        NA_integer_
      } else {
        vapply(features, function(f) sum(vapply(starcat$sets, function(s) f %in% s, logical(1))), integer(1))
      },
      row.names = NULL
    ),
    clusters = clusters$membership,
    similarity = similarity,
    graph = c(graph, list(merge_decisions = clusters$merge_decisions)),
    excluded_members = members[!members$eligible, , drop = FALSE],
    feature_sets = list(reference_features = features, source_feature_sets = starcat$sets),
    comparison_manifest = manifest,
    projection_manifest = projection,
    projection_ready = ready,
    source_table = sources,
    build_settings = list(
      algorithm_version = .reference_algorithm_version,
      metric = metric,
      min_shared_features = as.integer(min_shared_features),
      min_similarity = min_similarity,
      similarity_comparison = "strictly_greater",
      rank_max = if (is.null(rank_max)) NULL else as.integer(rank_max),
      reciprocal = reciprocal,
      best_per_target_dataset = starcat_mode,
      density_rule = density_rule,
      min_cluster_density = min_cluster_density,
      max_members_per_dataset = 1L,
      min_datasets = as.integer(min_datasets),
      tie_rule = "similarity_then_radix_member_id",
      consensus = if (starcat_mode) "mean_source_sd_spectra_no_renormalization" else "mean_profile_renormalized",
      sign_alignment = "signed_signed_only_to_first_member_then_largest_absolute_feature_positive"
    ),
    reference_id = NA_character_,
    schema_version = "1.0.0",
    provenance = list(
      package_version = as.character(utils::packageVersion("flashier.utils")),
      conceptual_source = "STARCAT (Kotliar, Curtis et al. 2025); starCAT 5399183 audited",
      source_basis_ids = stats::setNames(sources$basis_id, sources$dataset_id)
    )
  )
  if (keep_source_objects) {
    result$source_objects <- collected$views
  }
  if (!ready && !is.null(projection)) {
    result$provenance$projection_not_ready <- if (!nrow(constraints)) {
      "no_programs"
    } else {
      "heterogeneous_activity_support"
    }
  }
  class(result) <- c("factor_reference", "list")
  result$reference_id <- .reference_digest(result)
  .attach_provenance(
    result,
    "build_factor_reference",
    match.call(),
    .resolved_parameters("build_factor_reference", environment())
  )
}

#' Access a consensus reference
#'
#' Read parts of a [build_factor_reference()] result. `reference_programs()`
#' returns the features by programs matrix, `reference_members()` the member
#' table (one row per source factor, with its cluster, program and eligibility),
#' and `reference_diagnostics()` one audit table.
#'
#' Experimental: arguments and outputs may change in a later release.
#' @param reference A `factor_reference`.
#' @param programs Optional character program IDs. Default all.
#' @param what Table to return: `"programs"` (per-cluster support and
#'   heterogeneity), `"similarity"` (all pairs with reasons), `"candidates"`
#'   (directed neighbours), `"edges"` (accepted edges), `"merges"` (merge
#'   decisions), `"excluded"` (ineligible factors) or `"sources"`.
#' @return A matrix (`reference_programs()`) or a data frame.
#' @seealso [build_factor_reference()]
#' @author Mikhael Manurung.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: NEW_QOL / none. Use factor_provenance("reference_programs").
#' @examples
#' set.seed(1)
#' programs <- matrix(rnorm(60), 30, 2)
#' fits <- lapply(1:2, function(d) {
#'   X <- tcrossprod(matrix(rnorm(80), 40, 2), programs) + matrix(rnorm(1200, sd = 0.3), 40, 30)
#'   dimnames(X) <- list(paste0("d", d, "_s", 1:40), paste0("g", 1:30))
#'   fit_ebmf(X, "rows",
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
#'   min_similarity = 0.8
#' )
#' head(reference_programs(ref))
#' reference_diagnostics(ref, "merges")
#' @export
reference_programs <- function(reference, programs = NULL) {
  reference <- .validate_factor_reference(reference)
  keep <- .select_ids(programs, colnames(reference$programs), "program")
  reference$programs[, seq_len(ncol(reference$programs)) %in% keep, drop = FALSE]
}

#' @rdname reference_programs
#' @export
reference_members <- function(reference) {
  .validate_factor_reference(reference)$members_table
}

#' @rdname reference_programs
#' @export
reference_diagnostics <- function(
  reference,
  what = c("programs", "similarity", "candidates", "edges", "merges", "excluded", "sources")
) {
  reference <- .validate_factor_reference(reference)
  what <- match.arg(what)
  switch(
    what,
    programs = reference$programs_table,
    similarity = reference$similarity,
    candidates = reference$graph$directed_candidates,
    edges = reference$graph$accepted_edges,
    merges = reference$graph$merge_decisions,
    excluded = reference$excluded_members,
    sources = reference$source_table
  )
}

#' @export
print.factor_reference <- function(x, ...) {
  cat(
    "<factor_reference> ",
    x$comparison_manifest$mode,
    " mode, ",
    nrow(x$source_table),
    " datasets, ",
    ncol(x$programs),
    " programs on ",
    nrow(x$programs),
    " features\n",
    sep = ""
  )
  clusters <- x$programs_table
  cat(
    "  clusters: ",
    nrow(clusters),
    " (",
    sum(clusters$n_members == 1L),
    " singletons); excluded factors: ",
    nrow(x$excluded_members),
    "\n",
    sep = ""
  )
  cat("  projection ready: ", x$projection_ready, "\n  reference_id: ", substr(x$reference_id, 1, 16), "\n", sep = "")
  invisible(x)
}

#' @export
summary.factor_reference <- function(object, ...) {
  object <- .validate_factor_reference(object)
  programs_table <- object$programs_table
  structure(
    list(
      mode = object$comparison_manifest$mode,
      datasets = object$source_table$dataset_id,
      n_programs = ncol(object$programs),
      n_features = nrow(object$programs),
      status = table(factor(programs_table$status)),
      programs = programs_table[programs_table$status == "program", , drop = FALSE],
      projection_ready = object$projection_ready,
      reference_id = object$reference_id
    ),
    class = "summary.factor_reference"
  )
}

#' @export
print.summary.factor_reference <- function(x, ...) {
  cat("factor_reference summary (", x$mode, ")\n", sep = "")
  cat("  datasets:", paste(x$datasets, collapse = ", "), "\n")
  cat("  programs:", x$n_programs, "on", x$n_features, "features; projection ready:", x$projection_ready, "\n")
  print(x$status)
  if (nrow(x$programs)) {
    print(x$programs[c("program_id", "n_datasets", "mean_pair_similarity", "min_member_consensus_congruence")])
  }
  invisible(x)
}

#' srr standards evidence for the consensus reference layer
#'
#' @srrstats {G1.0} build_factor_reference() cites Kotliar et al. (2025, Nature Methods 22:1964-1980, doi:10.1038/s41592-025-02793-1) in @references (R/reference-class.R); the starcat2025 mode follows starCAT build_consensus_reference.py at 5399183 (flashier-starcat/07_SOURCE_AUDIT_AND_DECISIONS.md D021-D023).
#' @srrstats {G2.0} Scalar arguments of build_factor_reference() are length-checked: min_similarity, min_cluster_density (length 1, finite, bounded), rank_max, min_shared_features and min_datasets through .program_number(), reciprocal and keep_source_objects through .is_flag() (R/reference-class.R).
#' @srrstats {G2.3a} match.arg() restricts metric in build_factor_reference(), what in reference_diagnostics(), and feature_match and normalize_activity in project_reference() (R/reference-class.R, R/reference-projection.R).
#' @srrstats {G5.4} Correctness against independent computation: all-pairs congruence equals factor_replicability() congruence to 1e-12, and starcat2025 Pearson and consensus spectra equal a by-hand implementation of the STARCAT formula (tests/testthat/test-reference.R).
#' @srrstats {G5.6} Parameter recovery: noiseless planted shared programs are recovered with congruence 1 and dataset-specific programs stay singletons; projected activities correlate above 0.99 with the planted truth (tests/testthat/test-reference.R).
#' @srrstats {UL3.3} project_reference() in R/reference-projection.R estimates activities for new samples on fixed consensus programs without refitting the sources.
#' @srrstats {UL7.0} Fail-closed tests cover a missing comparison contract, missing min_similarity, mismatched feature_scale, nonzero offsets, invalid starcat2025 inputs and non-projection-ready references (tests/testthat/test-reference.R).
#' @noRd
NULL
