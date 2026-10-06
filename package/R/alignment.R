# Public psych congruence and clue assignment orchestration; no new solver.
.align_programs <- function(reference, target, min_congruence = NULL, allow_sign_flip = NULL,
                            min_shared_features = 2L, feature_map = NULL) {
  reference <- .validate_view(reference); target <- .validate_view(target)
  if (identical(reference$manifest$feature_scale, "unknown") || !identical(reference$manifest$feature_scale, target$manifest$feature_scale)) stop("Program comparison requires matching declared feature units")
  R <- reference$effects; T <- target$effects
  if (!is.null(feature_map)) {
    if (!is.character(feature_map)) stop("feature_map must be an explicit one-to-one character mapping")
    .validate_ids(names(feature_map), "Reference feature map")
    .validate_ids(unname(feature_map), "Target feature map")
    index <- which(rownames(T) %in% feature_map)
    T <- T[index, , drop = FALSE]
    rownames(T) <- names(feature_map)[match(rownames(T), feature_map)]
  }
  features <- intersect(rownames(R), rownames(T))
  if (length(min_shared_features) != 1L || !is.finite(min_shared_features) || min_shared_features < 1) stop("min_shared_features must be positive")
  if (length(features) < min_shared_features) stop("Insufficient shared named features for program comparison")
  possible_flip <- all(reference$manifest$constraints == "signed") && all(target$manifest$constraints == "signed")
  if (is.null(allow_sign_flip)) allow_sign_flip <- possible_flip
  if (!is.logical(allow_sign_flip) || length(allow_sign_flip) != 1L || is.na(allow_sign_flip)) stop("allow_sign_flip must be TRUE or FALSE")
  if (allow_sign_flip && !possible_flip) stop("Sign flips violate constrained or unknown prior support")
  if (!is.null(min_congruence) && (length(min_congruence) != 1L || !is.finite(min_congruence) || min_congruence < 0 || min_congruence > 1)) stop("min_congruence must lie between zero and one")
  .check_engine("psych", "factor.congruence"); .check_engine("clue", "solve_LSAP")
  rids <- colnames(R); tids <- colnames(T)
  if (is.null(rids)) rids <- character()
  if (is.null(tids)) tids <- character()
  R <- R[features, , drop = FALSE]; T <- T[features, , drop = FALSE]
  similarity <- congruence <- matrix(NA_real_, length(rids), length(tids), dimnames = list(rids, tids))
  valid_R <- which(colSums(R^2) > 0); valid_T <- which(colSums(T^2) > 0)
  if (length(valid_R) && length(valid_T)) {
    native <- psych::factor.congruence(R[, valid_R, drop = FALSE], T[, valid_T, drop = FALSE], digits = 15)
    congruence[valid_R, valid_T] <- native
  }
  similarity <- if (allow_sign_flip) abs(congruence) else congruence
  threshold <- if (is.null(min_congruence)) 0 else min_congruence
  benefits <- pmax(similarity - threshold, 0); benefits[!is.finite(benefits)] <- 0
  assignment <- rep(NA_integer_, length(rids)); margin <- NA_real_
  if (length(rids)) {
    augmented <- cbind(benefits, matrix(0, length(rids), length(rids)))
    chosen <- as.integer(clue::solve_LSAP(augmented, maximum = TRUE))
    real <- which(chosen <= length(tids))
    for (i in real) if (is.finite(similarity[i, chosen[i]]) && benefits[i, chosen[i]] > 0) assignment[i] <- chosen[i]
    selected <- which(!is.na(assignment))
    if (length(selected)) {
      optimum <- sum(augmented[cbind(seq_along(rids), chosen)])
      # ponytail: K native reassignments for the ambiguity margin; optimize for very large K.
      losses <- vapply(selected, function(i) {
        alternative <- augmented; alternative[i, assignment[i]] <- 0
        next_choice <- as.integer(clue::solve_LSAP(alternative, maximum = TRUE))
        optimum - sum(alternative[cbind(seq_along(rids), next_choice)])
      }, numeric(1))
      margin <- min(losses)
    }
  }
  matches <- data.frame(reference_factor = rids, target_factor = ifelse(is.na(assignment), NA_character_, tids[assignment]),
                        congruence = NA_real_, similarity = NA_real_, sign = NA_real_, matched = rep(if (is.null(min_congruence)) NA else FALSE, length(rids)), candidate = !is.na(assignment),
                        shared_features = rep(length(features), length(rids)), eligible = seq_along(rids) %in% valid_R & length(valid_T) > 0L, status = rep("unmatched_reference", length(rids)))
  for (i in which(!is.na(assignment))) {
    matches$congruence[i] <- congruence[i, assignment[i]]; matches$similarity[i] <- similarity[i, assignment[i]]
    matches$sign[i] <- if (allow_sign_flip && matches$congruence[i] < 0) -1 else 1
    matches$matched[i] <- if (is.null(min_congruence)) NA else TRUE
    matches$status[i] <- if (is.null(min_congruence)) "candidate_only" else "accepted"
  }
  unmatched <- setdiff(tids, matches$target_factor[!is.na(matches$target_factor)])
  if (length(unmatched)) matches <- rbind(matches, data.frame(reference_factor = NA_character_, target_factor = unmatched, congruence = NA_real_, similarity = NA_real_, sign = NA_real_, matched = if (is.null(min_congruence)) NA else FALSE, candidate = FALSE, shared_features = length(features), eligible = FALSE, status = "unmatched_target"))
  ties <- vapply(seq_along(rids), function(i) {
    values <- sort(similarity[i, is.finite(similarity[i, ])], decreasing = TRUE)
    if (length(values) < 2L) NA_real_ else values[1] - values[2]
  }, numeric(1))
  list(matches = matches, similarity_matrix = similarity, congruence_matrix = congruence,
       unmatched_reference = matches$reference_factor[!is.na(matches$reference_factor) & is.na(matches$target_factor)],
       unmatched_target = matches$target_factor[is.na(matches$reference_factor) & !is.na(matches$target_factor)],
       feature_coverage = data.frame(side = c("reference", "target"), total_features = c(nrow(reference$effects),nrow(target$effects)), shared_features = length(features), fraction = length(features)/c(nrow(reference$effects),nrow(target$effects))),
       analysis_metadata = list(reference_basis = reference$manifest, target_basis = target$manifest, shared_feature_ids = features,
         reference_feature_count = nrow(reference$effects), target_feature_count = nrow(target$effects),
         reference_coverage = length(features) / nrow(reference$effects), target_coverage = length(features) / nrow(target$effects),
         dropped_reference_features = setdiff(rownames(reference$effects), features), dropped_target_features = setdiff(rownames(target$effects), if (is.null(feature_map)) features else unname(feature_map[features])),
         min_congruence = min_congruence, allow_sign_flip = allow_sign_flip, solver = "clue::solve_LSAP", congruence_digits = 15L,
         tie_rule = "native_solver_in_declared_factor_order", row_top_two_margin = stats::setNames(ties, rids), global_assignment_margin = margin,
         margin_definition = "minimum_objective_loss_excluding_a_selected_real_edge", interpretation = "concordance_not_measurement_invariance"))
}

#' Compare independently fitted feature programs
#'
#' Tucker congruence on named shared features, with native rectangular assignment
#' and dummy unmatched nodes. Without a threshold all matches are candidates.
#' Correspondence does not establish measurement invariance or activity units.
#' @param reference,target Native fits with recorded side, or matched views.
#' @param feature_map Optional one-to-one reference-to-target feature mapping.
#' @param alignment Only congruence is supported.
#' @param min_congruence Optional acceptance threshold; NULL yields candidate-only.
#' @param min_shared_features Computational minimum overlap (not scientific sufficiency).
#' @param allow_sign_flip Optional declaration; signed support on both sides required.
#' @param bootstrap Donor-refit bootstrap is unsupported here; use explicit stability recipes.
#' @param refit_inputs Reserved for a separately validated biological bootstrap recipe.
#' @param control Reserved; unknown settings error.
#' @return A replicability_result containing matches, full similarities and coverage metadata.
#' @export
#' @details Origin: NEW_ORCHESTRATION / call_only.
#'   Credits psych::factor.congruence; clue::solve_LSAP (probe pending); use provenance("factor_replicability")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examplesIf requireNamespace("psych",quietly=TRUE) && requireNamespace("clue",quietly=TRUE)
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' factor_replicability(view,view,min_congruence=0.8)$matches
factor_replicability <- function(reference, target, feature_map = NULL, alignment = "congruence",
                                 min_congruence = NULL, min_shared_features = 2L,
                                 allow_sign_flip = NULL, bootstrap = FALSE, refit_inputs = NULL, control = list()) {
  if (alignment != "congruence" || length(control)) stop("Only declared native congruence alignment is supported")
  if (bootstrap) stop("Biological bootstrap requires explicit donor refit inputs and a separate validated recipe; feature bootstrap is not supported")
  result <- .align_programs(.resolve_view(reference), .resolve_view(target), min_congruence, allow_sign_flip, min_shared_features, feature_map)
  .attach_provenance(structure(c(list(schema_version = "1.0.0"), result), class = "replicability_result"), "factor_replicability", match.call(), .resolved_parameters("factor_replicability", environment()))
}
