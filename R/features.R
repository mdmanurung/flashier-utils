# Public conveniences around audited MIT singlecelljamboreeR helper adaptations.
#' Rank signed feature programs
#'
#' Uses singlecelljamboreeR rank_effects semantics (descending, average ties).
#' Positive/negative inputs preserve original estimates; display selection
#' filters sign explicitly. No probability or significance is implied.
#' @param fit Native fit or matched view.
#' @param factors Optional unique factor IDs.
#' @param direction absolute, positive or negative.
#' @param ties Only native average ties are supported.
#' @param representation Matched view or settings.
#' @return feature/factor/estimate/rank/direction/basis_id table.
#' @export
#' @author Peter Carbonetto; Matthew Stephens.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: ADAPTED / modified.
#'   Credits rank_effects; use provenance("rank_features")
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
#' head(rank_features(view,factors="F1",direction="absolute"))
rank_features <- function(fit, factors = NULL, direction = "absolute", ties = "average", representation = NULL) {
  if (ties != "average") stop("Only native average ties are supported")
  direction <- match.arg(direction, c("absolute", "positive", "negative"))
  view <- .resolve_view(fit, representation)
  if (!is.null(factors)) .require_biological_factors(view$manifest$factor_ids)
  B <- factor_effects(view, factors = factors)
  values <- switch(direction, absolute = abs(B), positive = B, negative = -B)
  ranks <- if (ncol(B)) .rank_effects(values) else B
  out <- .long_matrix(B, "feature")
  out$rank <- as.vector(ranks)
  out$direction <- rep(direction, nrow(out))
  out$basis_id <- rep(view$manifest$basis_id, nrow(out))
  attr(out, "analysis_metadata") <- view$manifest
  .attach_provenance(out, "rank_features", match.call(), .resolved_parameters("rank_features", environment()))
}

#' Compute least-extreme signed feature effects
#'
#' Adapted from singlecelljamboreeR compute_le_effects. The margin to the
#' most extreme same-direction competitor includes zero. By default all factors compete.
#' Distinctiveness depends on the recorded common scale and orientation.
#' @inheritParams rank_features
#' @param compare_factors Optional competition pool; NULL uses all factors.
#' @param engine singlecelljamboreeR selects the audited local implementation.
#' @param format long or matrix.
#' @return Original estimate and distinctiveness separately, or margin matrix.
#' @export
#' @author Peter Carbonetto; Matthew Stephens.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: ADAPTED / modified.
#'   Credits compute_le_effects; use provenance("factor_distinctiveness")
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
#' head(factor_distinctiveness(view))
factor_distinctiveness <- function(fit, representation = NULL, engine = "singlecelljamboreeR", direction = "both", format = "long", compare_factors = NULL) {
  if (engine != "singlecelljamboreeR") stop("Unsupported distinctiveness engine")
  direction <- match.arg(direction, c("both", "positive", "negative"))
  format <- match.arg(format, c("long", "matrix"))
  view <- .resolve_view(fit, representation)
  if (is.null(view$manifest$orientation) || is.null(view$manifest$scaling)) stop("Distinctiveness requires declared common scaling and orientation")
  B <- view$effects
  pool <- colnames(B)[.select_ids(compare_factors, colnames(B), "comparison factor")]
  margin <- .compute_le_effects(B)
  if (!is.null(compare_factors)) for (k in seq_len(ncol(B))) {
    competitors <- unique(c(colnames(B)[k], pool))
    margin[, k] <- .compute_le_effects(B[, competitors, drop = FALSE])[, 1L]
  }
  if (direction == "positive") margin[margin < 0] <- 0
  if (direction == "negative") margin[margin > 0] <- 0
  metadata <- c(view$manifest, list(competition_factors = pool,
    definition = "singlecelljamboreeR_least_extreme", code_reuse = "modified"))
  if (format == "matrix") {
    attr(margin, "analysis_metadata") <- metadata
    return(.attach_provenance(margin, "factor_distinctiveness", match.call(), .resolved_parameters("factor_distinctiveness", environment())))
  }
  out <- .long_matrix(B, "feature")
  out$distinctiveness <- as.vector(margin)
  out$direction <- rep(direction, nrow(out))
  out$basis_id <- rep(view$manifest$basis_id, nrow(out))
  attr(out, "analysis_metadata") <- metadata
  .attach_provenance(out, "factor_distinctiveness", match.call(), .resolved_parameters("factor_distinctiveness", environment()))
}

#' Select features by global cross-factor effect size
#'
#' Global selection for multi-factor heatmaps: rank each feature by
#' `max_k(abs(effect[gene, k]))` (or the signed extreme for
#' `direction = "positive"/"negative"`) over an explicitly selected factor
#' set, with stable ties (score, then feature ID). This is not the same as
#' [factor_features()], which selects per factor. Descriptive display
#' selection only; no inferential labeling.
#'
#' @param fit Native fit or representation.
#' @param factors Explicit factor set to compete features across.
#' @param n Maximum number of features.
#' @param direction absolute, positive or extreme-negative.
#' @param representation Optional representation name.
#' @return Data frame with feature, argmax factor, signed estimate at the
#'   argmax, selection score, rank, direction and basis_id.
#' @author Mikhael Manurung. Conceptual cross-factor display workflow: ImmGenT-GP.
#' @export
factor_features_global <- function(fit, factors, n = 20L, direction = c("absolute", "positive", "negative"), representation = NULL) {
  direction <- match.arg(direction)
  view <- .resolve_view(fit, representation)
  .require_biological_factors(view$manifest$factor_ids)
  fids <- view$manifest$factor_ids[.select_ids(factors, view$manifest$factor_ids, "factor")]
  if (!length(fids)) stop("Select at least one factor")
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) || n < 0 || n != as.integer(n)) stop("n must be a nonnegative integer")
  B <- factor_effects(view, factors = fids)
  raw <- switch(direction, absolute = apply(abs(B), 1L, max), positive = apply(B, 1L, max), negative = apply(B, 1L, min))
  eligible <- switch(direction, absolute = raw != 0, positive = raw > 0, negative = raw < 0)
  score <- if (direction == "negative") -raw else raw
  argmax <- if (direction == "negative") apply(B, 1L, which.min) else apply(abs(B), 1L, which.max)
  if (direction == "positive") argmax <- apply(B, 1L, which.max)
  ord <- order(-score, rownames(B))
  keep <- ord[eligible[ord]]
  keep <- utils::head(keep, n)
  out <- data.frame(feature = rownames(B)[keep],
                    factor = colnames(B)[argmax[keep]],
                    estimate = B[cbind(keep, argmax[keep])],
                    score = score[keep],
                    rank = seq_along(keep),
                    direction = rep(direction, length(keep)),
                    basis_id = rep(view$manifest$basis_id, length(keep)),
                    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  attr(out, "analysis_metadata") <- c(view$manifest, list(selection = "global_cross_factor_max", display_tie_rule = "score_then_feature_ID", selection_factors = colnames(B), argmax_tie_rule = "first_factor_in_selection_order"))
  .attach_provenance(out, "factor_features_global", match.call(), .resolved_parameters("factor_features_global", environment()))
}

#' Select largest and distinctive program features
#'
#' Returns two explicit ranking views when select is both. Cutoff ties are
#' ordered by feature ID, with signed nonzero eligibility. No biological labels.
#' @inheritParams rank_features
#' @param factor One unique factor ID.
#'   A fit with no biological factors cannot supply a feature selection.
#' @param select largest, distinctive or both.
#' @param compare_factors Optional competition pool; NULL uses all factors.
#' @param n Nonnegative maximum number per selection.
#' @return Tidy descriptive table with selection and tie rule metadata.
#' @export
#' @author Peter Carbonetto; Matthew Stephens.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: NEW_QOL / none.
#'   Credits rank_effects; compute_le_effects; use provenance("factor_features")
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
#' factor_features(view,"F1",select="both",n=3)
factor_features <- function(fit, factor, select = "both", n = 20L, direction = "absolute", representation = NULL, compare_factors = NULL) {
  select <- match.arg(select, c("largest", "distinctive", "both"))
  direction <- match.arg(direction, c("absolute", "positive", "negative"))
  if (length(factor) != 1L) stop("Select exactly one factor")
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) || n < 0 || n != as.integer(n)) stop("n must be a nonnegative integer")
  view <- .resolve_view(fit, representation)
  .require_biological_factors(view$manifest$factor_ids)
  .select_ids(factor, view$manifest$factor_ids, "factor")
  selections <- if (select == "both") c("largest", "distinctive") else select
  parts <- lapply(selections, function(selection) {
    out <- if (selection == "largest") rank_features(view, factors = factor, direction = direction) else {
      all <- factor_distinctiveness(view, compare_factors = compare_factors)
      all <- all[all$factor == factor, , drop = FALSE]
      score <- switch(direction, absolute = abs(all$distinctiveness), positive = all$distinctiveness, negative = -all$distinctiveness)
      all$rank <- rank(-score, ties.method = "average")
      all$direction <- rep(direction, nrow(all))
      all
    }
    value <- if (selection == "largest") out$estimate else out$distinctiveness
    eligible <- switch(direction, absolute = value != 0, positive = value > 0, negative = value < 0)
    out <- out[eligible, , drop = FALSE]
    out <- utils::head(out[order(out$rank, out$feature), , drop = FALSE], n)
    out$selection <- rep(selection, nrow(out))
    if (!"distinctiveness" %in% names(out)) out$distinctiveness <- rep(NA_real_, nrow(out))
    out[, c("feature", "factor", "estimate", "rank", "selection", "direction", "basis_id", "distinctiveness")]
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  attr(out, "analysis_metadata") <- c(view$manifest, list(display_tie_rule = "rank_then_feature_ID", competition_factors = colnames(view$effects)[.select_ids(compare_factors, colnames(view$effects), "comparison factor")]))
  .attach_provenance(out, "factor_features", match.call(), .resolved_parameters("factor_features", environment()))
}
