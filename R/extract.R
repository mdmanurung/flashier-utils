# Newly written selections in the shared matched basis.
.extract <- function(fit, side, sample_side = NULL, representation = NULL,
                     factors = NULL, ids = NULL, format = "matrix") {
  view <- .resolve_view(fit, representation, sample_side)
  format <- match.arg(format, c("matrix", "long"))
  mat <- view[[side]]
  rows <- .select_ids(ids, rownames(mat), side)
  columns <- .select_ids(factors, view$manifest$factor_ids, "factor")
  mat <- mat[rows, columns, drop = FALSE]
  if (format == "matrix") {
    attr(mat, "analysis_metadata") <- view$manifest
    return(mat)
  }
  out <- .long_matrix(mat, if (side == "activity") "sample" else "feature")
  out$basis_id <- rep(view$manifest$basis_id, nrow(out))
  attr(out, "analysis_metadata") <- view$manifest
  out
}

.long_matrix <- function(mat, id_col, value_col = "estimate") {
  out <- data.frame(id = rep(rownames(mat), times = ncol(mat)),
                    factor = rep(if (ncol(mat)) colnames(mat) else character(), each = nrow(mat)), value = as.vector(mat))
  names(out) <- c(id_col, "factor", value_col)
  out
}

#' Extract sample by factor activity
#' @param fit Native flash fit or matched representation.
#' @param sample_side Explicit native sample side if not recorded.
#' @param representation Matched view or named standardization settings.
#' @param factors Optional unique factor IDs in requested order.
#' @param samples Optional unique sample IDs in requested order.
#' @param format matrix or long; selection precedes reshaping.
#' @return Matrix or sample/factor/estimate table carrying basis metadata.
#' @export
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: NEW_QOL / none.
#'   Credits native means / ldf; use provenance("factor_activity")
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
#' factor_activity(view,factors="F1",samples=c("s2","s1"))
factor_activity <- function(fit, sample_side = NULL, representation = NULL,
                            factors = NULL, samples = NULL, format = "matrix") {
  .attach_provenance(.extract(fit, "activity", sample_side, representation, factors, samples, format), "factor_activity", match.call(), .resolved_parameters("factor_activity", environment()))
}

#' Extract feature by factor effects
#' @inheritParams factor_activity
#' @param features Optional unique feature IDs in requested order.
#' @return Matrix or feature/factor/estimate table carrying basis metadata.
#' @export
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: NEW_QOL / none.
#'   Credits native means / ldf; use provenance("factor_effects")
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
#' factor_effects(view,factors="F1",features=c("g2","g1"))
factor_effects <- function(fit, sample_side = NULL, representation = NULL,
                           factors = NULL, features = NULL, format = "matrix") {
  .attach_provenance(.extract(fit, "effects", sample_side, representation, factors, features, format), "factor_effects", match.call(), .resolved_parameters("factor_effects", environment()))
}

#' Extract native posterior standard deviations in a matched basis
#'
#' Native variational uncertainty is conditional on estimated priors. Posterior
#' SD is not a downstream regression standard error. Missing moments are NA.
#' @inheritParams factor_activity
#' @param side activity or effects.
#' @param ids Optional sample or feature IDs.
#' @param include_lfsr Include native posterior false sign summaries if available.
#' @param require_uncertainty Error when moments are unavailable.
#' @return Table with estimate, posterior.sd, uncertainty_type and basis_id.
#' @export
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: WRAPPER / call_only.
#'   Credits native posterior moments / public accessors; use provenance("factor_uncertainty")
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
#' factor_uncertainty(fit,"activity",factors="F1",ids=c("s1","s2"))
factor_uncertainty <- function(fit, side = c("activity", "effects"), representation = NULL,
                               sample_side = NULL, factors = NULL, ids = NULL,
                               include_lfsr = TRUE, require_uncertainty = FALSE) {
  side <- match.arg(side)
  if (!inherits(fit, "flash")) stop("Posterior uncertainty requires the native fit")
  view <- .resolve_view(fit, representation, sample_side)
  selected <- .extract(view, side, factors = factors, ids = ids, format = "long")
  prefix <- if ((side == "activity") == (view$manifest$sample_side == "rows")) "L" else "F"
  multiplier <- view$manifest[[paste0(side, "_multiplier")]]
  psd <- fit[[paste0(prefix, "_psd")]]
  id_col <- if (side == "activity") "sample" else "feature"
  lookup <- function(mat, transform = FALSE) {
    if (is.null(mat)) return(rep(NA_real_, nrow(selected)))
    if (!identical(dim(mat), dim(view[[side]]))) stop("Native uncertainty dimensions do not match basis")
    if (transform) mat <- sweep(mat, 2L, abs(multiplier), "*")
    index <- cbind(match(selected[[id_col]], rownames(view[[side]])), match(selected$factor, view$manifest$factor_ids))
    mat[index]
  }
  selected$posterior.sd <- lookup(psd, TRUE)
  if (require_uncertainty && anyNA(selected$posterior.sd)) stop("Native posterior uncertainty is unavailable")
  if (include_lfsr && !is.null(fit[[paste0(prefix, "_lfsr")]])) selected$lfsr <- lookup(fit[[paste0(prefix, "_lfsr")]])
  selected$uncertainty_type <- rep("conditional_variational_ebmf", nrow(selected))
  selected$status <- ifelse(is.na(selected$posterior.sd), "unavailable", "supported")
  selected$reason <- ifelse(is.na(selected$posterior.sd), "native moments unavailable", NA_character_)
  .attach_provenance(selected, "factor_uncertainty", match.call(), .resolved_parameters("factor_uncertainty", environment()))
}

#' Extract native factor proportions of variance explained
#'
#' Uses native second-moment and residual-variance semantics. This is not
#' unique variance decomposition or held-out prediction R squared.
#' @param fit Native flash fit or view carrying native PVE.
#' @param factors Optional unique factor IDs in requested order.
#' @return factor/pve/definition table.
#' @export
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: WRAPPER / none.
#'   Credits native pve field; use provenance("factor_pve")
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
#' factor_pve(fit)
factor_pve <- function(fit, factors = NULL) {
  if (inherits(fit, "factor_representation")) {
    view <- .validate_view(fit)
    ids <- view$manifest$factor_ids
    pve <- view$manifest$native_pve
  } else {
    if (!inherits(fit, "flash")) stop("Expected a native flash fit")
    ids <- attr(fit, "flashbridge_metadata")$factor_ids
    if (is.null(ids)) ids <- if (fit$n_factors == 0L) character() else paste0("F", seq_len(fit$n_factors))
    pve <- fit$pve
  }
  if (length(ids) && (is.null(pve) || length(pve) != length(ids))) stop("Native PVE is unavailable")
  index <- .select_ids(factors, ids, "factor")
  .attach_provenance(data.frame(factor = ids[index], pve = as.numeric(pve[index]), definition = rep("flashier_native", length(index))), "factor_pve", match.call(), .resolved_parameters("factor_pve", environment()))
}
