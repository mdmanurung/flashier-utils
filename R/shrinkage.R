# Public mashr fit and public ashr summary generics, registered for class mash.
#' Shrink comparable factor contrasts across conditions
#'
#' Requires a common activity basis and all eligible factor rows, never only
#' significant rows. V describes correlated estimation errors across conditions,
#' not true-effect covariance. Small factor counts limit covariance learning.
#' @param results A coefficient table with attached activity manifest.
#' @param condition Condition column name.
#' @param method Only mash is supported.
#' @param term Required comparable term value.
#' @param covariances canonical or a caller-supplied covariance list.
#' @param V Optional condition error-correlation matrix.
#' @param independent_conditions Explicit TRUE when using identity V.
#' @param control Named public mash settings; outputlevel must retain summaries.
#' @return Posterior table, native fit, input matrices and covariance metadata.
#' @export
#' @author Mikhael Manurung; Matthew Stephens; Sarah Urbut; Gao Wang; Yuxin Zou; Peter Carbonetto (mashr); ashr contributors.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: WRAPPER / call_only.
#'   Credits mashr public API (implementation audit pending); use provenance("shrink_factor_effects")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examplesIf requireNamespace("mashr",quietly=TRUE)
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' # Small-row illustration; does not establish covariance/lfsr calibration.
#' ids <- view$manifest$factor_ids
#' coefficients <- data.frame(factor=rep(ids,2),condition=rep(c("A","B"),each=length(ids)),
#'   term="treatment",estimate=seq_len(2*length(ids))/5,std.error=0.5,basis_id=view$manifest$basis_id)
#' attr(coefficients,"analysis_metadata") <- view$manifest
#' shrink_factor_effects(coefficients,"condition",term="treatment",independent_conditions=TRUE)$table
shrink_factor_effects <- function(results, condition, method = "mash", term, covariances = "canonical",
                                  V = NULL, independent_conditions = NULL, control = list()) {
  if (method != "mash" || missing(term)) stop("Choose mash and one comparable term explicitly")
  .check_engine("mashr", c("mash_set_data", "cov_canonical", "mash"))
  .check_engine("ashr", c("get_pm", "get_psd", "get_lfsr"))
  if (is.list(results) && !is.data.frame(results) && !is.null(results$table)) {
    manifest <- results$analysis_metadata
    table <- results$table
  } else { table <- results; manifest <- attr(table, "analysis_metadata") }
  required <- c("factor", "estimate", "std.error", "term", "basis_id", condition)
  if (!is.data.frame(table) || !all(required %in% names(table)) || !is.list(manifest)) stop("Supply factor coefficient table and its activity basis manifest")
  table <- table[table$term %in% term, , drop = FALSE]
  if (length(term) != 1L || !nrow(table) || length(unique(table$term)) != 1L) stop("Exactly one comparable term is required")
  if (length(unique(table$basis_id)) != 1L || !identical(unique(table$basis_id), manifest$basis_id)) stop("Condition effects have incompatible activity bases")
  ids <- manifest$factor_ids
  conditions <- unique(as.character(table[[condition]]))
  .validate_ids(conditions, "Condition")
  if (length(conditions) < 2L) stop("mash requires at least two conditions")
  if (any(!is.finite(table$estimate)) || any(!is.finite(table$std.error)) || any(table$std.error <= 0)) stop("mash requires finite estimates and positive SEs")
  if (anyDuplicated(paste(table$factor, table[[condition]], sep = "\r"))) stop("Duplicate factor/condition coefficients")
  if (!setequal(unique(table$factor), ids) || nrow(table) != length(ids)*length(conditions)) stop("Supply every eligible factor in every condition; selected significant rows or missing cells are unsupported")
  Bhat <- Shat <- matrix(NA_real_, length(ids), length(conditions), dimnames = list(ids, conditions))
  index <- cbind(match(table$factor, ids), match(table[[condition]], conditions))
  Bhat[index] <- table$estimate; Shat[index] <- table$std.error
  if (is.null(V)) {
    if (!isTRUE(independent_conditions)) stop("Supply V or explicitly independent_conditions=TRUE")
    V <- diag(length(conditions)); dimnames(V) <- list(conditions, conditions)
  } else {
    if (!is.matrix(V) || !is.numeric(V) || any(!is.finite(V)) || !identical(dim(V), c(length(conditions),length(conditions)))) stop("V must be a finite named condition correlation matrix")
    V <- V[.select_ids(conditions, rownames(V), "V condition"), .select_ids(conditions, colnames(V), "V condition"), drop = FALSE]
    if (max(abs(V-t(V))) > 1e-10 || max(abs(diag(V)-1)) > 1e-10 || min(eigen(V,symmetric=TRUE,only.values=TRUE)$values) <= 0) stop("V must be symmetric positive definite with unit diagonal")
  }
  if (!is.list(control) || (length(control) && is.null(names(control))) || any(!names(control) %in% setdiff(names(formals(mashr::mash)), c("data", "Ulist", "...")))) stop("Unsupported native mash settings")
  if (!is.null(control$outputlevel) && control$outputlevel < 2L) stop("mash outputlevel must retain posterior summaries")
  data <- mashr::mash_set_data(Bhat, Shat, V = V)
  U <- if (identical(covariances, "canonical")) mashr::cov_canonical(data) else covariances
  if (!is.list(U) || !length(U)) stop("covariances must be canonical or an explicit native covariance list")
  native <- do.call(mashr::mash, c(list(data = data, Ulist = U), if (!"verbose" %in% names(control)) list(verbose = FALSE) else list(), control))
  pm <- ashr::get_pm(native); psd <- ashr::get_psd(native); lfsr <- ashr::get_lfsr(native)
  out <- data.frame(factor = rep(ids, length(conditions)), condition = rep(conditions, each = length(ids)), term = term,
    estimate = as.vector(pm), posterior.sd = as.vector(psd), lfsr = as.vector(lfsr), basis_id = manifest$basis_id)
  metadata <- c(manifest, list(engine = "mashr", engine_version = as.character(utils::packageVersion("mashr")), error_correlation = V,
    independent_conditions = independent_conditions, eligible_factor_ids = ids, covariance_family = if (is.character(covariances)) covariances else "caller_supplied",
    control = control, small_row_caveat = "few factors and dependent rows can limit empirical Bayes learning; sensitivity required", uncertainty_type = "mash_posterior_conditional_on_first_stage"))
  attr(out, "analysis_metadata") <- metadata
  .attach_provenance(list(table = out, native = native, Bhat = Bhat, Shat = Shat, covariances = U, analysis_metadata = metadata), "shrink_factor_effects", match.call(), .resolved_parameters("shrink_factor_effects", environment()))
}
