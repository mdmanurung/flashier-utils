# Thin wrappers around public flashier operations; no private slot edits.
.factor_ids <- function(fit) {
  if (!inherits(fit, "flash")) stop("Expected a native flash fit")
  ids <- attr(fit, "flashbridge_metadata")$factor_ids
  if (is.null(ids)) ids <- if (fit$n_factors == 0) character() else paste0("F", seq_len(fit$n_factors))
  ids
}

.native_changed <- function(original, updated, ids, operation) {
  metadata <- attr(original, "flashbridge_metadata")
  if (is.null(metadata)) metadata <- list()
  metadata$factor_ids <- ids
  metadata$fit_id <- if (operation == "reorder" && !is.null(metadata$fit_id)) metadata$fit_id else
    digest::digest(list(.raw_pair(updated), Sys.time(), operation), algo = "sha256")
  metadata$parent_fit_id <- attr(original, "flashbridge_metadata")$fit_id
  metadata$operation <- operation
  attr(updated, "flashbridge_metadata") <- metadata
  updated
}

#' Fit factors constrained to be nonnegative
#'
#' Native point-exponential priors on both sides, with finite greedy Kmax.
#' No cycling or NNLM initialization is advertised. Observed noisy values may
#' be negative: the constraint concerns latent factors, not observations.
#' @inheritParams fit_ebmf
#' @param initialization Only native_greedy is supported.
#' @param ... Other explicit fit_ebmf arguments, excluding ebnm_fn.
#' @return Native flash fit; attained rank may be lower than max_factors.
#' @export
#' @author Peter Carbonetto; Matthew Stephens.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: ADAPTED / call_only.
#'   Credits flashier_nmf workflow; bounded native recipe; use provenance("fit_nonnegative_ebmf")
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
#' nonnegative <- fit_nonnegative_ebmf(abs(X),"rows",max_factors=2,seed=4,verbose=0)
#' factor_activity(nonnegative)
fit_nonnegative_ebmf <- function(X, sample_side, max_factors = 50L, initialization = "native_greedy",
                                seed = NULL, control = list(), ...) {
  if (initialization != "native_greedy") stop("Only bounded native_greedy initialization is supported; NNLM/cycling is unavailable")
  args <- list(...)
  if ("ebnm_fn" %in% names(args)) stop("Nonnegative fitting sets native point-exponential priors")
  .attach_provenance(do.call(fit_ebmf, c(list(X = X, sample_side = sample_side, max_factors = max_factors,
    ebnm_fn = ebnm::ebnm_point_exponential, seed = seed, control = control), args)), "fit_nonnegative_ebmf", match.call(), .resolved_parameters("fit_nonnegative_ebmf", environment()))
}

#' Backfit named existing factors
#' @param fit Native fit with retained native data.
#' @param factors Optional factor IDs; defaults to all.
#' @param control Named public flash_backfit settings.
#' @param seed Optional scoped seed.
#' @return New native fit with retained factor identity and new numerical identity.
#' @export
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: WRAPPER / call_only.
#'   Credits flash_backfit; use provenance("refit_factors")
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
#' refitted <- refit_factors(fit,control=list(maxiter=10,verbose=0),seed=4)
#' factor_pve(refitted)
refit_factors <- function(fit, factors = NULL, control = list(), seed = NULL) {
  ids <- .factor_ids(fit)
  kset <- .select_ids(factors, ids, "factor")
  if (!length(kset)) stop("No factors to refit")
  if (!is.list(control) || (length(control) && (is.null(names(control)) || anyDuplicated(names(control)))) || any(!names(control) %in% setdiff(names(formals(flashier::flash_backfit)), c("flash", "kset")))) stop("Unknown or protected backfit control")
  result <- .with_seed(seed, do.call(flashier::flash_backfit, c(list(flash = fit, kset = kset), control)))
  .attach_provenance(.native_changed(fit, result, ids, "backfit"), "refit_factors", match.call(), .resolved_parameters("refit_factors", environment()))
}

#' Fix or unfix a named native factor side
#' @param fit Native fit with recorded sample side.
#' @param factors Required factor IDs.
#' @param side activity or effects.
#' @param ids Optional sample/feature IDs for partial fixing.
#' @param use_fixed_in_ebnm Native handling of fixed values in prior estimation.
#' @param action fix or unfix; native unfix operates on the selected pair.
#' @return New native fit retaining named factor identity.
#' @export
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: WRAPPER / call_only.
#'   Credits flash_factors_fix; flash_factors_unfix; use provenance("fix_factors")
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
#' fixed <- fix_factors(fit,"F1","effects")
#' factor_effects(fixed,factors="F1",features="g1")
fix_factors <- function(fit, factors, side, ids = NULL, use_fixed_in_ebnm = NULL, action = "fix") {
  action <- match.arg(action, c("fix", "unfix"))
  side <- match.arg(side, c("activity", "effects"))
  all_ids <- .factor_ids(fit)
  kset <- .select_ids(factors, all_ids, "factor")
  if (!length(kset)) stop("No factors selected")
  sample_side <- .resolve_sample_side(fit)
  mode <- if ((side == "activity") == (sample_side == "rows")) "loadings" else "factors"
  if (action == "unfix") {
    if (!is.null(ids) || !is.null(use_fixed_in_ebnm)) stop("Native unfix applies to the entire selected factor pair; partial IDs unsupported")
    result <- flashier::flash_factors_unfix(fit, kset = kset)
  } else {
    view <- standardize_factors(fit, scaling = "raw")
    fixed_idx <- if (is.null(ids)) NULL else .select_ids(ids, rownames(view[[side]]), side)
    result <- flashier::flash_factors_fix(fit, kset = kset, which_dim = mode, fixed_idx = fixed_idx, use_fixed_in_ebnm = use_fixed_in_ebnm)
  }
  .attach_provenance(.native_changed(fit, result, all_ids, action), "fix_factors", match.call(), .resolved_parameters("fix_factors", environment()))
}

#' Remove named factors without refitting
#' @param fit Native flash fit.
#' @param factors Unique factor IDs to remove.
#' @return Native fit retaining remaining factor IDs.
#' @export
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: WRAPPER / call_only.
#'   Credits flash_factors_remove; use provenance("remove_factors")
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
#' remaining <- remove_factors(fit,"F1")
#' colnames(factor_activity(remaining))
remove_factors <- function(fit, factors) {
  ids <- .factor_ids(fit)
  kset <- .select_ids(factors, ids, "factor")
  if (!length(kset)) stop("No factors selected")
  .attach_provenance(.native_changed(fit, flashier::flash_factors_remove(fit, kset), ids[-kset], "remove"), "remove_factors", match.call(), .resolved_parameters("remove_factors", environment()))
}

#' Reorder native factors by exact named permutation
#' @param fit Native flash fit.
#' @param order Exact permutation of retained factor IDs.
#' @return Native fit with stable identities in requested order.
#' @export
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: WRAPPER / call_only.
#'   Credits flash_factors_reorder; use provenance("reorder_factors")
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
#' reordered <- reorder_factors(fit,rev(view$manifest$factor_ids))
#' colnames(factor_activity(reordered))
reorder_factors <- function(fit, order) {
  ids <- .factor_ids(fit)
  kset <- .select_ids(order, ids, "factor")
  if (length(kset) != length(ids)) stop("order must be an exact permutation of all factors")
  .attach_provenance(.native_changed(fit, flashier::flash_factors_reorder(fit, kset), order, "reorder"), "reorder_factors", match.call(), .resolved_parameters("reorder_factors", environment()))
}
