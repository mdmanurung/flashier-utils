# Thin wrappers around public flashier operations; no private slot edits.
.factor_ids <- function(fit) {
  if (!inherits(fit, "flash")) stop("Expected a native flash fit")
  ids <- attr(fit, "flashbridge_metadata")$factor_ids
  if (is.null(ids)) ids <- if (fit$n_factors == 0) character() else paste0("F", seq_len(fit$n_factors))
  ids
}

.native_factor_indices <- function(fit) {
  metadata <- attr(fit, "flashbridge_metadata")
  indices <- metadata$native_indices
  if (is.null(indices)) indices <- seq_along(.factor_ids(fit))
  if (!is.numeric(indices) || length(indices) != length(.factor_ids(fit)) || anyNA(indices) ||
      any(indices < 1 | indices > fit$n_factors | indices != floor(indices)) || anyDuplicated(indices)) {
    stop("Invalid native factor mapping")
  }
  if (!is.null(metadata$native_factor_ids)) {
    .validate_ids(metadata$native_factor_ids, "Native factor")
    if (!is.null(metadata$intercept_index) &&
        (!identical(metadata$intercept_index, 1L) ||
         !identical(metadata$native_factor_ids[metadata$intercept_index], "Intercept") ||
         !identical(unname(metadata$factor_roles[metadata$intercept_index]), "offset"))) {
      stop("Invalid protected intercept role")
    }
    if (length(metadata$native_factor_ids) != fit$n_factors ||
        !identical(metadata$native_factor_ids[indices], .factor_ids(fit)) ||
        !identical(names(metadata$factor_roles), metadata$native_factor_ids) ||
        anyNA(metadata$factor_roles) ||
        any(!metadata$factor_roles %in% c("offset", "factor", "initialized")) ||
        any(metadata$factor_roles[indices] == "offset") ||
        !setequal(indices, setdiff(seq_len(fit$n_factors), metadata$intercept_index)) ||
        !setequal(which(metadata$factor_roles == "offset"), metadata$intercept_index)) {
      stop("Invalid native factor identity/role mapping")
    }
  }
  if (!is.null(metadata$factor_support) || !is.null(metadata$support_version)) {
    .factor_support(list(factor_ids = .factor_ids(fit), native_indices = indices,
      factor_roles = metadata$factor_roles, constraints = metadata$prior_support,
      support_version = metadata$support_version, factor_support = metadata$factor_support))
  }
  indices
}

.native_changed <- function(original, updated, ids, operation) {
  .native_factor_indices(original)
  metadata <- attr(original, "flashbridge_metadata")
  if (is.null(metadata)) metadata <- list()
  metadata$factor_ids <- ids
  if (!is.null(metadata$native_indices)) {
    metadata$native_indices <- setdiff(seq_len(updated$n_factors), metadata$intercept_index)
  }
  if (!is.null(metadata$native_factor_ids)) {
    metadata$native_factor_ids <- c(if (!is.null(metadata$intercept_index)) "Intercept", ids)
    metadata$factor_roles <- metadata$factor_roles[metadata$native_factor_ids]
  }
  if (!is.null(metadata$factor_support)) {
    records <- metadata$factor_support[match(ids, metadata$factor_support$factor), , drop = FALSE]
    rownames(records) <- NULL
    records$native_index <- metadata$native_indices
    metadata$factor_support <- records
    metadata$prior_support <- .support_summary(records, metadata$prior_support)
  }
  # every edit, including reorder, gets a new identity so views of the original fit are refused
  metadata$fit_id <- digest::digest(list(.raw_pair(updated), Sys.time(), operation), algo = "sha256")
  metadata$parent_fit_id <- attr(original, "flashbridge_metadata")$fit_id
  metadata$operation <- operation
  if (operation != "reorder") metadata$refit_policy <- "callback_required"
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
#' @keywords internal
#' @author Peter Carbonetto; Matthew Stephens.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: ADAPTED / call_only.
#'   Credits flashier_nmf workflow; bounded native recipe; use provenance("fit_nonnegative_ebmf")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examplesIf FALSE
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

#' Fit nonnegative activities with signed feature effects
#'
#' A thin wrapper around native FLASH priors. Activity and effect priors are
#' declared in semantic sample/feature order, then mapped to native axes once.
#' Observations may be negative. No centering or sign flipping is performed.
#' @inheritParams fit_ebmf
#' @param activity_prior Recognized nonnegative native EBNM prior.
#' @param effect_prior Signed native EBNM prior, or a custom function whose
#'   support remains unknown and whose stability refits require a callback.
#' @param initialization Only native_greedy is supported.
#' @param ... Other fit_ebmf arguments, excluding ebnm_fn.
#' @return A native flash fit with nonnegative activity and declared effect support.
#' @author Mikhael Manurung. Native FLASH: Jason Willwerscheid, Peter Carbonetto,
#'   Wei Wang and Matthew Stephens; native EBNM authors are listed in credits.
#' @details Origin: WRAPPER / call_only. No upstream implementation source copied.
#'   Use provenance("fit_seminonnegative_ebmf") for source pins and authors.
#' @keywords internal
#' @examplesIf FALSE
#' set.seed(3)
#' A <- matrix(rexp(24), 12, 2)
#' B <- matrix(rnorm(16), 8, 2)
#' X <- tcrossprod(A, B) + matrix(rnorm(96, sd = .2), 12, 8)
#' dimnames(X) <- list(paste0("s", 1:12), paste0("g", 1:8))
#' fit <- fit_seminonnegative_ebmf(X, "rows", max_factors = 2,
#'   backfit = TRUE, seed = 4, verbose = 0, feature_scale = "simulated_intensity")
#' factor_activity(fit)
fit_seminonnegative_ebmf <- function(X, sample_side,
                                   activity_prior = ebnm::ebnm_point_exponential,
                                   effect_prior = ebnm::ebnm_point_laplace,
                                   initialization = "native_greedy", ...) {
  sample_side <- match.arg(sample_side, c("rows", "columns"))
  if (!identical(initialization, "native_greedy")) stop("Only native_greedy initialization is supported")
  if (!is.function(activity_prior) || .prior_support(activity_prior) != "nonnegative") {
    stop("activity_prior must be a recognized nonnegative native EBNM function")
  }
  if (!is.function(effect_prior)) stop("effect_prior must be an EBNM function")
  effect_support <- .prior_support(effect_prior)
  if (effect_support == "nonnegative") stop("Use fit_nonnegative_ebmf for nonnegative effects")
  args <- list(...)
  if ("ebnm_fn" %in% names(args)) stop("Use activity_prior and effect_prior instead of ebnm_fn")
  priors <- list(activity_prior, effect_prior)
  if (sample_side == "columns") priors <- rev(priors)
  fit <- do.call(fit_ebmf, c(list(X = X, sample_side = sample_side, ebnm_fn = priors), args))
  if (effect_support == "unknown") attr(fit, "flashbridge_metadata")$refit_policy <- "callback_required"
  .attach_provenance(fit, "fit_seminonnegative_ebmf", match.call(),
    .resolved_parameters("fit_seminonnegative_ebmf", environment()))
}

#' Backfit named existing factors
#' @param fit Native fit with retained native data.
#' @param factors Optional factor IDs; defaults to all.
#'   Default refitting can update an intercept-only fit. Explicit adaptive
#'   selections require at least one biological factor.
#' @param control Named public flash_backfit settings.
#' @param seed Optional scoped seed.
#' @return New native fit with retained factor identity and new numerical identity.
#' @keywords internal
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: WRAPPER / call_only.
#'   Credits flash_backfit; use provenance("refit_factors")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examplesIf FALSE
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
  if (!is.null(factors) || is.null(attr(fit, "flashbridge_metadata")$intercept_index)) .require_biological_factors(ids)
  kset <- .native_factor_indices(fit)[.select_ids(factors, ids, "factor")]
  if (is.null(factors)) kset <- c(attr(fit, "flashbridge_metadata")$intercept_index, kset)
  if (!length(kset)) stop("No factors to refit")
  if (!is.list(control) || (length(control) && (is.null(names(control)) || anyDuplicated(names(control)))) || any(!names(control) %in% setdiff(names(formals(flashier::flash_backfit)), c("flash", "kset")))) stop("Unknown or protected backfit control")
  result <- .with_seed(seed, do.call(flashier::flash_backfit, c(list(flash = fit, kset = kset), control)))
  .attach_provenance(.native_changed(fit, result, ids, "backfit"), "refit_factors", match.call(), .resolved_parameters("refit_factors", environment()))
}

#' Fix or unfix a named native factor side
#' @param fit Native fit with recorded sample side.
#' @param factors Required factor IDs.
#' @param side activity or effects.
#' @param ids Optional sample/feature IDs for partial fixing; NULL fixes every entry of the side.
#' @param use_fixed_in_ebnm Native handling of fixed values in prior estimation.
#' @param action fix or unfix; native unfix operates on the selected pair.
#' @param sample_side Explicit rows or columns for a bare native fit.
#' @return New native fit retaining named factor identity.
#' @keywords internal
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: WRAPPER / call_only.
#'   Credits flash_factors_fix; flash_factors_unfix; use provenance("fix_factors")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examplesIf FALSE
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
fix_factors <- function(fit, factors, side, ids = NULL, use_fixed_in_ebnm = NULL, action = "fix", sample_side = NULL) {
  if (is.null(factors)) stop("factors is required; supply the factor IDs to ", action)
  action <- match.arg(action, c("fix", "unfix"))
  side <- match.arg(side, c("activity", "effects"))
  all_ids <- .factor_ids(fit)
  kset <- .native_factor_indices(fit)[.select_ids(factors, all_ids, "factor")]
  if (!length(kset)) stop("No factors selected")
  sample_side <- .resolve_sample_side(fit, sample_side)
  mode <- if ((side == "activity") == (sample_side == "rows")) "loadings" else "factors"
  if (action == "unfix") {
    if (!is.null(ids) || !is.null(use_fixed_in_ebnm)) stop("Native unfix applies to the entire selected factor pair; partial IDs unsupported")
    result <- flashier::flash_factors_unfix(fit, kset = kset)
  } else {
    view <- standardize_factors(fit, sample_side = sample_side, scaling = "raw")
    axis_ids <- rownames(view[[side]])
    if (!is.null(ids) && !length(ids)) stop("ids must name at least one ", side, " ID; use NULL to fix every entry")
    # native fixed_idx is a full-axis logical mask, not integer positions
    fixed_idx <- if (is.null(ids)) NULL else seq_along(axis_ids) %in% .select_ids(ids, axis_ids, side)
    result <- flashier::flash_factors_fix(fit, kset = kset, which_dim = mode, fixed_idx = fixed_idx, use_fixed_in_ebnm = use_fixed_in_ebnm)
  }
  result <- .native_changed(fit, result, all_ids, action)
  attr(result, "flashbridge_metadata")$sample_side <- sample_side
  .attach_provenance(result, "fix_factors", match.call(), .resolved_parameters("fix_factors", environment()))
}

#' Remove named factors without refitting
#' @param fit Native flash fit.
#' @param factors Unique factor IDs to remove.
#' @return Native fit retaining remaining factor IDs.
#' @keywords internal
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: WRAPPER / call_only.
#'   Credits flash_factors_remove; use provenance("remove_factors")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examplesIf FALSE
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
  if (is.null(factors)) stop("factors is required; supply the factor IDs to remove")
  ids <- .factor_ids(fit)
  selected <- .select_ids(factors, ids, "factor")
  kset <- .native_factor_indices(fit)[selected]
  if (!length(kset)) stop("No factors selected")
  .attach_provenance(.native_changed(fit, flashier::flash_factors_remove(fit, kset), ids[-selected], "remove"), "remove_factors", match.call(), .resolved_parameters("remove_factors", environment()))
}

#' Reorder native factors by exact named permutation
#' @param fit Native flash fit.
#' @param order Exact permutation of retained factor IDs.
#' @return Native fit with stable identities in requested order.
#' @keywords internal
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: WRAPPER / call_only.
#'   Credits flash_factors_reorder; use provenance("reorder_factors")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examplesIf FALSE
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
  native <- .native_factor_indices(fit)
  kset <- c(setdiff(seq_len(fit$n_factors), native), native[kset])
  .attach_provenance(.native_changed(fit, flashier::flash_factors_reorder(fit, kset), order, "reorder"), "reorder_factors", match.call(), .resolved_parameters("reorder_factors", environment()))
}
