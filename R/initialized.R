# Public native pipeline only; FLASH authors credited below and in provenance.

#' Fit named initial factors, optionally fixing one semantic side
#'
#' Initial factors are backfitted before adding adaptive factors. Their IDs are
#' retained even when native estimation zeros them. Only newly added factors
#' are automatically nullchecked and removed. No intercept is added.
#' @inheritParams fit_ebmf
#' @param activity_init Named finite samples by factors matrix.
#' @param effect_init Named finite features by factors matrix, with identical
#'   factor names and order. Rows of both matrices are matched to X by ID.
#' @param activity_prior Native EBNM function for activities.
#' @param effect_prior Native EBNM function for effects.
#' @param fixed_side none, activity or effects; fixes that side of all initial factors.
#' @param add_factors Nonnegative upper bound on additional greedy factors.
#' @param nullcheck_new Whether to nullcheck newly added factors.
#' @param control Named S_dim, greedy, backfit and nullcheck settings; kset is
#'   protected for backfit/nullcheck so initialized factors cannot be bypassed.
#' @return Native flash fit with stable initialized IDs, role metadata and input
#'   hashes. Stability refits require an explicit callback recreating initialization.
#' @author Mikhael Manurung. Native FLASH: Jason Willwerscheid, Peter Carbonetto,
#'   Wei Wang and Matthew Stephens. EBNM authors are listed in credits.
#' @details Origin: WRAPPER / call_only; no upstream implementation source copied.
#'   Priors use semantic sample/feature order and map to native axes once.
#'   Recognized nonnegative priors require nonnegative initial values; custom
#'   prior support remains unknown. Neither initialized side may have a zero column.
#' @examplesIf FALSE
#' A <- cbind(program = seq(-1, 1, length.out = 12))
#' B <- cbind(program = c(2, -1, 0, 1, 3, -2, 4, 0))
#' rownames(A) <- paste0("s", 1:12)
#' rownames(B) <- paste0("g", 1:8)
#' X <- tcrossprod(A, B)
#' fit <- fit_initialized_ebmf(X, "rows", A, B, fixed_side = "effects",
#'   seed = 4, verbose = 0, feature_scale = "simulated_intensity")
#' factor_effects(fit, representation = list(scaling = "raw"))
#' @keywords internal
fit_initialized_ebmf <- function(X, sample_side, activity_init, effect_init,
                                 activity_prior = ebnm::ebnm_point_normal,
                                 effect_prior = ebnm::ebnm_point_laplace,
                                 fixed_side = c("none", "activity", "effects"),
                                 add_factors = 0L, nullcheck_new = TRUE,
                                 S = NULL, var_type = 0L, id_policy = "require",
                                 feature_scale = "unknown", independent_unit = NULL,
                                 seed = NULL, verbose = 1L, control = list()) {
  sample_side <- match.arg(sample_side, c("rows", "columns"))
  fixed_side <- match.arg(fixed_side)
  generated <- c(rows = is.null(rownames(X)), columns = is.null(colnames(X)))
  X <- .validate_matrix(X, id_policy)
  .program_number(add_factors, "add_factors", integer = TRUE)
  if (!is.logical(nullcheck_new) || length(nullcheck_new) != 1L || is.na(nullcheck_new)) stop("nullcheck_new must be TRUE or FALSE")
  if (!is.function(activity_prior) || !is.function(effect_prior)) stop("Initial priors must be native EBNM functions")
  if (!is.matrix(activity_init) || !is.matrix(effect_init)) stop("Initial factors must be named numeric matrices")
  activity_init <- .validate_matrix(activity_init, allow_missing = FALSE)
  effect_init <- .validate_matrix(effect_init, allow_missing = FALSE)
  ids <- colnames(activity_init)
  if (!identical(ids, colnames(effect_init))) stop("Initial factor names and order must match")
  if ("Intercept" %in% ids) stop("Intercept is a reserved native factor ID")
  samples <- if (sample_side == "rows") rownames(X) else colnames(X)
  features <- if (sample_side == "rows") colnames(X) else rownames(X)
  if (!setequal(rownames(activity_init), samples) || !setequal(rownames(effect_init), features)) stop("Initial sample/feature IDs must exactly match X")
  activity_init <- activity_init[samples, , drop = FALSE]
  effect_init <- effect_init[features, , drop = FALSE]
  pair <- list(activity_init, effect_init)
  priors <- list(activity_prior, effect_prior)
  support <- stats::setNames(vapply(priors, .prior_support, character(1)), c("activity", "effects"))
  for (side in 1:2) {
    if (any(colSums(abs(pair[[side]])) == 0)) stop("Initial ", names(support)[side], " columns must be nonzero")
    if (support[side] == "nonnegative" && any(pair[[side]] < 0)) stop("Initial ", names(support)[side], " violates nonnegative prior support")
  }
  if (!is.list(control)) stop("control must be a named pipeline list")
  for (key in c("greedy", "backfit", "nullcheck")) {
    if (!is.null(control[[key]]) && !is.list(control[[key]])) stop("Pipeline control must be a named list")
  }
  if (!is.null(control$backfit$kset) || !is.null(control$nullcheck$kset)) stop("Initialized backfit/nullcheck kset is protected")
  if (sample_side == "columns") {
    pair <- rev(pair)
    priors <- rev(priors)
  }
  metadata <- NULL
  roles <- stats::setNames(rep("initialized", length(ids)), ids)
  fit <- .with_seed(seed, {
    object <- fit_ebmf(X, sample_side, ebnm_fn = priors, S = S, var_type = var_type,
      max_factors = 0, nullcheck = FALSE, id_policy = id_policy, feature_scale = feature_scale,
      independent_unit = independent_unit, verbose = verbose, control = control)
    metadata <- attr(object, "flashbridge_metadata")
    object <- flashier::flash_set_verbose(object, verbose = verbose)
    object <- flashier::flash_factors_init(object, pair, ebnm_fn = priors)
    initial <- seq_along(ids)
    if (fixed_side != "none") {
      mode <- if ((fixed_side == "activity") == (sample_side == "rows")) "loadings" else "factors"
      object <- flashier::flash_factors_fix(object, initial, which_dim = mode)
    }
    object <- do.call(flashier::flash_backfit, c(list(flash = object, verbose = verbose), control$backfit))
    if (add_factors > 0) {
      object <- do.call(flashier::flash_greedy, c(list(flash = object, Kmax = add_factors, ebnm_fn = priors, verbose = verbose), control$greedy))
      n_new <- object$n_factors - length(initial)
      next_id <- 1L
      for (i in seq_len(n_new)) {
        while (paste0("F", next_id) %in% ids) next_id <- next_id + 1L
        id <- paste0("F", next_id)
        ids <- c(ids, id)
        roles[id] <- "factor"
      }
      object <- do.call(flashier::flash_backfit, c(list(flash = object, verbose = verbose), control$backfit))
      if (nullcheck_new && n_new > 0) {
        checked <- .nullcheck_adaptive(object, seq_len(object$n_factors)[-initial],
          if (is.null(control$nullcheck)) list() else control$nullcheck, verbose)
        object <- checked$fit
        ids <- ids[checked$retained]
        roles <- roles[ids]
      }
    }
    object
  })
  metadata$fit_id <- digest::digest(list(.raw_pair(fit), sample_side, Sys.time()), algo = "sha256")
  metadata$factor_ids <- metadata$native_factor_ids <- ids
  metadata$native_indices <- seq_along(ids)
  metadata$factor_roles <- roles
  metadata$support_version <- "1.0.0"
  metadata$factor_support <- .support_records(ids, metadata$native_indices,
    metadata$prior_support, roles)
  metadata$generated_ids <- generated & id_policy == "generate"
  metadata$refit_policy <- "callback_required"
  metadata$seed <- seed
  metadata$initialization <- list(activity_sha256 = digest::digest(activity_init, algo = "sha256"),
    effect_sha256 = digest::digest(effect_init, algo = "sha256"), factor_ids = colnames(activity_init))
  metadata$configuration <- list(S = S, var_type = var_type, fixed_side = fixed_side,
    add_factors = add_factors, nullcheck_new = nullcheck_new, id_policy = id_policy,
    verbose = verbose, control = control)
  attr(fit, "flashbridge_metadata") <- metadata
  .attach_provenance(fit, "fit_initialized_ebmf", match.call(),
    .resolved_parameters("fit_initialized_ebmf", environment()))
}
