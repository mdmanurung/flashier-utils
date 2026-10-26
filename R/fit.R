# Public flashier calls only; native prior list order is never reinterpreted.
.prior_support <- function(fn) {


  signed <- c("ebnm_normal", "ebnm_point_normal", "ebnm_point_laplace", "ebnm_unimodal",
    "ebnm_unimodal_symmetric", "ebnm_normal_scale_mixture", "ebnm_horseshoe",
    "ebnm_npmle", "ebnm_flat", "ebnm_ash")
  nonnegative <- c("ebnm_point_exponential", "ebnm_unimodal_nonnegative", "ebnm_generalized_binary")
  for (name in c(signed, nonnegative)) {
    if (name %in% getNamespaceExports("ebnm") && identical(fn, getExportedValue("ebnm", name)))
      return(if (name %in% signed) "signed" else "nonnegative")
  }
  "unknown"
}

# FU-0: fail fast when the installed flashier lacks a capability the fit
# pipeline will need LATE. flash_fit_get_p2m/flash_factors_remove are only
# called inside .nullcheck_adaptive() after an expensive fit + nullcheck, so
# a missing symbol there would waste the whole run. Checked once at fit
# entry; trivial extractors are deliberately not preflighted.
.ns_has_function <- function(pkg, sym) {
  exists(sym, envir = asNamespace(pkg), mode = "function")
}
.flashier_preflight <- function(need_nullcheck_moments = FALSE) {
  core <- c("flash_init", "flash_set_verbose", "flash_add_intercept",
            "flash_greedy", "flash_backfit", "flash_nullcheck")
  late <- c("flash_fit_get_p2m", "flash_factors_remove")
  need <- if (isTRUE(need_nullcheck_moments)) c(core, late) else core
  missing <- need[!vapply(need, .ns_has_function, logical(1), pkg = "flashier")]
  if (length(missing)) {
    stop(sprintf(paste0(
      "The installed flashier (%s) lacks native capabilities required for this ",
      "fit: %s. Install the pinned upstream flashier ",
      "(willwerscheid/flashier@fd44811ab22d61ee905e1964175571c452fe7683) or ",
      "disable the requesting feature (e.g. nullcheck = FALSE)."),
      as.character(utils::packageVersion("flashier")), paste(missing, collapse = ", ")),
      call. = FALSE)
  }
  invisible(TRUE)
}

.nullcheck_adaptive <- function(fit, kset, control = list(), verbose = 0) {
  retained <- seq_len(fit$n_factors)
  if (!is.null(control$kset)) {
    if (any(!control$kset %in% kset)) stop("Nullcheck cannot select protected factors")
    kset <- control$kset
  }
  remove <- if (is.null(control$remove)) TRUE else control$remove
  if (!is.logical(remove) || length(remove) != 1L || is.na(remove)) stop("remove must be TRUE or FALSE")
  control$kset <- control$remove <- NULL
  if (!length(kset)) return(list(fit = fit, retained = retained))
  fit <- do.call(flashier::flash_nullcheck,
    c(list(flash = fit, kset = kset, remove = FALSE, verbose = verbose), control))
  if (remove) {
    moments <- lapply(1:2, function(n) flashier::flash_fit_get_p2m(flashier::flash_fit(fit), n))
    for (k in sort(unique(kset), decreasing = TRUE)) {
      if (all(moments[[1]][, k] == 0) && all(moments[[2]][, k] == 0)) {
        fit <- flashier::flash_factors_remove(fit, kset = k)
        retained <- retained[-k]
      }
    }
  }
  list(fit = fit, retained = retained)
}

#' Fit a named matrix using native empirical Bayes factorization
#'
#' Uses public flashier fitting calls. max_factors is an upper bound. Priors in a
#' list are in native row/column order, independent of sample_side. A supplied
#' seed restores the caller's RNG state; NULL consumes it. No preprocessing.
#' @param X Named numeric matrix or double sparse Matrix; NA is observed missingness.
#' @param sample_side Required rows or columns identifying samples.
#' @param ebnm_fn Native EBNM function or native row/column list.
#' @param S Native observation standard errors.
#' @param var_type Native residual variance setting.
#' @param max_factors Nonnegative upper bound on native greedy rank.
#' @param backfit Run native backfit.
#' @param nullcheck Run native null check.
#' @param id_policy require IDs, or explicitly generate them.
#' @param feature_scale Meaning and units of the supplied matrix.
#' @param independent_unit Biological unit description; recorded, not inferred.
#' @param seed Optional scoped seed.
#' @param verbose Native verbosity.
#' @param control Named pipeline settings: S_dim, greedy, backfit, nullcheck lists.
#' @param intercept none (default) or feature. A feature intercept has fixed
#'   sample-side ones and a signed native point-normal feature offset. It is
#'   excluded from semantic factor IDs; max_factors counts adaptive factors.
#' @return A native flash object with flashbridge_metadata attribute. Sparse input
#'   with explicit S is unsupported by the audited native engine and refused.
#' @importClassesFrom Matrix sparseMatrix dMatrix
#' @export
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: WRAPPER / call_only.
#'   Credits flash; use provenance("fit_ebmf")
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
#' fit
fit_ebmf <- function(X, sample_side = c("rows", "columns"), ebnm_fn = ebnm::ebnm_point_normal,
                     S = NULL, var_type = 0L, max_factors = 50L, backfit = FALSE, nullcheck = TRUE,
                     id_policy = "require", feature_scale = "unknown", independent_unit = NULL,
                     seed = NULL, verbose = 1L, intercept = c("none", "feature"), control = list()) {
  if (missing(sample_side)) stop("sample_side must be explicitly chosen")
  sample_side <- match.arg(sample_side)
  intercept <- match.arg(intercept)
  generated <- c(rows = is.null(rownames(X)), columns = is.null(colnames(X)))
  X <- .validate_matrix(X, id_policy)
  .check_S(S)
  if (!is.null(var_type) && (!is.numeric(var_type) || !length(var_type) || !all(var_type %in% 0:2) || anyDuplicated(var_type))) stop("var_type must be NULL, 0, 1, 2 or c(1, 2)")
  if (inherits(X, "sparseMatrix") && !is.null(S)) stop("Sparse observations with explicit S are unsupported by native flashier (S4 subscript); supply explicitly dense input within your memory budget")
  if (length(max_factors) != 1L || !is.numeric(max_factors) || !is.finite(max_factors) || max_factors < 0 || max_factors != as.integer(max_factors)) stop("max_factors must be a nonnegative integer")
  if (!is.list(control) || (length(control) && (is.null(names(control)) || anyDuplicated(names(control))))) stop("control must be a uniquely named list")
  if (any(!names(control) %in% c("S_dim", "greedy", "backfit", "nullcheck"))) stop("Unknown control setting")
  validate_control <- function(values, fn, protected) {
    if (is.null(values)) return(list())
    if (!is.list(values) || (length(values) && (is.null(names(values)) || anyDuplicated(names(values))))) stop("Pipeline control must be a uniquely named list")
    if (any(!names(values) %in% setdiff(names(formals(fn)), protected))) stop("Unknown or protected native pipeline argument")
    values
  }
  greedy <- validate_control(control$greedy, flashier::flash_greedy, c("flash", "Kmax", "ebnm_fn", "verbose"))
  refit <- validate_control(control$backfit, flashier::flash_backfit, c("flash", "verbose"))
  check <- validate_control(control$nullcheck, flashier::flash_nullcheck, c("flash", "verbose"))
  if (intercept == "feature" && !is.null(check$kset)) stop("Intercept fitting selects adaptive nullcheck factors; kset is protected")
  .flashier_preflight(need_nullcheck_moments = isTRUE(nullcheck))
  priors <- if (is.list(ebnm_fn)) ebnm_fn else rep(list(ebnm_fn), 2L)
  if (length(priors) != 2L) stop("Native prior list must contain row and column functions")
  intercept_ids <- NULL
  fit <- .with_seed(seed, {
    if (intercept == "feature") {
      object <- flashier::flash_init(X, S = S, var_type = var_type, S_dim = control$S_dim)
      object <- flashier::flash_set_verbose(object, verbose = verbose)
      object <- flashier::flash_add_intercept(object, rowwise = sample_side == "rows")
      if (max_factors > 0L) {
        object <- do.call(flashier::flash_greedy, c(list(flash = object, Kmax = max_factors,
          ebnm_fn = ebnm_fn, verbose = verbose), greedy))
      }
      intercept_ids <- if (object$n_factors > 1L) paste0("F", seq_len(object$n_factors - 1L)) else character()
      if (backfit) object <- do.call(flashier::flash_backfit, c(list(flash = object, verbose = verbose), refit))
      if (nullcheck) {
        checked <- .nullcheck_adaptive(object, seq_len(object$n_factors)[-1L], check, verbose)
        object <- checked$fit
        intercept_ids <- intercept_ids[checked$retained[checked$retained != 1L] - 1L]
      }
      object
    } else if (max_factors == 0L) {
      flashier::flash_init(X, S = S, var_type = var_type, S_dim = control$S_dim)
    } else if (!length(control)) {
      flashier::flash(X, S = S, ebnm_fn = ebnm_fn, var_type = var_type,
                     greedy_Kmax = max_factors, backfit = backfit, nullcheck = nullcheck, verbose = verbose)
    } else {
      object <- flashier::flash_init(X, S = S, var_type = var_type, S_dim = control$S_dim)
      object <- do.call(flashier::flash_greedy, c(list(flash = object, Kmax = max_factors, ebnm_fn = ebnm_fn, verbose = verbose), greedy))
      if (backfit) object <- do.call(flashier::flash_backfit, c(list(flash = object, verbose = verbose), refit))
      if (nullcheck) object <- do.call(flashier::flash_nullcheck, c(list(flash = object, verbose = verbose), check))
      object
    }
  })
  support <- vapply(priors, .prior_support, character(1))
  if (sample_side == "columns") support <- rev(support)
  names(support) <- c("activity", "effects")
  raw <- .raw_pair(fit)
  attr(fit, "flashbridge_metadata") <- c(list(
    fit_id = digest::digest(list(raw, sample_side, Sys.time()), algo = "sha256"),
    factor_ids = if (fit$n_factors == 0L) character() else paste0("F", seq_len(fit$n_factors)),
    sample_side = sample_side, row_ids = rownames(X), column_ids = colnames(X),
    generated_ids = generated & id_policy == "generate", feature_scale = feature_scale,
    independent_unit = independent_unit, preprocessing = list(), prior_support = support,
    # A caller greedy initializer is code the recorded configuration cannot audit.
    recipe_version = "1.0.0", refit_policy = if (is.null(control$greedy$init_fn)) "ordinary" else "callback_required",
    prior_identifiers = vapply(priors, function(fn) {
      exports <- getNamespaceExports("ebnm")
      matches <- exports[vapply(exports, function(n) identical(fn, getExportedValue("ebnm", n)), logical(1))]
      if (length(matches)) paste(matches, collapse = ";") else "custom_callback_required"
    }, character(1)),
    seed = seed, rng_kind = RNGkind(),
    configuration = list(S = S, var_type = var_type, max_factors = max_factors,
                         backfit = backfit, nullcheck = nullcheck, id_policy = id_policy,
                         verbose = verbose, control = control)), .engine_metadata())
  if (intercept == "feature") {
    metadata <- attr(fit, "flashbridge_metadata")
    metadata$factor_ids <- intercept_ids
    metadata$native_indices <- seq_len(fit$n_factors)[-1L]
    metadata$intercept_index <- 1L
    metadata$native_factor_ids <- c("Intercept", intercept_ids)
    metadata$factor_roles <- stats::setNames(c("offset", rep("factor", length(intercept_ids))), metadata$native_factor_ids)
    metadata$offset_model <- "feature"
    metadata$configuration$intercept <- "feature"
    metadata$support_version <- "1.0.0"
    metadata$factor_support <- .support_records(metadata$factor_ids,
      metadata$native_indices, metadata$prior_support, metadata$factor_roles)
    attr(fit, "flashbridge_metadata") <- metadata
  }
  .attach_provenance(fit, "fit_ebmf", match.call(), .resolved_parameters("fit_ebmf", environment()))
}

.multistart_initializer <- function(method, seed, constraints) {
  # Lazy promises would retain the caller's data and all candidates in saved fits.
  force(method)
  force(seed)
  force(constraints)
  trace <- new.env(parent = emptyenv())
  trace$hashes <- character()
  callback <- function(flash) {
    values <- switch(method,
      native_random = flashier::flash_greedy_init_default(flash,
        sign_constraints = constraints, seed = seed),
      softImpute = flashier::flash_greedy_init_softImpute(flash, seed = seed),
      irlba = flashier::flash_greedy_init_irlba(flash, seed = seed))
    trace$hashes <- c(trace$hashes, digest::digest(values, algo = "sha256"))
    values
  }
  list(callback = callback, trace = trace)
}

#' Compare seeded native empirical Bayes optimization starts
#'
#' Each candidate uses a distinct native initializer seed. Native-random starts
#' preserve audited prior sign constraints. softImpute and IRLBA starts require
#' signed priors on both native sides; IRLBA also requires dense complete input.
#' Candidate failures are retained. ELBO ranks optimization solutions under the
#' fitted model and does not measure biological interpretability or uncertainty.
#' @inheritParams fit_ebmf
#' @param n_starts Integer number of candidates, at least two.
#' @param init_method native_random, softImpute or irlba.
#' @param init_seeds Optional unique nonnegative integer initializer seeds, one
#'   per candidate. Supply these or a master seed.
#' @param seed Optional scoped master seed used to generate initializer seeds.
#' @param select elbo chooses the largest finite native ELBO; exact ties choose
#'   the smallest start ID. none leaves best_fit and selected_start unset.
#' @param keep_fits Retain all successful native candidates; the selected best_fit
#'   is available even when FALSE.
#' @param max_factors Positive bound on each candidate's greedy rank.
#' @param control Native pipeline settings; greedy init_fn is owned by this wrapper.
#' @return A multistart_result with best_fit, selected_start, candidates, optional
#'   fits and analysis_metadata. All-fail errors retain their candidates table.
#'   Candidate digests identify numerical coordinates and exclude elapsed time
#'   and time-based fit IDs. Initializer traces include native validation probes.
#'   Stability refits require an explicit callback for the full selection recipe.
#' @author Mikhael Manurung; native FLASH authors Jason Willwerscheid,
#'   Peter Carbonetto, Wei Wang and Matthew Stephens. No native source is copied.
#'   Optional initializers: Trevor Hastie, Rahul Mazumder and Balasubramanian
#'   Narasimhan (softImpute); Jim Baglama, Lothar Reichel and B. W. Lewis (irlba).
#' @keywords internal
#' @examplesIf FALSE
#' set.seed(3)
#' X <- matrix(rnorm(192), 16, 12,
#'   dimnames = list(paste0("s", 1:16), paste0("g", 1:12)))
#' starts <- fit_multistart_ebmf(X, "rows", n_starts = 2, max_factors = 1,
#'   seed = 4, verbose = 0)
#' starts$candidates[, c("start", "init_seed", "status", "rank", "elbo", "selected")]
fit_multistart_ebmf <- function(X, sample_side, ebnm_fn = ebnm::ebnm_point_normal,
                              S = NULL, var_type = 0L, max_factors = 50L,
                              backfit = FALSE, nullcheck = TRUE, id_policy = "require",
                              feature_scale = "unknown", independent_unit = NULL,
                              n_starts = 10L, init_method = c("native_random", "softImpute", "irlba"),
                              init_seeds = NULL, seed = NULL, select = c("elbo", "none"),
                              keep_fits = TRUE, verbose = 1L, control = list()) {
  if (missing(sample_side)) stop("sample_side must be explicitly chosen")
  sample_side <- match.arg(sample_side, c("rows", "columns"))
  init_method <- match.arg(init_method)
  select <- match.arg(select)
  .program_number(n_starts, "n_starts", TRUE, 2)
  .program_number(max_factors, "max_factors", TRUE, 1)
  if (!is.logical(keep_fits) || length(keep_fits) != 1L || is.na(keep_fits)) stop("keep_fits must be TRUE or FALSE")
  generated <- c(rows = is.null(rownames(X)), columns = is.null(colnames(X)))
  X <- .validate_matrix(X, id_policy)
  priors <- if (is.list(ebnm_fn)) ebnm_fn else rep(list(ebnm_fn), 2L)
  if (length(priors) != 2L || !all(vapply(priors, is.function, logical(1)))) stop("Native prior list must contain two functions")
  support <- vapply(priors, .prior_support, character(1))
  if (any(support == "unknown")) stop("Multistart requires recognized native prior support")
  constraints <- ifelse(support == "nonnegative", 1, 0)
  if (init_method != "native_random" && any(support != "signed")) stop("softImpute and IRLBA initializers require signed priors on both sides")
  if (init_method == "irlba" && (inherits(X, "sparseMatrix") || anyNA(X))) stop("IRLBA initialization requires dense complete input")
  if (init_method != "native_random") .check_engine(init_method, init_method)
  if (!is.list(control) || (length(control) && (is.null(names(control)) || anyDuplicated(names(control))))) stop("control must be a uniquely named list")
  if (!is.null(control$greedy) && !is.list(control$greedy)) stop("greedy control must be a named list")
  if (!is.null(control$greedy$init_fn)) stop("Multistart owns greedy init_fn; caller initializer conflicts")
  if (is.null(init_seeds)) {
    if (is.null(seed)) stop("Supply seed or unique init_seeds")
    init_seeds <- .with_seed(seed, sample.int(.Machine$integer.max, n_starts))
  }
  if (!is.numeric(init_seeds) || length(init_seeds) != n_starts || anyNA(init_seeds) ||
      any(!is.finite(init_seeds)) || any(init_seeds < 0 | init_seeds > .Machine$integer.max | init_seeds != floor(init_seeds)) || anyDuplicated(init_seeds)) {
    stop("init_seeds must contain one unique nonnegative integer per start")
  }
  init_seeds <- as.integer(init_seeds)
  fits <- vector("list", n_starts)
  candidates <- data.frame(start = seq_len(n_starts), init_seed = init_seeds,
    status = "failed", reason = NA_character_, rank = NA_integer_, elbo = NA_real_,
    fit_id = NA_character_, candidate_digest = NA_character_,
    initialization_digest = NA_character_, initializer_calls = 0L,
    elapsed_seconds = NA_real_, selected = FALSE)
  for (i in seq_len(n_starts)) {
    initializer <- .multistart_initializer(init_method, init_seeds[i], constraints)
    settings <- control
    settings$greedy$init_fn <- initializer$callback
    started <- proc.time()[[3L]]
    tryCatch({
      fitted <- fit_ebmf(X, sample_side, ebnm_fn, S, var_type, max_factors,
        backfit, nullcheck, id_policy, feature_scale, independent_unit,
        seed = init_seeds[i], verbose = verbose, control = settings)
      candidates$rank[i] <- fitted$n_factors
      candidates$fit_id[i] <- attr(fitted, "flashbridge_metadata")$fit_id
      elbo <- flashier::flash_fit_get_elbo(flashier::flash_fit(fitted))
      if (!is.numeric(elbo) || length(elbo) != 1L || !is.finite(elbo)) stop("Candidate native ELBO is not finite")
      candidates$elbo[i] <- elbo
      candidates$candidate_digest[i] <- digest::digest(list(.raw_pair(fitted),
        rank = fitted$n_factors, elbo = elbo), algo = "sha256")
      metadata <- attr(fitted, "flashbridge_metadata")
      metadata$refit_policy <- "callback_required"
      metadata$generated_ids <- generated & id_policy == "generate"
      metadata$multistart <- list(start = i, init_seed = init_seeds[i],
        init_method = init_method, n_starts = n_starts, selection = select)
      attr(fitted, "flashbridge_metadata") <- metadata
      fits[[i]] <- .attach_provenance(fitted, "fit_multistart_ebmf", match.call(),
        .resolved_parameters("fit_multistart_ebmf", environment()))
      candidates$status[i] <- "supported"
    }, error = function(e) candidates$reason[i] <<- conditionMessage(e))
    candidates$initializer_calls[i] <- length(initializer$trace$hashes)
    candidates$initialization_digest[i] <- digest::digest(initializer$trace$hashes, algo = "sha256")
    candidates$elapsed_seconds[i] <- proc.time()[[3L]] - started
  }
  successful <- which(candidates$status == "supported")
  if (!length(successful)) {
    error <- simpleError(paste("All multistart candidates failed:",
      paste(paste0(candidates$start, ": ", candidates$reason), collapse = "; ")))
    error$candidates <- candidates
    stop(error)
  }
  selected <- if (select == "elbo") successful[which.max(candidates$elbo[successful])] else NA_integer_
  best <- if (is.na(selected)) NULL else fits[[selected]]
  if (!is.na(selected)) candidates$selected[selected] <- TRUE
  result <- structure(list(best_fit = best, selected_start = selected, candidates = candidates,
    fits = if (keep_fits) fits else NULL,
    analysis_metadata = list(sample_side = sample_side, init_method = init_method,
      initializer_engine = if (init_method == "native_random") "flashier" else init_method,
      initializer_version = as.character(utils::packageVersion(if (init_method == "native_random") "flashier" else init_method)),
      init_seeds = init_seeds, native_sign_constraints = constraints,
      selection = select, refit_policy = "callback_required", engine = "flashier",
      feature_scale = feature_scale, independent_unit = independent_unit,
      interpretation = "optimization_sensitivity_not_inferential_uncertainty")), class = "multistart_result")
  .attach_provenance(result, "fit_multistart_ebmf", match.call(),
    .resolved_parameters("fit_multistart_ebmf", environment()))
}
