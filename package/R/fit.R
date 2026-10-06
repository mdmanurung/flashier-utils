# Public flashier calls only; native prior list order is never reinterpreted.
.prior_support <- function(fn) {
  signed <- c("ebnm_normal", "ebnm_point_normal", "ebnm_point_laplace", "ebnm_unimodal", "ebnm_ash")
  nonnegative <- c("ebnm_point_exponential", "ebnm_point_gamma", "ebnm_exponential")
  for (name in c(signed, nonnegative)) {
    if (name %in% getNamespaceExports("ebnm") && identical(fn, getExportedValue("ebnm", name)))
      return(if (name %in% signed) "signed" else "nonnegative")
  }
  "unknown"
}

#' Fit a named matrix using native empirical Bayes factorization
#'
#' Delegates to flashier::flash. max_factors is an upper bound. Priors in a
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
#' @return A native flash object with flashbridge_metadata attribute. Sparse input
#'   with explicit S is unsupported by the audited native engine and refused.
#' @importClassesFrom Matrix sparseMatrix dMatrix
#' @export
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
                     seed = NULL, verbose = 1L, control = list()) {
  if (missing(sample_side)) stop("sample_side must be explicitly chosen")
  sample_side <- match.arg(sample_side)
  generated <- c(rows = is.null(rownames(X)), columns = is.null(colnames(X)))
  X <- .validate_matrix(X, id_policy)
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
  fit <- .with_seed(seed, {
    if (max_factors == 0L) {
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
  priors <- if (is.list(ebnm_fn)) ebnm_fn else rep(list(ebnm_fn), 2L)
  if (length(priors) != 2L) stop("Native prior list must contain row and column functions")
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
    prior_identifiers = vapply(priors, function(fn) {
      exports <- getNamespaceExports("ebnm")
      matches <- exports[vapply(exports, function(n) identical(fn, getExportedValue("ebnm", n)), logical(1))]
      if (length(matches)) paste(matches, collapse = ";") else "custom_callback_required"
    }, character(1)),
    seed = seed, rng_kind = RNGkind(),
    configuration = list(S = S, var_type = var_type, max_factors = max_factors,
                         backfit = backfit, nullcheck = nullcheck, id_policy = id_policy,
                         verbose = verbose, control = control)), .engine_metadata())
  .attach_provenance(fit, "fit_ebmf", match.call(), .resolved_parameters("fit_ebmf", environment()))
}
