# Call-only adapter for fastTopics grade-of-membership differential expression.
.fasttopics_pin <- "5cf9a448b17a7eb8f139cc916e22cb3ad1766b70"

#' Grade-of-membership differential expression with fixed memberships
#'
#' Delegates to [fastTopics::de_analysis()], which fits one Poisson model per
#' feature with the supplied memberships held fixed, then summarizes log2 fold
#' changes by MCMC. Counts and memberships are validated before the native call;
#' nothing is normalized here.
#' @param counts Named count matrix or double sparse Matrix: finite, nonnegative,
#'   integer-valued, with a positive library size for every sample.
#' @param memberships Named numeric matrix of topic or program memberships for
#'   the same samples: finite, nonnegative, at least two columns, and every row
#'   summing to one within 1e-8. EBMF activities are not memberships; normalize
#'   them yourself only if that is scientifically justified.
#' @param sample_side Sample axis of `counts`: rows or columns. Memberships always
#'   have samples in rows.
#' @param seed Required nonnegative integer; the MCMC uses R's random numbers.
#'   The caller's random number state is restored afterwards.
#' @param shrink_method "ash" (default) or "none", as in the native function.
#' @param lfc_stat "le" (least-extreme, default), "vsnull", or one membership
#'   column name, as in the native function.
#' @param pseudocount Positive native pseudocount.
#' @param fit_method Native Poisson fitting method.
#' @param control Named native control settings from
#'   [fastTopics::de_analysis_control_default()]; `nc` is fixed at 1.
#' @param keep_native Retain the native result.
#' @return A list with a long `statistics` table (feature, factor, log2_fc,
#'   lower, upper, z, neg_log10_p, lfsr, svalue), `analysis_metadata`, and the
#'   native result when requested.
#' @details Statistics are conditional on the fixed memberships, so membership
#'   uncertainty is ignored. `lower`/`upper` are highest-posterior-density
#'   intervals at `control$conf.level` (native default 0.68), not 95% equal-tail
#'   intervals. With ash shrinkage `neg_log10_p` is unavailable; without it
#'   `lfsr` and `svalue` are unavailable. Native MCMC warnings, such as zero
#'   acceptance rates, are kept in `analysis_metadata$warnings`. No calibration
#'   claim is made for small inputs.
#' @author Mikhael Manurung (adapter). fastTopics: Peter Carbonetto, Kevin Luo,
#'   Kushal Dey and Matthew Stephens, with contributors listed in its
#'   DESCRIPTION. No engine source is copied.
#' @references \url{https://github.com/stephenslab/fastTopics};
#'   \url{https://stephenslab.github.io/fastTopics/reference/de_analysis.html}
#' @export
#' @examples
#' \dontrun{
#' de <- test_gom_de(counts, memberships, "rows", seed = 1)
#' subset(de$statistics, lfsr < 0.01)
#' }
test_gom_de <- function(counts, memberships, sample_side, seed, shrink_method = c("ash", "none"),
                        lfc_stat = "le", pseudocount = 0.01, fit_method = "scd",
                        control = list(), keep_native = FALSE) {
  if (missing(seed)) stop("Choose seed explicitly; GoM DE uses Monte Carlo sampling")
  sample_side <- match.arg(sample_side, c("rows", "columns"))
  shrink_method <- match.arg(shrink_method)
  X <- .validate_matrix(counts, allow_missing = FALSE)
  if (sample_side == "columns") X <- Matrix::t(X)
  values <- if (inherits(X, "sparseMatrix")) methods::slot(X, "x") else X
  if (any(values < 0) || any(values != round(values))) stop("counts must be nonnegative integers")
  if (any(Matrix::rowSums(X) <= 0)) stop("Every sample needs a positive library size")
  if (!is.matrix(memberships) || !is.numeric(memberships)) stop("memberships must be a numeric matrix")
  if (ncol(memberships) < 2L) stop("memberships need at least two columns")
  .validate_ids(rownames(memberships), "Membership sample")
  .validate_ids(colnames(memberships), "Membership factor")
  if (any(!is.finite(memberships)) || any(memberships < 0)) stop("memberships must be finite and nonnegative")
  if (!setequal(rownames(memberships), rownames(X))) stop("memberships and counts must cover the same sample IDs")
  memberships <- memberships[rownames(X), , drop = FALSE]
  if (any(abs(rowSums(memberships) - 1) > 1e-8)) stop("Every membership row must sum to one; EBMF activities are not normalized here")
  if (length(pseudocount) != 1L || !is.finite(pseudocount) || pseudocount <= 0) stop("pseudocount must be positive")
  if (!is.logical(keep_native) || length(keep_native) != 1L || is.na(keep_native)) stop("keep_native must be TRUE or FALSE")
  .check_engine("fastTopics", c("de_analysis", "de_analysis_control_default"))
  if (utils::packageVersion("fastTopics") < "0.7.46") stop("test_gom_de needs fastTopics >= 0.7-46 (stephenslab/fastTopics@", .fasttopics_pin, ")")
  defaults <- fastTopics::de_analysis_control_default()
  if (!is.list(control) || (length(control) && (is.null(names(control)) || anyDuplicated(names(control)) ||
      any(!names(control) %in% names(defaults))))) stop("Unknown fastTopics de_analysis control setting")
  control$nc <- 1
  if (!is.matrix(X)) X <- methods::as(X, "CsparseMatrix")
  warnings <- character()
  native <- .with_seed(seed, tryCatch(withCallingHandlers(
    fastTopics::de_analysis(memberships, X, pseudocount = pseudocount, fit.method = fit_method,
      shrink.method = shrink_method, lfc.stat = lfc_stat, control = control, verbose = FALSE),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w))),
    error = function(e) stop("fastTopics de_analysis failed: ", conditionMessage(e), call. = FALSE)))
  column <- function(x) if (is.matrix(x)) as.vector(x) else rep(NA_real_, length(native$postmean))
  stats <- data.frame(feature = rep(rownames(native$postmean), ncol(native$postmean)),
    factor = rep(colnames(native$postmean), each = nrow(native$postmean)),
    log2_fc = as.vector(native$postmean), lower = column(native$lower), upper = column(native$upper),
    z = column(native$z), neg_log10_p = column(native$lpval), lfsr = column(native$lfsr),
    svalue = column(native$svalue))
  info <- utils::packageDescription("fastTopics")
  out <- list(statistics = stats, native = if (keep_native) native else NULL,
    analysis_metadata = list(engine = "fastTopics", inference = "conditional_on_fixed_memberships", lfc_scale = "log2",
      lfc_stat = lfc_stat, shrink_method = shrink_method, pseudocount = pseudocount, fit_method = fit_method,
      control = utils::modifyList(defaults, control), interval = "hpd", seed = seed,
      unavailable = if (shrink_method == "ash") "neg_log10_p (ash shrinkage)" else "lfsr and svalue (no shrinkage)",
      acceptance_rates = native$ar, warnings = warnings, sample_ids = rownames(X),
      engine_versions = list(fastTopics = info$Version),
      source_commit_if_known = list(fastTopics = if (is.null(info$RemoteSha)) "unknown" else info$RemoteSha)))
  .attach_provenance(out, "test_gom_de", match.call(), .resolved_parameters("test_gom_de", environment()))
}
