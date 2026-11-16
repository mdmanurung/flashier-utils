#' flashier.utils: Transparent Interpretation of Empirical Bayes Matrix Factorizations
#'
#' Fit or import a 'flashier' matrix factorization, standardize it into a named
#' view (samples by factors activity, features by factors effects), then interpret,
#' compare, validate and transfer the factors. Start with
#' [fit_ebmf()] and [standardize_factors()]; the site's
#' Get started tutorial walks through the workflow.
#'
#' Every exported function's reference page names its original authors;
#' [factor_provenance()] returns the audited source registry. Full credits:
#' \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#'
#' @section Conventions:
#' Character arguments that take one of a fixed set of values, such as
#' `sample_side` or `scaling`, are case-sensitive and checked with
#' [match.arg()]. Group columns in metadata may be factors or character;
#' factor level order is ignored and only the grouping is used. Only dimnames
#' of an input matrix are kept; other attributes are dropped. Outputs hold no
#' `NA` or `Inf` unless a function's documentation says so.
#'
#' @section Glossary:
#' \describe{
#'   \item{view}{A named `factor_representation`: samples by factors activity and
#'     features by factors effects that share factor columns and reconstruct
#'     the fit.}
#'   \item{activity, effects}{The sample-side and feature-side loadings of each
#'     factor.}
#'   \item{sample_side}{Whether samples are the `"rows"` or the `"columns"` of the
#'     fitted matrix.}
#'   \item{scaling}{The coordinate choice for factors: `"ldf"` (unit-norm
#'     factors), `"raw"` (flashier's own) or `"rms"` (unit root-mean-square
#'     activity). A scaling changes units, not the fitted model.}
#'   \item{PVE}{Proportion of variance explained as flashier computes it. It is
#'     not a unique variance decomposition and not held-out prediction.}
#'   \item{lfsr}{Local false sign rate, the posterior probability that the sign
#'     of an effect is wrong (Stephens 2017). Reported only for the EB-SNMF
#'     point estimates that carry it.}
#'   \item{specificity}{A descriptive one-versus-rest measure of how
#'     concentrated a factor's activity is in a group; see [factor_specificity()].}
#'   \item{congruence}{Tucker congruence, the cosine similarity of two factors'
#'     effect vectors on shared named features; see [factor_replicability()].}
#' }
#'
#' @keywords internal
#' @author Mikhael Manurung; native engine Jason Willwerscheid, Peter Carbonetto,
#'   Wei Wang and Matthew Stephens; see CREDITS.md for every upstream author.
"_PACKAGE"

#' srr standards evidence for flashier.utils-package.R
#'
#' @srrstats {G1.2} README.Rmd has a Lifecycle section: experimental status, which APIs may change, and that NEWS.md lists every breaking change and renamed functions stay as deprecated aliases for one release.
#' @srrstats {G1.3} The package help page (R/flashier.utils-package.R) has a Glossary section defining view, activity, effects, sample_side, scaling, PVE, lfsr, specificity and congruence.
#' @srrstats {G1.4} roxygen2 documents the exported functions, with markdown enabled and the srr roclet set in the Roxygen field of DESCRIPTION. See the full block above R/fit.R fit_ebmf().
#' @srrstats {G2.3} All fixed-choice character arguments use match.arg() (G2.3a) and the package help page states they are case-sensitive (G2.3b).
#' @srrstats {G2.3b} The Conventions section of the package help page (R/flashier.utils-package.R) states that fixed-choice character arguments are case-sensitive and checked with match.arg().
#' @noRd
NULL
