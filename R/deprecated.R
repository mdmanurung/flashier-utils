# Deprecated aliases kept for one release after the 2026-11 renames.
#' Deprecated function names
#'
#' `provenance()` and `test_factors()` were renamed to [factor_provenance()]
#' and [test_factor_association()], because other CRAN packages export the
#' same names. The old names still work for one release and warn.
#' @param ... Arguments passed to the renamed function.
#' @return The renamed function's result.
#' @name flashier.utils-deprecated
#' @author Mikhael Manurung.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @keywords internal
NULL

#' @rdname flashier.utils-deprecated
#' @export
provenance <- function(...) {
  .Deprecated("factor_provenance", package = "flashier.utils")
  factor_provenance(...)
}

#' @rdname flashier.utils-deprecated
#' @export
test_factors <- function(...) {
  .Deprecated("test_factor_association", package = "flashier.utils")
  test_factor_association(...)
}
