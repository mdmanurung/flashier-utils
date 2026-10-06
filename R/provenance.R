# Newly written registry accessor; no upstream source copied.
#' Inspect code origins and runtime dependencies
#'
#' The source_commit column is the audited source pin. runtime_source_commit is
#' read from an installed package's RemoteSha and is unknown when absent.
#' @param function_name Optional public function names.
#' @param result Optional result carrying analysis_metadata.
#' @param include_dependencies Include installed dependency versions and SHAs.
#' @return A data frame of source and reuse records with runtime metadata attributes.
#' @export
#' @author Mikhael Manurung.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: NEW_QOL / none.
#'   Credits new registry interface; use provenance("provenance")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examples
#' provenance(c("standardize_factors","backproject_contrast"))
provenance <- function(function_name = NULL, result = NULL, include_dependencies = TRUE) {
  requested <- function_name
  path <- system.file("provenance.csv", package = "flashier.utils")
  if (!nzchar(path)) stop("Installed provenance registry is missing")
  out <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(requested)) {
    if (any(!requested %in% out[["function"]])) stop("Unknown provenance function")
    out <- out[match(requested, out[["function"]]), , drop = FALSE]
  }
  if (include_dependencies) {
    packages <- c("flashier", "ebnm", "Matrix", "digest", "ggplot2")
    attr(out, "dependencies") <- do.call(rbind, lapply(packages, function(p) {
      d <- utils::packageDescription(p)
      data.frame(package = p, version = d$Version,
                 runtime_source_commit = if (is.null(d$RemoteSha)) "unknown" else d$RemoteSha)
    }))
  }
  attr(out, "analysis_metadata") <- if (is.list(result) && !is.null(result$analysis_metadata))
    result$analysis_metadata else attr(result, "analysis_metadata")
  result_metadata <- attr(out,"analysis_metadata")
  out <- .attach_provenance(out, "provenance", match.call(), .resolved_parameters("provenance", environment()))
  attr(out,"accessor_metadata") <- attr(out,"analysis_metadata")
  attr(out,"analysis_metadata") <- result_metadata
  out
}

.engine_metadata <- function() {
  packages <- c("flashier", "ebnm")
  list(engine_versions = stats::setNames(lapply(packages, function(p) {
    as.character(utils::packageVersion(p))
  }), packages), source_commit_if_known = stats::setNames(lapply(packages, function(p) {
    sha <- utils::packageDescription(p)$RemoteSha
    if (is.null(sha)) "unknown" else sha
  }), packages))
}
