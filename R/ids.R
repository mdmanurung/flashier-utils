# Newly written integrity checks; matrices are never densified here.
.validate_ids <- function(ids, label) {
  if (is.null(ids) || !is.character(ids) || anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
    stop(label, " IDs must be unique, nonempty strings; supply explicit IDs")
  }
  ids
}

.validate_matrix <- function(X, id_policy = "require", allow_missing = TRUE) {
  id_policy <- match.arg(id_policy, c("require", "generate"))
  sparse <- inherits(X, "sparseMatrix") && inherits(X, "dMatrix")
  if (!(is.matrix(X) && is.numeric(X)) && !sparse) {
    stop("X must be a numeric matrix or a double sparse Matrix")
  }
  if (length(dim(X)) != 2L || any(dim(X) < 1L)) {
    stop("X requires two nonempty dimensions")
  }
  values <- if (sparse) methods::slot(X, "x") else X
  if (any(is.infinite(values)) || any(is.nan(values))) {
    stop("X contains Inf or NaN")
  }
  if (!allow_missing && anyNA(values)) {
    stop("Missing values are not supported for this operation")
  }
  for (side in 1:2) {
    ids <- dimnames(X)[[side]]
    if (is.null(ids) && id_policy == "generate") {
      ids <- paste0(if (side == 1L) "row" else "col", seq_len(dim(X)[side]))
      if (side == 1L) rownames(X) <- ids else colnames(X) <- ids
    }
    .validate_ids(ids, if (side == 1L) "Row" else "Column")
  }
  X
}

.match_ids <- function(
  ids,
  metadata,
  id_col = "sample",
  extra = "drop",
  na_action = "fail",
  required_columns = names(metadata)
) {
  .validate_ids(ids, "Required sample")
  if (!is.data.frame(metadata) || !id_col %in% names(metadata)) {
    stop("Metadata requires an explicit ", id_col, " column")
  }
  mids <- .validate_ids(as.character(metadata[[id_col]]), "Metadata sample")
  extra <- match.arg(extra, c("fail", "drop"))
  na_action <- match.arg(na_action, c("fail", "omit"))
  if (!all(required_columns %in% names(metadata))) {
    stop("Required metadata column is missing")
  }
  missing <- setdiff(ids, mids)
  if (length(missing)) {
    stop("Missing metadata IDs: ", paste(missing, collapse = ", "))
  }
  extras <- setdiff(mids, ids)
  if (length(extras) && extra == "fail") {
    stop("Extra metadata IDs: ", paste(extras, collapse = ", "), "; the table was built for a different sample set")
  }
  index <- match(ids, mids)
  aligned <- metadata[index, , drop = FALSE]
  bad <- !stats::complete.cases(aligned[, required_columns, drop = FALSE])
  if (any(bad) && na_action == "fail") {
    stop("Missing metadata values for IDs: ", paste(ids[bad], collapse = ", "))
  }
  exclusions <- data.frame(
    sample = c(extras, ids[bad]),
    reason = c(rep("extra", length(extras)), rep("missing_metadata", sum(bad)))
  )
  list(
    metadata = aligned[!bad, , drop = FALSE],
    matching = data.frame(sample = ids[!bad], matrix_row = which(!bad), metadata_row = index[!bad]),
    exclusions = exclusions
  )
}

.select_ids <- function(requested, available, label) {
  if (is.null(requested)) {
    return(seq_along(available))
  }
  if (is.factor(requested)) {
    requested <- as.character(requested)
  }
  if (!is.character(requested)) {
    stop(label, " selections must be character IDs, not positions or logicals")
  }
  .validate_ids(requested, label)
  absent <- setdiff(requested, available)
  if (length(absent)) {
    stop("Unknown ", label, " IDs: ", paste(absent, collapse = ", "))
  }
  match(requested, available)
}

.require_biological_factors <- function(ids) {
  if (!length(ids)) {
    stop("No factors: no biological factors available")
  }
  invisible(NULL)
}

.tuple_groups <- function(metadata) {
  columns <- lapply(metadata, function(column) {
    column <- as.factor(column)
    levels(column) <- as.character(seq_len(nlevels(column)))
    column
  })
  do.call(interaction, c(columns, list(drop = TRUE, lex.order = TRUE)))
}

.with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  if (
    length(seed) != 1L ||
      !is.numeric(seed) ||
      !is.finite(seed) ||
      seed < 0 ||
      seed > .Machine$integer.max ||
      seed != as.integer(seed)
  ) {
    stop("seed must be a nonnegative integer")
  }
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) {
    old <- get(".Random.seed", envir = .GlobalEnv)
  }
  kind <- RNGkind()
  on.exit(
    {
      do.call(RNGkind, as.list(kind))
      if (existed) {
        assign(".Random.seed", old, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  set.seed(seed)
  force(code)
}

# flashier accepts one of 0, 1, 2, or c(1, 2) in that order (c(2, 1) fails inside flashier);
# NULL estimates no residual variance, so it needs a supplied S.
.check_var_type <- function(var_type, S) {
  valid <- is.null(var_type) ||
    (is.numeric(var_type) && length(var_type) == 1L && var_type %in% 0:2) ||
    (is.numeric(var_type) && identical(as.numeric(var_type), c(1, 2)))
  if (!valid) {
    stop("var_type must be NULL, 0, 1, 2 or c(1, 2)", call. = FALSE)
  }
  if (is.null(var_type) && is.null(S)) stop("var_type = NULL estimates no residual variance; supply S", call. = FALSE)
}

.check_S <- function(S) {
  if (!is.null(S) && (!is.numeric(S) || any(!is.finite(S)) || any(S <= 0))) {
    stop("S must be finite and positive")
  }
  invisible(S)
}

# TRUE or FALSE, nothing else.
.is_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)

# A list whose entries, if any, carry unique names.
.is_named_list <- function(x) is.list(x) && (!length(x) || (!is.null(names(x)) && !anyDuplicated(names(x))))

#' srr standards evidence for ids.R
#'
#' @srrstats {G2.0} Length assertions exist for scalar inputs: R/programs.R .program_number(), R/ids.R .with_seed() (seed), R/fit.R fit_ebmf() (max_factors), R/ids.R .is_flag(). Coverage is by helper, not for every parameter of every function.
#' @srrstats {G2.1} Type assertions: R/ids.R .validate_matrix() (numeric matrix or double sparse Matrix), R/ids.R .check_S(), R/programs.R .program_number() (finite integer or number), R/ids.R .is_flag() (logical flag). Not every parameter is type-checked.
#' @srrstats {G2.4} Explicit conversions are used where input is coerced; see the lettered sub-items. Examples: R/ids.R .match_ids() and .select_ids(), R/batch-correction.R correct_feature_batch(), R/fit.R fit_multistart_ebmf(), R/extract.R factor_pve().
#' @srrstats {G2.4c} Character conversion uses as.character() (not paste) for metadata IDs in R/ids.R .match_ids() and for requested selections in R/ids.R .select_ids().
#' @srrstats {G2.4d} Factor conversion uses as.factor() inside R/ids.R .tuple_groups(), R/batch-correction.R correct_feature_batch() and R/batch-correction.R harmony_assignments(). It is used internally for grouping only. No user-facing conversion argument exists.
#' @srrstats {G2.4e} Factor-to-character conversion goes through as.character() in R/ids.R .select_ids() (is.factor, then as.character) before the character-only check in the same function.
#' @srrstats {G2.6} One-dimensional ID input is normalised: R/ids.R .select_ids() turns factors into character and rejects other non-character input. R/metadata.R .weights() requires numeric, named, finite, nonnegative weights. Partial: not all vector classes are handled.
#' @srrstats {G2.8} Fit input goes through one gate, R/ids.R .validate_matrix(), called from R/fit.R fit_ebmf(). Downstream functions resolve fits and views to one factor_representation class in R/representation.R .resolve_view().
#' @srrstats {G2.10} Metadata columns are taken with [[ ]] and rows with drop = FALSE in R/ids.R .match_ids(), so extraction does not depend on data.frame, tibble or data.table class.
#' @srrstats {G2.13} Missing data are checked before analysis. R/ids.R .validate_matrix() screens X for Inf/NaN and optionally NA. R/ids.R .match_ids() applies complete.cases to metadata rows. NA in X is documented as observed missingness in the R/fit.R fit_ebmf() docs.
#' @srrstats {G2.14a} Error on missing data is the default. R/ids.R .match_ids() stops on incomplete metadata when na_action = 'fail'. R/ids.R .validate_matrix() stops when allow_missing = FALSE, as R/batch-correction.R correct_feature_batch() requests.
#' @srrstats {G2.16} Inf, -Inf and NaN in X are rejected with an error in R/ids.R .validate_matrix(). Activities and weights must be finite: see R/representation.R .make_representation() and R/metadata.R .weights(). Error-only: there is no option to ignore or remove them, and the contract is not documented.
#' @srrstats {UL1.1} Input checks sit in separate helpers, .validate_ids() and .validate_matrix() in R/ids.R, called from fit_ebmf() in R/fit.R and fit_multistart_ebmf(). Messages are informative but come from untyped bare stop() calls.
#' @srrstats {UL1.2} Default id_policy = "require" makes .validate_matrix() in R/ids.R stop on missing dimnames. The "generate" path creates row and column names silently, but fit_ebmf() in R/fit.R records them in the generated_ids metadata. No message is issued.
#' @srrstats {UL3.2} fit_ebmf() in R/fit.R has id_policy = "generate", which creates row1.../col1... labels for unlabelled input in .validate_matrix() (R/ids.R). The default "require" stops on missing dimnames.
#' @noRd
NULL
