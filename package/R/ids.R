# Newly written integrity checks; matrices are never densified here.
.validate_ids <- function(ids, label) {
  if (is.null(ids) || !is.character(ids) || anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids))
    stop(label, " IDs must be unique, nonempty strings; supply explicit IDs")
  ids
}

.validate_matrix <- function(X, id_policy = "require", allow_missing = TRUE) {
  id_policy <- match.arg(id_policy, c("require", "generate"))
  sparse <- inherits(X, "sparseMatrix") && inherits(X, "dMatrix")
  if (!(is.matrix(X) && is.numeric(X)) && !sparse)
    stop("X must be a numeric matrix or a double sparse Matrix")
  if (length(dim(X)) != 2L || any(dim(X) < 1L)) stop("X requires two nonempty dimensions")
  values <- if (sparse) methods::slot(X, "x") else X
  if (any(is.infinite(values)) || any(is.nan(values))) stop("X contains Inf or NaN")
  if (!allow_missing && anyNA(values)) stop("Missing values are not supported for this operation")
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

.match_ids <- function(ids, metadata, id_col = "sample", extra = "fail",
                       na_action = "fail", required_columns = names(metadata)) {
  .validate_ids(ids, "Required sample")
  if (!is.data.frame(metadata) || !id_col %in% names(metadata)) stop("Metadata requires an explicit ", id_col, " column")
  mids <- .validate_ids(as.character(metadata[[id_col]]), "Metadata sample")
  extra <- match.arg(extra, c("fail", "drop"))
  na_action <- match.arg(na_action, c("fail", "omit"))
  if (!all(required_columns %in% names(metadata))) stop("Required metadata column is missing")
  missing <- setdiff(ids, mids)
  if (length(missing)) stop("Missing metadata IDs: ", paste(missing, collapse = ", "))
  extras <- setdiff(mids, ids)
  if (length(extras) && extra == "fail") stop("Extra metadata IDs: ", paste(extras, collapse = ", "), "; use extra='drop'")
  index <- match(ids, mids)
  aligned <- metadata[index, , drop = FALSE]
  bad <- !stats::complete.cases(aligned[, required_columns, drop = FALSE])
  if (any(bad) && na_action == "fail") stop("Missing metadata values for IDs: ", paste(ids[bad], collapse = ", "))
  exclusions <- data.frame(sample = c(extras, ids[bad]), reason = c(rep("extra", length(extras)), rep("missing_metadata", sum(bad))))
  list(metadata = aligned[!bad, , drop = FALSE],
       matching = data.frame(sample = ids[!bad], matrix_row = which(!bad), metadata_row = index[!bad]),
       exclusions = exclusions)
}

.select_ids <- function(requested, available, label) {
  if (is.null(requested)) return(seq_along(available))
  .validate_ids(requested, label)
  absent <- setdiff(requested, available)
  if (length(absent)) stop("Unknown ", label, " IDs: ", paste(absent, collapse = ", "))
  match(requested, available)
}

.with_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) || seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) stop("seed must be a nonnegative integer")
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv)
  kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (existed) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)
  force(code)
}
