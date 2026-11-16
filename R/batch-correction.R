# Harmony's per-cluster ridge correction, rewritten in R for a feature matrix.
# Algorithm: Korsunsky et al. (2019); gene-matrix use: Kotliar's cNMF Preprocess.
# harmony.cpp (2.1.0) was read to confirm the update order; no code is copied.
# Parity is tested against the installed harmony package (2.0.5 when checked).

#' Remove batch effects from a feature matrix, cluster by cluster
#'
#' Removes batch shifts from a samples by features matrix before
#' [fit_ebmf()], using Harmony's mixture-of-experts ridge correction. Each
#' sample is softly assigned to clusters of similar cell state. Within each
#' cluster, every feature is regressed on batch with a ridge penalty, weighted
#' by the assignments, and the fitted batch terms are subtracted. Batch shifts
#' may therefore differ between cell states, which a single per-batch centering
#' cannot model. The input is not clipped, so signed data such as shifted log
#' expression stay signed. Nonnegative input can gain negative values, and
#' zeros in sparse input become nonzero.
#' @param X Named numeric matrix or double sparse Matrix with no missing values,
#'   for example shifted log expression. The result is dense.
#' @param sample_side `"rows"` or `"columns"` of `X`.
#' @param metadata Data frame with one row per sample, keyed by `id_col`. Extra
#'   rows are dropped and recorded.
#' @param batch_col Metadata column holding the batch to remove; at least two
#'   levels. Factor level order is ignored beyond the order of the penalty columns.
#' @param assignments Either a samples by clusters matrix of soft assignments,
#'   with sample row names and cluster column names, nonnegative and each row
#'   summing to one within 1e-5 (as from [harmony_assignments()]); or the name
#'   of a metadata column of hard cluster labels, such as cell types.
#' @param protect Optional one-sided formula of metadata covariates to keep,
#'   such as `~ condition`. They enter each cluster's regression without
#'   penalty and are not subtracted. A covariate fully determined by batch
#'   within a cluster (for example a condition applied to whole batches) keeps
#'   the batch shift aligned with it; the function warns in that case.
#' @param id_col Metadata column holding sample IDs.
#' @param lambda Ridge penalty on the batch terms. `NULL` (default) uses
#'   Harmony's automatic penalty, `alpha` times the expected number of samples
#'   of each batch in each cluster. A positive number applies one penalty
#'   everywhere. A clusters by batches matrix, named like `assignments` columns
#'   and the batch levels, gives each its own penalty.
#' @param alpha Positive scale of the automatic penalty (Harmony default 0.2).
#'   Used only when `lambda` is `NULL`. The ridge leaves part of each batch
#'   shift in place: about `alpha * E / (O + alpha * E)`, where `O` is the
#'   batch's weighted sample count in the cluster and `E` its expected count.
#'   In one test with a planted shift of 2 the default left a gap of 0.33.
#'   Batches with fewer samples than expected in a cluster are corrected less.
#'   Use a small `lambda` such as 1e-6 for near-complete removal.
#' @param batch_prop_cutoff A batch is corrected within a cluster only when its
#'   mean assignment to that cluster exceeds this value (Harmony default 1e-5).
#'   A cluster with fewer than two such batches is left unchanged.
#' @return The corrected matrix, with the same orientation and IDs as `X`. Its
#'   `analysis_metadata` attribute records the batch levels, the penalty matrix
#'   used, skipped clusters, protected columns, dropped metadata rows and
#'   SHA-256 hashes of the input and assignments.
#' @details For cluster \eqn{k} with assignment weights \eqn{r}, design
#'   \eqn{\Phi} (intercept, protected covariates, batch indicators) and penalty
#'   \eqn{\Lambda_k} (zero for the intercept and protected columns), the
#'   coefficients are \eqn{W_k = (\Phi^T \mathrm{diag}(r) \Phi + \Lambda_k)^{-1}
#'   \Phi^T \mathrm{diag}(r) X}, always computed from the uncorrected `X`. The
#'   batch rows of \eqn{\mathrm{diag}(r) \Phi W_k} are subtracted, and the
#'   corrections add up over clusters. Without `protect` this reproduces
#'   Harmony's own correction of its input embedding.
#'
#'   Batch here should be technical. If `batch_col` holds donors, and donors
#'   are nested in a condition, the condition's effect within each cell state
#'   is removed too unless `protect` names it. A protected covariate that is
#'   constant within a cluster's corrected samples is an error.
#'
#'   The assignments are learned from the same data, so the correction is part
#'   of the fit. For [factor_stability()] or [factor_holdout()], run it inside
#'   the `fit_fn` callback on each resample; correcting once beforehand
#'   overstates stability. Memory: the result is a dense n by p double matrix,
#'   and `sample_side = "columns"` holds one transposed copy as well.
#'
#'   The automatic penalty (`alpha`) and `batch_prop_cutoff` come from the Harmony
#'   2.x software, not from the 2019 paper.
#' @references Korsunsky I, Millard N, Fan J, et al. (2019). Fast, sensitive and
#'   accurate integration of single-cell data with Harmony. *Nature Methods*
#'   16, 1289-1296. \doi{10.1038/s41592-019-0619-0}
#'
#'   Kotliar D, Veres A, Nagy MA, et al. (2019). Identifying gene expression
#'   programs of cell-type identity and cellular activity with single-cell
#'   RNA-Seq. *eLife* 8, e43803. \doi{10.7554/eLife.43803}
#' @seealso [harmony_assignments()] for soft assignments from an embedding.
#' @author Mikhael Manurung (R implementation). Harmony algorithm: Ilya
#'   Korsunsky, Soumya Raychaudhuri and coauthors; gene-matrix adaptation:
#'   Dylan Kotliar (cNMF). No upstream code is copied.
#' @export
#' @examples
#' set.seed(1)
#' state <- rep(c("alpha", "beta"), each = 20)
#' batch <- rep(c("b1", "b2"), times = 20)
#' X <- matrix(rnorm(40 * 6), 40, 6) + outer(batch == "b2" & state == "beta", rep(2, 6))
#' dimnames(X) <- list(paste0("s", 1:40), paste0("g", 1:6))
#' metadata <- data.frame(sample = rownames(X), batch = batch, state = state)
#' corrected <- correct_feature_batch(X, "rows", metadata, "batch", assignments = "state")
#' # batch b2 minus b1 within the beta state, before and after
#' beta <- state == "beta"
#' colMeans(X[beta & batch == "b2", ]) - colMeans(X[beta & batch == "b1", ])
#' colMeans(corrected[beta & batch == "b2", ]) - colMeans(corrected[beta & batch == "b1", ])
correct_feature_batch <- function(
  X,
  sample_side,
  metadata,
  batch_col,
  assignments,
  protect = NULL,
  id_col = "sample",
  lambda = NULL,
  alpha = 0.2,
  batch_prop_cutoff = 1e-5
) {
  sample_side <- match.arg(sample_side, c("rows", "columns"))
  X <- .validate_matrix(X, allow_missing = FALSE)
  input_sha256 <- digest::digest(X, algo = "sha256")
  X <- as.matrix(X)
  if (sample_side == "columns") {
    X <- t(X)
  }
  samples <- rownames(X)
  if (!is.character(batch_col) || length(batch_col) != 1L) {
    stop("batch_col must name one metadata column")
  }
  if (!is.null(protect) && !inherits(protect, "formula")) {
    stop("protect must be NULL or a one-sided formula such as ~ condition")
  }
  hard_labels <- is.character(assignments) && length(assignments) == 1L
  required <- unique(c(id_col, batch_col, if (hard_labels) assignments, if (!is.null(protect)) all.vars(protect)))
  matched <- .match_ids(samples, metadata, id_col = id_col, extra = "drop", required_columns = required)
  aligned <- matched$metadata

  # batch indicators in as.factor() level order, as Harmony builds them
  batch <- droplevels(as.factor(aligned[[batch_col]]))
  if (nlevels(batch) < 2L) {
    stop("batch_col needs at least two batches")
  }
  onehot <- outer(as.integer(batch), seq_len(nlevels(batch)), "==") * 1
  colnames(onehot) <- levels(batch)
  batch_sizes <- colSums(onehot)

  # soft assignments: a matrix keyed by sample, or one-hot hard labels
  if (hard_labels) {
    cluster_labels <- droplevels(as.factor(aligned[[assignments]]))
    cluster_weights <- outer(as.integer(cluster_labels), seq_len(nlevels(cluster_labels)), "==") * 1
    colnames(cluster_weights) <- levels(cluster_labels)
    assignment_source <- paste0("metadata column ", assignments)
  } else {
    if (!is.matrix(assignments) || !is.numeric(assignments)) {
      stop("assignments must be a samples by clusters numeric matrix or one metadata column name")
    }
    assignments <- .validate_matrix(assignments, allow_missing = FALSE)
    if (!setequal(rownames(assignments), samples)) {
      stop("assignments row names must be exactly the samples of X")
    }
    cluster_weights <- assignments[samples, , drop = FALSE]
    if (any(cluster_weights < 0) || any(abs(rowSums(cluster_weights) - 1) > 1e-5)) {
      stop("assignments must be nonnegative with each row summing to one")
    }
    assignment_source <- "matrix"
  }
  rownames(cluster_weights) <- samples
  clusters <- colnames(cluster_weights)

  # protected covariates enter without penalty and are never subtracted
  protected <- if (is.null(protect)) {
    matrix(0, length(samples), 0L)
  } else {
    # na.fail: a formula that gives NA or NaN must stop, not drop rows by position
    design <- stats::model.matrix(protect, stats::model.frame(protect, data = aligned, na.action = stats::na.fail))
    design[, colnames(design) != "(Intercept)", drop = FALSE]
  }

  # penalty per cluster and batch
  expected <- outer(colSums(cluster_weights), batch_sizes / length(samples))
  dimnames(expected) <- list(clusters, levels(batch))
  penalty <- if (is.null(lambda)) {
    if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) {
      stop("alpha must be one positive number")
    }
    alpha * expected
  } else if (is.matrix(lambda)) {
    if (!identical(rownames(lambda), clusters) || !identical(colnames(lambda), levels(batch))) {
      stop("lambda matrix needs row names equal to the clusters and column names equal to the batch levels")
    }
    if (!is.numeric(lambda) || anyNA(lambda) || any(lambda <= 0)) {
      stop("lambda must be positive")
    }
    lambda
  } else {
    if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) || lambda <= 0) {
      stop("lambda must be NULL, one positive number or a clusters by batches matrix")
    }
    expected * 0 + lambda
  }
  .program_number(batch_prop_cutoff, "batch_prop_cutoff")

  # per-cluster ridge from the uncorrected matrix; corrections add up
  corrected <- X
  skipped <- character()
  for (k in clusters) {
    r <- cluster_weights[, k]
    kept <- which(colSums(onehot * r) / batch_sizes > batch_prop_cutoff)
    if (length(kept) < 2L) {
      skipped <- c(skipped, k)
      next
    }
    rows <- which(rowSums(onehot[, kept, drop = FALSE]) > 0)
    # a protected column with no weight in this cluster (an absent factor level) is dropped here
    protected_k <- protected[rows, colSums(abs(protected[rows, , drop = FALSE]) * r[rows]) > 0, drop = FALSE]
    n_unpenalized <- 1L + ncol(protected_k)
    design_k <- cbind(1, protected_k, onehot[rows, kept, drop = FALSE])
    weighted <- design_k * r[rows]
    unpenalized <- design_k[, seq_len(n_unpenalized), drop = FALSE] * sqrt(r[rows])
    if (qr(unpenalized)$rank < n_unpenalized) {
      stop(
        "Protected covariates are rank deficient within the corrected samples of cluster ",
        k,
        " (constant, collinear with each other, or an intercept-free factor coding)"
      )
    }
    # a protected column fully explained by batch keeps the batch shift aligned with it
    if (ncol(protected_k) > 0L) {
      batch_part <- onehot[rows, kept, drop = FALSE] * sqrt(r[rows])
      protected_part <- protected_k * sqrt(r[rows])
      residual <- qr.resid(qr(batch_part), protected_part)
      nested <- colnames(protected_k)[colSums(residual^2) <= 1e-8 * pmax(colSums(protected_part^2), 1e-12)]
      if (length(nested) > 0L) {
        warning(
          "Protected covariate(s) ",
          paste(nested, collapse = ", "),
          " are fully determined by batch in cluster ",
          k,
          "; the batch shift aligned with them is kept, not removed",
          call. = FALSE
        )
      }
    }
    ridge <- diag(c(rep(0, n_unpenalized), penalty[k, kept]), ncol(design_k))
    coefficients <- solve(crossprod(weighted, design_k) + ridge, crossprod(weighted, X[rows, , drop = FALSE]))
    coefficients[seq_len(n_unpenalized), ] <- 0
    corrected[rows, ] <- corrected[rows, , drop = FALSE] - weighted %*% coefficients
  }
  if (length(skipped) == length(clusters)) {
    warning(
      "Every cluster has fewer than two batches above batch_prop_cutoff, so nothing was corrected. ",
      "Check that the assignments are not the batch labels.",
      call. = FALSE
    )
  } else if (length(skipped) > 0L) {
    warning(
      length(skipped),
      " cluster(s) with fewer than two batches above batch_prop_cutoff were left unchanged: ",
      paste(utils::head(skipped, 5L), collapse = ", "),
      if (length(skipped) > 5L) ", ...",
      call. = FALSE
    )
  }
  if (sample_side == "columns") {
    corrected <- t(corrected)
  }

  attr(corrected, "analysis_metadata") <- list(
    sample_side = sample_side,
    batch_col = batch_col,
    batch_levels = levels(batch),
    assignment_source = assignment_source,
    cluster_ids = clusters,
    skipped_clusters = skipped,
    penalty = penalty,
    penalty_source = if (is.null(lambda)) "automatic_alpha_times_expected" else "supplied",
    alpha = if (is.null(lambda)) alpha else NA_real_,
    batch_prop_cutoff = batch_prop_cutoff,
    protected_columns = colnames(protected),
    exclusions = matched$exclusions,
    input_sha256 = input_sha256,
    assignments_sha256 = digest::digest(cluster_weights, algo = "sha256")
  )
  .attach_provenance(
    corrected,
    "correct_feature_batch",
    match.call(),
    .resolved_parameters("correct_feature_batch", environment())
  )
}

#' Soft cluster assignments from Harmony
#'
#' Runs [harmony::RunHarmony()] on a low-dimensional embedding, such as
#' principal components of shifted log expression, and returns its soft cluster
#' assignments and penalties for [correct_feature_batch()]. Needs the optional
#' `harmony` package.
#' @param embedding Named numeric samples by dimensions matrix, such as PCs.
#' @param metadata Data frame with one row per sample, keyed by `id_col`.
#' @param batch_col Metadata column holding the batch; at least two levels.
#' @param seed Nonnegative integer. Required; Harmony's k-means start uses R's
#'   random numbers. The caller's random number state is restored afterwards.
#' @param id_col Metadata column holding sample IDs.
#' @param ... Further arguments to [harmony::RunHarmony()], such as `nclust`,
#'   `theta` or `.options`.
#' @return A list with `$assignments` (samples by clusters, rows summing to
#'   one), `$lambda` (clusters by batch levels, the penalties Harmony used),
#'   `$embedding` (Harmony's corrected embedding) and `$analysis_metadata`.
#' @details Pass `$assignments` and `$lambda` to [correct_feature_batch()] to
#'   apply Harmony's correction to the full feature matrix. Applied to
#'   `embedding` itself, that reproduces `$embedding` to Harmony's
#'   single-precision arithmetic. Only one batch variable is supported.
#' @references Korsunsky I, Millard N, Fan J, et al. (2019). Fast, sensitive and
#'   accurate integration of single-cell data with Harmony. *Nature Methods*
#'   16, 1289-1296. \doi{10.1038/s41592-019-0619-0}
#' @author Mikhael Manurung (adapter). harmony: Ilya Korsunsky, Nikolaos
#'   Patikas, Soumya Raychaudhuri and contributors listed in its DESCRIPTION.
#' @export
#' @examplesIf requireNamespace("harmony", quietly = TRUE)
#' set.seed(1)
#' batch <- rep(c("b1", "b2"), each = 50)
#' pcs <- matrix(rnorm(100 * 5), 100, 5) + outer(batch == "b2", rep(1, 5))
#' dimnames(pcs) <- list(paste0("s", 1:100), paste0("PC", 1:5))
#' metadata <- data.frame(sample = rownames(pcs), batch = batch)
#' h <- harmony_assignments(pcs, metadata, "batch", seed = 1, nclust = 3)
#' dim(h$assignments)
harmony_assignments <- function(embedding, metadata, batch_col, seed, id_col = "sample", ...) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("harmony_assignments() needs the harmony package")
  }
  if (!is.matrix(embedding) || !is.numeric(embedding)) {
    stop("embedding must be a named numeric samples by dimensions matrix")
  }
  embedding <- .validate_matrix(embedding, allow_missing = FALSE)
  .program_number(seed, "seed", integer = TRUE)
  if (!is.character(batch_col) || length(batch_col) != 1L) {
    stop("batch_col must name one metadata column")
  }
  samples <- rownames(embedding)
  matched <- .match_ids(samples, metadata, id_col = id_col, extra = "drop", required_columns = c(id_col, batch_col))
  aligned <- matched$metadata
  # Harmony builds one batch per factor level, so drop levels with no samples
  aligned[[batch_col]] <- droplevels(as.factor(aligned[[batch_col]]))
  batch_levels <- levels(aligned[[batch_col]])
  if (length(batch_levels) < 2L) {
    stop("batch_col needs at least two batches")
  }

  # run Harmony and keep its fitted object
  fitted <- .with_seed(seed, {
    harmony::RunHarmony(
      embedding,
      meta_data = aligned,
      vars_use = batch_col,
      return_object = TRUE,
      verbose = FALSE,
      ...
    )
  })

  # assignments, penalties and corrected embedding, named by sample and cluster
  clusters <- paste0("H", seq_len(fitted$K))
  assignments <- t(fitted$getR())
  dimnames(assignments) <- list(samples, clusters)
  lambda <- fitted$getLambda()[, -1L, drop = FALSE]
  dimnames(lambda) <- list(clusters, batch_levels)
  corrected <- t(fitted$getZcorr())
  dimnames(corrected) <- dimnames(embedding)

  out <- list(
    assignments = assignments,
    lambda = lambda,
    embedding = corrected,
    analysis_metadata = list(
      engine = "harmony",
      engine_version = as.character(utils::packageVersion("harmony")),
      batch_col = batch_col,
      batch_levels = batch_levels,
      n_clusters = fitted$K,
      exclusions = matched$exclusions,
      seed = seed
    )
  )
  .attach_provenance(
    out,
    "harmony_assignments",
    match.call(),
    .resolved_parameters("harmony_assignments", environment())
  )
}
