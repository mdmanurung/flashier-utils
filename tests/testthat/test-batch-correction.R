# Fixture: two cell states, three batches; batch shifts differ by state.
batch_fixture <- function() {
  .with_seed(11, {
    n <- 120
    state <- rep(c("alpha", "beta"), each = n / 2)
    batch <- rep(c("b1", "b2", "b3"), length.out = n)
    condition <- rep(c("ctrl", "case"), times = c(70, 50))
    shift <- outer(batch == "b2" & state == "beta", c(3, 0, 1, 2)) +
      outer(batch == "b3" & state == "alpha", c(0, 2, 0, 1))
    X <- matrix(rnorm(n * 4), n, 4) + outer(state == "beta", c(4, 1, 0, 2)) + shift
    dimnames(X) <- list(paste0("s", seq_len(n)), paste0("g", 1:4))
    metadata <- data.frame(sample = rownames(X), batch = batch, state = state, condition = condition)
    list(X = X, metadata = metadata)
  })
}

batch_means <- function(X, metadata, state) {
  in_state <- metadata$state == state
  rowsum(X[in_state, , drop = FALSE], metadata$batch[in_state]) / as.vector(table(metadata$batch[in_state]))
}

test_that("one cluster with a vanishing penalty equals limma::removeBatchEffect", {
  skip_if_not_installed("limma")
  fx <- batch_fixture()
  one_cluster <- matrix(1, nrow(fx$X), 1, dimnames = list(rownames(fx$X), "all"))
  corrected <- correct_feature_batch(fx$X, "rows", fx$metadata, "batch", one_cluster, lambda = 1e-6)
  expected <- t(limma::removeBatchEffect(t(fx$X), batch = fx$metadata$batch, design = matrix(1, nrow(fx$X), 1)))
  expect_equal(corrected, expected, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("a protected covariate matches removeBatchEffect with that design", {
  skip_if_not_installed("limma")
  fx <- batch_fixture()
  one_cluster <- matrix(1, nrow(fx$X), 1, dimnames = list(rownames(fx$X), "all"))
  corrected <- correct_feature_batch(
    fx$X,
    "rows",
    fx$metadata,
    "batch",
    one_cluster,
    protect = ~condition,
    lambda = 1e-6
  )
  design <- stats::model.matrix(~condition, fx$metadata)
  expected <- t(limma::removeBatchEffect(t(fx$X), batch = fx$metadata$batch, design = design))
  expect_equal(corrected, expected, tolerance = 1e-5, ignore_attr = TRUE)
  expect_identical(attr(corrected, "analysis_metadata")$protected_columns, "conditionctrl")
})

test_that("with a vanishing penalty, hard labels equalize state-specific batch means that global centering leaves", {
  fx <- batch_fixture()
  corrected <- correct_feature_batch(fx$X, "rows", fx$metadata, "batch", "state", lambda = 1e-6)
  one_cluster <- matrix(1, nrow(fx$X), 1, dimnames = list(rownames(fx$X), "all"))
  global <- correct_feature_batch(fx$X, "rows", fx$metadata, "batch", one_cluster, lambda = 1e-6)
  for (state in c("alpha", "beta")) {
    means <- batch_means(corrected, fx$metadata, state)
    expect_lt(max(abs(sweep(means, 2, means[1, ]))), 1e-4)
  }
  global_means <- batch_means(global, fx$metadata, "beta")
  expect_gt(max(abs(sweep(global_means, 2, global_means[1, ]))), 1)
})

test_that("per-state means are kept: cell-state differences survive", {
  fx <- batch_fixture()
  corrected <- correct_feature_batch(fx$X, "rows", fx$metadata, "batch", "state", lambda = 1e-6)
  gap_after <- colMeans(corrected[fx$metadata$state == "beta", ]) - colMeans(corrected[fx$metadata$state == "alpha", ])
  # unweighted mean over batches of the per-batch state means, as the ridge limit keeps
  beta_mean <- colMeans(batch_means(fx$X, fx$metadata, "beta"))
  alpha_mean <- colMeans(batch_means(fx$X, fx$metadata, "alpha"))
  expect_equal(gap_after, beta_mean - alpha_mean, tolerance = 1e-4)
})

test_that("sample_side columns returns the transpose of the rows result", {
  fx <- batch_fixture()
  by_rows <- correct_feature_batch(fx$X, "rows", fx$metadata, "batch", "state")
  by_columns <- correct_feature_batch(t(fx$X), "columns", fx$metadata, "batch", "state")
  expect_equal(by_columns, t(by_rows), ignore_attr = TRUE)
  expect_identical(dimnames(by_columns), rev(dimnames(fx$X)))
})

test_that("a cluster with fewer than two batches is left unchanged and recorded", {
  fx <- batch_fixture()
  fx$metadata$state[fx$metadata$batch == "b1"] <- "only_b1"
  expect_warning(
    corrected <- correct_feature_batch(fx$X, "rows", fx$metadata, "batch", "state"),
    "left unchanged: only_b1"
  )
  untouched <- fx$metadata$state == "only_b1"
  expect_identical(unclass(corrected)[untouched, ], fx$X[untouched, ])
  expect_identical(attr(corrected, "analysis_metadata")$skipped_clusters, "only_b1")
})

test_that("the automatic penalty is alpha times the expected batch counts", {
  fx <- batch_fixture()
  corrected <- correct_feature_batch(fx$X, "rows", fx$metadata, "batch", "state", alpha = 0.5)
  penalty <- attr(corrected, "analysis_metadata")$penalty
  # 60 samples per state, 40 per batch out of 120
  expect_equal(penalty, matrix(0.5 * 60 * 40 / 120, 2, 3, dimnames = list(c("alpha", "beta"), c("b1", "b2", "b3"))))
})

test_that("metadata rows are matched by ID, not position", {
  fx <- batch_fixture()
  shuffled <- fx$metadata[rev(seq_len(nrow(fx$metadata))), ]
  expect_equal(
    correct_feature_batch(fx$X, "rows", shuffled, "batch", "state"),
    correct_feature_batch(fx$X, "rows", fx$metadata, "batch", "state"),
    ignore_attr = TRUE
  )
})

test_that("invalid inputs stop before any correction", {
  fx <- batch_fixture()
  soft <- matrix(0.4, nrow(fx$X), 2, dimnames = list(rownames(fx$X), c("H1", "H2")))
  expect_error(correct_feature_batch(fx$X, "rows", fx$metadata, "batch", soft), "summing to one")
  one_batch <- transform(fx$metadata, batch = "b1")
  expect_error(correct_feature_batch(fx$X, "rows", one_batch, "batch", "state"), "at least two batches")
  expect_error(correct_feature_batch(fx$X, "rows", fx$metadata, "batch", "state", lambda = -1), "positive")
  wrong_names <- matrix(1, 2, 3, dimnames = list(c("alpha", "beta"), c("x", "y", "z")))
  expect_error(correct_feature_batch(fx$X, "rows", fx$metadata, "batch", "state", lambda = wrong_names), "batch levels")
  expect_error(correct_feature_batch(fx$X, "rows", fx$metadata[-1, ], "batch", "state"), "Missing metadata IDs")
  expect_error(correct_feature_batch(fx$X, "rows", fx$metadata, "batch", "state", alpha = 0), "alpha")
  with_missing <- fx$X
  with_missing[1, 1] <- NA
  expect_error(correct_feature_batch(with_missing, "rows", fx$metadata, "batch", "state"), "Missing values")
})

test_that("a protected covariate constant within a cluster is an error", {
  fx <- batch_fixture()
  fx$metadata$condition <- ifelse(fx$metadata$state == "beta", "case", "ctrl")
  expect_error(
    correct_feature_batch(fx$X, "rows", fx$metadata, "batch", "state", protect = ~condition),
    "rank deficient within the corrected samples of cluster alpha"
  )
})

test_that("applied to Harmony's input, the correction reproduces Harmony's embedding", {
  skip_if_not_installed("harmony")
  fx <- batch_fixture()
  pcs <- fx$X
  colnames(pcs) <- paste0("PC", 1:4)
  h <- harmony_assignments(pcs, fx$metadata, "batch", seed = 3, nclust = 4)
  expect_equal(unname(rowSums(h$assignments)), rep(1, nrow(pcs)), tolerance = 1e-5)
  expect_identical(colnames(h$lambda), c("b1", "b2", "b3"))
  corrected <- correct_feature_batch(pcs, "rows", fx$metadata, "batch", h$assignments, lambda = h$lambda)
  # Harmony works in single precision, so agreement is to about 1e-5 absolute
  expect_equal(corrected, h$embedding, tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("harmony_assignments is reproducible under a seed", {
  skip_if_not_installed("harmony")
  fx <- batch_fixture()
  pcs <- fx$X
  colnames(pcs) <- paste0("PC", 1:4)
  first <- harmony_assignments(pcs, fx$metadata, "batch", seed = 3, nclust = 4)
  second <- harmony_assignments(pcs, fx$metadata, "batch", seed = 3, nclust = 4)
  expect_identical(first$assignments, second$assignments)
})

test_that("harmony_assignments drops batch levels with no samples", {
  skip_if_not_installed("harmony")
  fx <- batch_fixture()
  pcs <- fx$X
  colnames(pcs) <- paste0("PC", 1:4)
  fx$metadata$batch <- factor(fx$metadata$batch, levels = c("b0", "b1", "b2", "b3"))
  h <- harmony_assignments(pcs, fx$metadata, "batch", seed = 3, nclust = 4)
  expect_identical(colnames(h$lambda), c("b1", "b2", "b3"))
})

test_that("assignments equal to the batch labels warn that nothing was corrected", {
  fx <- batch_fixture()
  expect_warning(
    corrected <- correct_feature_batch(fx$X, "rows", fx$metadata, "batch", "batch"),
    "nothing was corrected"
  )
  expect_identical(unclass(corrected)[,], fx$X[,])
})

test_that("a protected covariate fully determined by batch warns and keeps the shift", {
  fx <- batch_fixture()
  fx$metadata$nested <- ifelse(fx$metadata$batch == "b3", "case", "ctrl")
  one_cluster <- matrix(1, nrow(fx$X), 1, dimnames = list(rownames(fx$X), "all"))
  expect_warning(
    corrected <- correct_feature_batch(
      fx$X,
      "rows",
      fx$metadata,
      "batch",
      one_cluster,
      protect = ~nested,
      lambda = 1e-6
    ),
    "fully determined by batch"
  )
  gap <- function(x) colMeans(x[fx$metadata$batch == "b3", ]) - colMeans(x[fx$metadata$batch != "b3", ])
  expect_equal(gap(corrected), gap(fx$X), tolerance = 1e-3, ignore_attr = TRUE)
})

test_that("a protected factor level absent from one cluster is dropped there, not an error", {
  fx <- batch_fixture()
  fx$metadata$group <- rep(c("ctrl", "case"), length.out = nrow(fx$metadata))
  fx$metadata$group[fx$metadata$state == "beta" & seq_len(nrow(fx$metadata)) %% 4 == 0] <- "other"
  corrected <- correct_feature_batch(fx$X, "rows", fx$metadata, "batch", "state", protect = ~group, lambda = 1e-6)
  expect_equal(dim(corrected), dim(fx$X))
})

test_that("a protect formula that gives NaN stops instead of dropping rows", {
  fx <- batch_fixture()
  fx$metadata$v <- rnorm(nrow(fx$metadata))
  fx$metadata$v[5] <- -100
  expect_error(
    suppressWarnings(
      correct_feature_batch(fx$X, "rows", fx$metadata, "batch", "state", protect = ~ log(v + 60))
    ),
    "missing values"
  )
})

test_that("Harmony parity holds with unbalanced batches and a cluster that drops a batch", {
  skip_if_not_installed("harmony")
  set.seed(5)
  n <- 300
  batch <- sample(c("run1", "run2", "run3"), n, replace = TRUE, prob = c(0.6, 0.3, 0.1))
  pcs <- matrix(rnorm(n * 4), n, 4, dimnames = list(paste0("s", seq_len(n)), paste0("PC", 1:4))) +
    outer(batch == "run2", c(1.5, 0, 1, 0))
  metadata <- data.frame(sample = rownames(pcs), batch = batch)
  h <- harmony_assignments(pcs, metadata, "batch", seed = 3, nclust = 8)
  # unbalanced batches make the automatic penalty differ across batches within a cluster
  expect_gt(max(apply(h$lambda, 1, function(x) diff(range(x)))), 0.1)
  corrected <- suppressWarnings(
    correct_feature_batch(pcs, "rows", metadata, "batch", h$assignments, lambda = h$lambda)
  )
  expect_equal(corrected, h$embedding, tolerance = 1e-3, ignore_attr = TRUE)
})
