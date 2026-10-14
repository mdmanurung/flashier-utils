test_that("zero-factor summaries preserve their schemas and sample population", {
  view <- fixture_view()
  empty <- .make_representation(view$activity[, FALSE, drop = FALSE],
    view$effects[, FALSE, drop = FALSE],
    utils::modifyList(view$manifest, list(native_indices = integer())))
  metadata <- data.frame(sample = rownames(view$activity), group = c("A", "A", "B"), age = 1:3)
  summary <- summarize_factor_activity(empty, metadata, "group")
  expect_s3_class(summary, "data.frame")
  expect_equal(nrow(summary), 0L)
  expect_identical(summary$factor, character())
  expect_type(summary$mean, "double")
  expect_equal(nrow(rank_factors_by_metadata(empty, metadata, "age")), 0L)
  diagnostics <- factor_diagnostics(empty)
  expect_identical(diagnostics$factors$factor, character())
  expect_equal(diagnostics$samples$programs, rep(0, 3))
  expect_equal(nrow(summarize_factors(empty)), 0L)
  expect_equal(nrow(rank_features(empty)), 0L)
  measurements <- cbind(marker = c(s1 = 1, s2 = 2, s3 = 3))
  correlations <- correlate_factor_markers(empty, measurements)
  expect_s3_class(correlations, "data.frame")
  expect_equal(nrow(correlations), 0L)
  expect_type(correlations$correlation, "double")
  expect_error(compare_factor_groups(empty, metadata, "group", c("A", "B")), "No factors")
  expect_error(plot_factor_activity(empty), "No factors")
  if (requireNamespace("pROC", quietly = TRUE)) {
    specificity <- factor_specificity(empty, metadata, "group")
    expect_s3_class(specificity, "data.frame")
    expect_equal(nrow(specificity), 0L)
    expect_type(specificity$auc, "double")
  }
})

test_that("named native zero-rank fits need no sampler for empty intervals", {
  fit <- fit_ebmf(fixture_matrix(), "rows", max_factors = 0, verbose = 0)
  result <- factor_intervals(fit, "effects", nsamp = 2)
  expect_identical(result$factor, character())
  expect_identical(result$lower, numeric())
})
