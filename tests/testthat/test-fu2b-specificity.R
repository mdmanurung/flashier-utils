# FU-2B: grouped specificity validation — missingness, ties, unit overlap,
# hand/pROC oracles, class imbalance. No API change; existing return columns
# are pinned. Thresholds/direction are learned on the same labels: the result
# is in-sample descriptive, never predictive (analysis_metadata$inference).

# 5-sample fixture with a hand-computable Mann-Whitney AUC.
fu2b_view <- function() {
  A <- matrix(c(1, 2.5, 2.5, 4, 0.5), 5, 1,
              dimnames = list(paste0("s", 1:5), "F1"))
  B <- matrix(1:4, 4, 1, dimnames = list(paste0("g", 1:4), "F1"))
  .make_representation(A, B, list(fit_id = "fu2b_fixture", schema_version = "1.0.0",
      sample_side = "rows", native_indices = 1L, scaling = "raw",
      orientation = "as_fit", constraints = c(activity = "signed", effects = "signed"),
      feature_scale = "arcsinh_intensity"))
}

test_that("NA group labels fail with an explicit message, not an any(NA) exception", {
  skip_if_not_installed("pROC")
  v <- fu2b_view()
  md <- data.frame(sample = paste0("s", 1:5), group = c("pos", "neg", NA, "pos", "neg"))
  expect_error(factor_specificity(v, md, "group"), "Missing metadata values")
  md$group[3] <- "pos"
  md$unit <- c("u1", "u1", "u2", NA, "u2")
  expect_error(factor_specificity(v, md, "group", unit_col = "unit"), "Missing metadata values")
})

test_that("AUC matches a hand-computed Mann-Whitney oracle, both directions", {
  skip_if_not_installed("pROC")
  v <- fu2b_view()
  md <- data.frame(sample = paste0("s", 1:5),
                   group = c("pos", "neg", "pos", "pos", "neg"))
  # pos activities 1, 2.5, 4 vs neg 2.5, 0.5: W=4, T=1, L=1 -> (4 + 0.5)/6
  r <- factor_specificity(v, md, "group", groups = "pos")
  expect_equal(r$auc[1], 0.75)
  rn <- factor_specificity(v, md, "group", groups = "pos", direction = "negative")
  expect_equal(rn$auc[1], 0.25)
  # pROC agrees with the hand oracle (engine cross-check)
  roc <- pROC::roc(md$group == "pos", v$activity[, 1], levels = c(FALSE, TRUE),
                   direction = "<", quiet = TRUE)
  expect_equal(as.numeric(pROC::auc(roc)), 0.75)
})

test_that("threshold ties are counted and the smallest finite threshold wins", {
  skip_if_not_installed("pROC")
  A <- matrix(c(1, 2, 2, 3), 4, 1, dimnames = list(paste0("s", 1:4), "F1"))
  B <- matrix(1:4, 4, 1, dimnames = list(paste0("g", 1:4), "F1"))
  v <- .make_representation(A, B, list(fit_id = "fu2b_tie", schema_version = "1.0.0",
      sample_side = "rows", native_indices = 1L, scaling = "raw",
      orientation = "as_fit", constraints = c(activity = "signed", effects = "signed"),
      feature_scale = "arcsinh_intensity"))
  md <- data.frame(sample = paste0("s", 1:4), group = c("neg", "neg", "pos", "pos"))
  r <- factor_specificity(v, md, "group", groups = "pos")
  roc <- pROC::roc(md$group == "pos", c(1, 2, 2, 3), levels = c(FALSE, TRUE),
                   direction = "<", quiet = TRUE)
  best <- pROC::coords(roc, "best", best.method = "closest.topleft", transpose = FALSE)
  expect_equal(r$threshold_ties[1], nrow(best))
  expect_gt(r$threshold_ties[1], 1)  # this fixture genuinely ties
  expect_equal(r$threshold[1], min(best$threshold[is.finite(best$threshold)]))
})

test_that("units spanning both classes are visibly double-counted; counts stay descriptive", {
  skip_if_not_installed("pROC")
  A <- matrix(c(1, 4, 2, 3, 0.5, 3.5), 6, 1, dimnames = list(paste0("s", 1:6), "F1"))
  B <- matrix(1:4, 4, 1, dimnames = list(paste0("g", 1:4), "F1"))
  v <- .make_representation(A, B, list(fit_id = "fu2b_units", schema_version = "1.0.0",
      sample_side = "rows", native_indices = 1L, scaling = "raw",
      orientation = "as_fit", constraints = c(activity = "signed", effects = "signed"),
      feature_scale = "arcsinh_intensity"))
  # u1 and u2 each span both classes -> not independent groups
  md <- data.frame(sample = paste0("s", 1:6),
                   group = c("pos", "neg", "pos", "neg", "neg", "pos"),
                   unit = c("u1", "u1", "u2", "u2", "u3", "u3"))
  r <- factor_specificity(v, md, "group", groups = "pos", unit_col = "unit")
  expect_equal(r$n_units[1], 3)
  expect_equal(r$n_positive_units[1], 3)    # u1, u2, u3 all have a positive sample
  expect_equal(r$n_negative_units[1], 3)    # u1, u2, u3 all have a negative sample
  expect_gt(r$n_positive_units[1] + r$n_negative_units[1], r$n_units[1])  # overlap visible
  expect_identical(attr(r, "analysis_metadata")$inference, "descriptive")
})

test_that("extreme class imbalance stays typed and hand-checkable", {
  skip_if_not_installed("pROC")
  n <- 20
  A <- matrix(c(25, seq_len(n - 1)), n, 1, dimnames = list(paste0("s", 1:n), "F1"))
  B <- matrix(1:4, 4, 1, dimnames = list(paste0("g", 1:4), "F1"))
  v <- .make_representation(A, B, list(fit_id = "fu2b_imb", schema_version = "1.0.0",
      sample_side = "rows", native_indices = 1L, scaling = "raw",
      orientation = "as_fit", constraints = c(activity = "signed", effects = "signed"),
      feature_scale = "arcsinh_intensity"))
  md <- data.frame(sample = paste0("s", 1:n), group = c("pos", rep("neg", n - 1)))
  r <- factor_specificity(v, md, "group", groups = "pos")
  expect_equal(r$status[1], "available")
  expect_equal(r$n_positive[1], 1)
  expect_equal(r$n_negative[1], n - 1)
  expect_equal(r$auc[1], 1)  # the single positive outranks every negative
})

test_that("return schema is preserved exactly (no silent column drift)", {
  skip_if_not_installed("pROC")
  v <- fu2b_view()
  md <- data.frame(sample = paste0("s", 1:5), group = c("pos", "neg", "pos", "pos", "neg"))
  r <- factor_specificity(v, md, "group", groups = "pos", unit_col = NULL)
  expect_named(r, c("factor", "group", "auc", "threshold", "sensitivity",
                    "specificity", "threshold_ties", "n_samples", "n_positive",
                    "n_negative", "n_units", "n_positive_units", "n_negative_units",
                    "status"))
})
