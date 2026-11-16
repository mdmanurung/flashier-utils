#' srr standards evidence for the tests
#'
#' @srrstats {G5.4a} Methods are checked against trivial and independent cases: a known rotation of a rank-2 view (rotated rank-2 test in tests/testthat/test-fu1-stability.R) and a from-scratch flashier pipeline (tests/testthat/test-fnd-oracles.R). The correlation path is checked against stats::cor in tests/testthat/test-correlation.R.
#' @srrstats {G5.4b} The package tests against earlier implementations: pinned singlecelljamboreeR helpers (tests/testthat/helper-upstream.R, used in tests/testthat/test-features.R) and stats::cor (tests/testthat/test-correlation.R; R/correlation.R factor_cor()). Coverage is limited to these functions.
#' @srrstats {G5.6} A simulated truth is recovered. tests/testthat/test-fu3a-rna-adt.R builds ADT effects from known activities and compares projected effects with that oracle. tests/testthat/test-fu1-stability.R recovers a rotated rank-2 subspace through alignment. Coverage targets projection and alignment, not every function.
#' @srrstats {G5.6a} Recovery tests use tolerances, not exact equality. tests/testthat/test-fu3a-rna-adt.R requires per-protein correlation above 0.98 and RMSE below 0.15. tests/testthat/test-fnd-oracles.R compares native and wrapped results with tolerance = 1e-12.
#' @srrstats {G5.8a} tests/testthat/test-srr-edges.R: fit_ebmf() refuses matrices with zero rows or zero columns ('X requires two nonempty dimensions').
#' @srrstats {G5.8b} tests/testthat/test-srr-edges.R: character, complex and data frame input are refused ('X must be a numeric matrix or a double sparse Matrix').
#' @srrstats {G5.8d} tests/testthat/test-srr-edges.R: a matrix with more features than samples fits, standardizes and reconstructs the fitted values.
#' @srrstats {G5.9} tests/testthat/test-srr-edges.R covers noise of 1e-10 (G5.9a) and different seeds (G5.9b).
#' @srrstats {G5.9a} tests/testthat/test-srr-edges.R: noise of sd 1e-10 added to the fixture matrix changes the fitted L_pm by less than 1e-6.
#' @srrstats {G5.9b} tests/testthat/test-srr-edges.R: fits with seeds 12 and 13 match both factors by factor_replicability() at congruence 0.99.
#' @srrstats {G5.12} tests/README.md describes how to run the tests, the fixture sizes, the optional-engine skips and the approximate runtime.
#' @noRd
NULL

# Edge cases for rOpenSci standards G5.8 (input extremes) and G5.9 (noise and seeds).
srr_fit <- function(X, seed = 12, sample_side = "rows") {
  fit_ebmf(
    X,
    sample_side,
    ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace),
    max_factors = 2L,
    nullcheck = FALSE,
    backfit = TRUE,
    verbose = 0L,
    seed = seed
  )
}

test_that("a matrix with a zero-extent dimension is refused", {
  empty_rows <- matrix(numeric(), 0, 2, dimnames = list(NULL, c("g1", "g2")))
  expect_error(fit_ebmf(empty_rows, "rows"), "nonempty")
  empty_cols <- matrix(numeric(), 3, 0, dimnames = list(paste0("s", 1:3), NULL))
  expect_error(fit_ebmf(empty_cols, "rows"), "nonempty")
})

test_that("non-numeric matrices are refused even with dimnames", {
  X <- fixture_matrix()
  expect_error(fit_ebmf(matrix(as.character(X), 12, 8, dimnames = dimnames(X)), "rows"), "numeric matrix")
  expect_error(fit_ebmf(matrix(complex(real = X), 12, 8, dimnames = dimnames(X)), "rows"), "numeric matrix")
  expect_error(fit_ebmf(as.data.frame(X), "rows"), "numeric matrix")
})

test_that("more features than samples fits and standardizes", {
  X <- t(fixture_matrix())
  expect_gt(ncol(X), nrow(X))
  view <- standardize_factors(srr_fit(X))
  expect_identical(rownames(view$activity), rownames(X))
  expect_identical(rownames(view$effects), colnames(X))
  expect_equal(unname(tcrossprod(view$activity, view$effects)), unname(stats::fitted(srr_fit(X))), tolerance = 1e-8)
})

test_that("noise of 1e-10 changes the fitted loadings by less than 1e-6", {
  X <- fixture_matrix()
  noisy <- X + .with_seed(3, matrix(rnorm(length(X), sd = 1e-10), nrow(X), ncol(X)))
  expect_equal(srr_fit(noisy)$L_pm, srr_fit(X)$L_pm, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("different seeds recover the same two factors", {
  skip_if_not_installed("psych")
  skip_if_not_installed("clue")
  v12 <- standardize_factors(srr_fit(fixture_matrix(), seed = 12), feature_scale = "intensity")
  v13 <- standardize_factors(srr_fit(fixture_matrix(), seed = 13), feature_scale = "intensity")
  expect_equal(sum(factor_replicability(v12, v13, min_congruence = 0.99)$matches$matched), 2)
})
