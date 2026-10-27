program_fixture <- function(B, A = NULL, activity = "signed", effects = "signed") {
  if (is.null(A)) A <- matrix(seq_len(3 * ncol(B)), 3, ncol(B))
  if (is.null(rownames(B))) rownames(B) <- paste0("g", seq_len(nrow(B)))
  if (is.null(colnames(B)) || any(!nzchar(colnames(B))) || anyDuplicated(colnames(B)))
    colnames(B) <- if (ncol(B)) paste0("F", seq_len(ncol(B))) else character()
  dimnames(A) <- list(paste0("s", seq_len(nrow(A))), colnames(B))
  list(activity = A, effects = B, support = list(activity = activity, effects = effects),
    feature_scale = "simulated", preprocessing = list(scope = "fixture", transformation = "none", compatibility = "fixture"))
}

testthat::test_that("observed poles survive legal sign flips and source labels do not choose support", {
  r <- program_fixture(matrix(c(1, 0), 2, 1))
  t <- program_fixture(matrix(c(1, 0), 2, 1), activity = "nonnegative", effects = "nonnegative")
  before <- serialize(list(r, t), NULL)
  a <- compare_prior_fits(r, t)
  flipped <- r; flipped$activity <- -r$activity; flipped$effects <- -r$effects
  b <- compare_prior_fits(flipped, t)
  testthat::expect_equal(a$coverage$feasible_residual, b$coverage$feasible_residual, tolerance = 1e-10)
  testthat::expect_equal(a$coverage$feasible_residual, c(0, 0), tolerance = 1e-10)
  testthat::expect_equal(b$coverage$pole, c(-1, 1))
  r$engine <- "arbitrary-renamed-model"
  testthat::expect_equal(compare_prior_fits(r, t)$coverage, a$coverage)
  testthat::expect_identical(serialize(list(within(r, rm(engine)), t), NULL), before)
  bad <- t; bad$activity <- -bad$activity
  testthat::expect_error(compare_prior_fits(r, bad), "nonnegative")
  bad <- t; bad$effects <- -bad$effects
  testthat::expect_error(compare_prior_fits(r, bad), "nonnegative")
})

testthat::test_that("cancelling faces cannot pass on a tiny gradient alone", {
  for (eps in c(1e-9, 1e-13)) {
    B <- cbind(c(1, 0), c(-1, eps), c(0, -1))
    oracle <- c(1/eps, 1/eps, 0)
    testthat::expect_equal(drop(B %*% oracle), c(0, 1), tolerance = 1e-12)
    x <- .program_bounded_projection(B, c(0, 1), rep(TRUE, 3))
    testthat::expect_true(all(is.na(x$residual)) || sum(x$residual^2) < 1e-8)
  }
})

testthat::test_that("empty sources, roles and projection-only evidence fail safely", {
  empty <- program_fixture(matrix(numeric(), 2, 0))
  x <- compare_prior_fits(empty, empty)
  testthat::expect_equal(nrow(x$coverage), 0L)
  testthat::expect_s3_class(ggplot2::ggplot_build(plot(x, "coverage")), "ggplot_built")
  r <- program_fixture(diag(2))
  bad <- r; bad$roles <- c(F1 = "factor", F1 = "initialized")
  testthat::expect_error(compare_prior_fits(r, bad), "roles")
  bad$roles <- list(F1 = "factor", F2 = "factor")
  testthat::expect_error(compare_prior_fits(r, bad), "roles")
  s <- .program_source(r)
  h <- structure(list(metrics = data.frame(), analysis_metadata = list(input_digest = "input",
    evaluation_digest = "same", evaluation_type = "projection")), class = "holdout_result")
  ev <- list(reference = list(holdout = h, source_basis_id = s$basis_id),
    target = list(holdout = h, source_basis_id = s$basis_id))
  testthat::expect_error(compare_prior_fits(r, r, evidence = ev), "masked")
})

testthat::test_that("span and mixed-support cone have independent analytic answers", {
  r <- program_fixture(matrix(c(-1, -1), 2, 1))
  t <- program_fixture(diag(2), activity = c(F1 = "signed", F2 = "nonnegative"))
  x <- compare_prior_fits(r, t)
  row <- x$coverage[x$coverage$source == "reference", ]
  testthat::expect_equal(row$span_residual, 0, tolerance = 1e-12)
  testthat::expect_equal(row$feasible_residual, .5, tolerance = 1e-8)
  testthat::expect_lte(row$solver_kkt, 1e-7)
  t <- program_fixture(cbind(c(1, 0), c(1, 0)), activity = "nonnegative")
  r <- program_fixture(matrix(c(1, 0), 2, 1))
  x <- compare_prior_fits(r, t)
  row <- x$coverage[x$coverage$source == "reference", ]
  testthat::expect_equal(row$target_numerical_rank, 1L)
  testthat::expect_equal(row$feasible_residual, 0, tolerance = 1e-10)
  t <- program_fixture(matrix(c(1, 1, 1, -1), 2, 2))
  testthat::expect_equal(compare_prior_fits(r, t)$coverage$span_residual, c(0, .5, .5), tolerance = 1e-12)
  r <- program_fixture(matrix(c(0, 0, 1), 3, 1))
  t <- program_fixture(rbind(diag(2), c(0, 0)))
  testthat::expect_equal(compare_prior_fits(r, t)$coverage$span_residual, rep(1, 3), tolerance = 1e-12)
})

testthat::test_that("inactive targets, zero loading and unknown supports remain explicit", {
  r <- program_fixture(matrix(c(1, 0), 2, 1))
  t <- r; t$activity[,] <- 0
  x <- compare_prior_fits(r, t)
  testthat::expect_equal(x$coverage$span_residual, 1)
  testthat::expect_equal(x$coverage$target_active_factors, 0L)
  testthat::expect_equal(x$exclusions$reason, "inactive_or_zero_shared_loading")
  t <- r; t$support$activity <- "unknown"
  x <- compare_prior_fits(r, t)
  testthat::expect_true(is.na(x$coverage$feasible_residual[x$coverage$source == "reference"]))
  testthat::expect_identical(x$evidence$state, "not_assessed")
  t <- r; t$effects[,] <- 0
  testthat::expect_true(nrow(compare_prior_fits(r, t)$exclusions) == 1L)
})

testthat::test_that("mapping, units, scale and ordering contracts are checked", {
  r <- program_fixture(matrix(c(1, 2, 0, 0, 1, 2), 3, 2))
  t <- r; t$activity <- t$activity * 4; t$effects <- t$effects / 4
  expected <- compare_prior_fits(r, r)
  actual <- compare_prior_fits(r, t)
  testthat::expect_equal(actual$coverage$span_residual, expected$coverage$span_residual, tolerance = 1e-12)
  t$activity <- t$activity[, 2:1, drop = FALSE]; t$effects <- t$effects[3:1, 2:1, drop = FALSE]
  testthat::expect_equal(compare_prior_fits(r, t)$similarity_matrix[, 2:1], expected$similarity_matrix, tolerance = 1e-12)
  bad <- t; bad$feature_scale <- "other_units"
  testthat::expect_error(compare_prior_fits(r, bad), "units differ")
  testthat::expect_error(compare_prior_fits(r, r, feature_map = c(g1 = "g1", g2 = "g1")), "unique")
  bad <- r; rownames(bad$effects) <- c("g1", "g1", "g3")
  testthat::expect_error(compare_prior_fits(r, bad), "unique")
  t <- r; t$effects <- t$effects[1:2, , drop = FALSE]
  testthat::expect_equal(compare_prior_fits(r, t)$feature_coverage$fraction, c(2/3, 1))
  testthat::expect_error(compare_prior_fits(r, r, evidence = list(reference = list(stability = list()))), "exact source basis")
})

testthat::test_that("a signed direction can split into two feasible nonnegative-score directions", {
  b <- c(1, -1)
  r <- program_fixture(matrix(b, 2, 1), A = matrix(c(-2, 0, 1), 3, 1))
  t <- program_fixture(cbind(b, -b), activity = "nonnegative")
  x <- compare_prior_fits(r, t)
  testthat::expect_equal(x$coverage$span_residual, rep(0, 4), tolerance = 1e-12)
  testthat::expect_equal(x$coverage$feasible_residual, rep(0, 4), tolerance = 1e-8)
  testthat::expect_equal(x$coverage$pole[x$coverage$source == "reference"], c(-1, 1))
  testthat::expect_equal(sum(x$coverage$pole_score_fraction[x$coverage$source == "reference"]), 1)
})

testthat::test_that("supported poles and initialized roles do not invent observed usage", {
  r <- program_fixture(matrix(c(1, 0), 2, 1))
  r$roles <- c(F1 = "initialized")
  t <- program_fixture(matrix(c(1, 0), 2, 1), activity = "nonnegative", effects = "nonnegative")
  x <- compare_prior_fits(r, t, activity_poles = "supported")
  d <- x$coverage[x$coverage$source == "reference", ]
  testthat::expect_equal(d$pole, c(-1, 1))
  testthat::expect_equal(d$feasible_residual, c(1, 0), tolerance = 1e-12)
  testthat::expect_equal(d$pole_score_fraction, c(0, 1))
  testthat::expect_true(all(d$role == "initialized"))
  r$support$activity <- "unknown"
  testthat::expect_true("unknown_source_support" %in% compare_prior_fits(r, t, activity_poles = "supported")$exclusions$reason)
})

testthat::test_that("native fits and views are comparison sources", {
  X <- fixture_matrix()
  fit <- fit_ebmf(X, "rows", ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace),
    max_factors = 2L, backfit = TRUE, nullcheck = FALSE, seed = 2, verbose = 0, feature_scale = "intensity")
  view <- standardize_factors(fit)
  x <- compare_prior_fits(fit, view)
  testthat::expect_equal(x$analysis_metadata$reference$basis_id, view$manifest$basis_id)
  testthat::expect_equal(x$coverage$span_residual, rep(0, nrow(x$coverage)), tolerance = 1e-10)
  testthat::expect_true(all(diag(abs(x$similarity_matrix)) > 0.999))
  testthat::expect_equal(x$evidence$state, "not_assessed")
  for (type in c("similarity", "coverage")) testthat::expect_s3_class(ggplot2::ggplot_build(plot(x, type)), "ggplot_built")
  testthat::expect_s3_class(ggplot2::ggplot_build(plot(x, "residual", factor = "F1")), "ggplot_built")
  testthat::expect_output(print(x), "reference factors")
})

testthat::test_that("holdout evidence binds only when inputs, folds and masks are paired", {
  X <- fixture_matrix(); md <- data.frame(sample = rownames(X), donor = rep(paste0("d", 1:6), each = 2))
  folds <- data.frame(sample = md$sample, fold = c(rep("train", 8), rep("test", 4)))
  mask <- matrix(FALSE, 12, 8, dimnames = dimnames(X)); mask[9:12, c(1, 3)] <- TRUE
  callback <- function(X, metadata, sample_side, seed, settings, observation_aux = list()) {
    fit_ebmf(X, sample_side, ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace), max_factors = 2L,
      backfit = TRUE, nullcheck = FALSE, seed = seed, verbose = 0, feature_scale = "intensity")
  }
  settings <- list(sample_side = "rows", preprocessing_manifest = list(scope = "identity"), evaluated_folds = "test")
  h1 <- factor_holdout(X, md, "donor", folds, callback, evaluation = "masked", mask = mask, seed = 4, control = settings)
  h2 <- factor_holdout(X, md, "donor", folds, callback, evaluation = "masked", mask = mask, seed = 5, control = settings)
  other <- mask; other[9:12, 2] <- TRUE
  h3 <- factor_holdout(X, md, "donor", folds, callback, evaluation = "masked", mask = other, seed = 4, control = settings)
  v1 <- standardize_factors(callback(X, md, "rows", 4, settings))
  v2 <- standardize_factors(callback(X, md, "rows", 5, settings))
  evidence <- function(a, b) list(reference = list(holdout = a, source_basis_id = v1$manifest$basis_id),
    target = list(holdout = b, source_basis_id = v2$manifest$basis_id))
  x <- compare_prior_fits(v1, v2, evidence = evidence(h1, h2))
  testthat::expect_equal(x$evidence$state, "driver_bound_paired_evaluation_attached")
  testthat::expect_equal(x$evidence$prediction$input_digest, h1$analysis_metadata$input_digest)
  testthat::expect_error(compare_prior_fits(v1, v2, evidence = evidence(h1, h3)), "paired")
  unbound <- compare_prior_fits(v1, v2, evidence = list(reference = list(holdout = h1), target = list(holdout = h2)))
  testthat::expect_equal(unbound$evidence$state, "unbound_holdout_not_assessed")
})
