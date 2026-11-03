test_that("backproject_contrast refuses undeclared or malformed coefficient inputs", {
  view <- fixture_view()
  basis <- view$manifest
  beta <- c(F1 = 2, F2 = -1)
  run <- function(x, ...) backproject_contrast(view, x, "coefficients", factor_basis = basis, ...)
  expect_error(backproject_contrast(view, beta, factor_basis = basis), "input must explicitly declare")
  expect_error(run(beta, feature_block_size = 0), "feature_block_size")
  expect_error(run(beta, level = 1), "level must be between")
  expect_error(backproject_contrast(view, beta, "coefficients", factor_basis = list(fit_id = "x")), "Supply factor_basis")
  inconsistent <- basis
  inconsistent$feature_scale <- "other_units"
  expect_error(backproject_contrast(view, beta, "coefficients", factor_basis = inconsistent), "inconsistent")
  expect_error(run(c(F1 = 2, F2 = NA)), "finite, complete")
  expect_error(run(c(F1 = 2, F2 = -1, F9 = 1)), "Unknown coefficient factors")
})

test_that("backproject_contrast selects rows of a coefficient coefficients explicitly", {
  view <- fixture_view()
  basis <- view$manifest
  coefficients <- data.frame(factor = rep(c("F1", "F2"), 2), contrast = rep(c("up", "down"), each = 2), estimate = c(2, -1, -2, 1))
  run <- function(x, ...) backproject_contrast(view, x, "coefficients", factor_basis = basis, ...)
  expect_error(run(coefficients[, c("factor", "contrast")]), "requires factor and estimate")
  expect_error(run(coefficients), "Explicit contrast selection")
  expect_error(run(coefficients, contrasts = "up", term = "x"), "no term column")
  expect_error(run(coefficients, contrasts = "sideways"), "No coefficients match")
  expect_error(run(rbind(coefficients, coefficients[1, ]), contrasts = "up"), "Duplicate factor/contrast")
  up <- run(coefficients, contrasts = "up")
  expect_equal(up$effects$estimate, c(1, 0, -4, 2))
  both <- run(coefficients, contrasts = c("up", "down"))
  expect_equal(both$effects$estimate[both$effects$contrast == "down"], -c(1, 0, -4, 2))
})

test_that("backproject_contrast checks the full coefficient covariance", {
  view <- fixture_view()
  basis <- view$manifest
  beta <- c(F1 = 2, F2 = -1)
  run <- function(Sigma) backproject_contrast(view, beta, "coefficients", factor_basis = basis, coefficient_covariance = Sigma)
  expect_error(run(list(diag(2), diag(2))), "one full covariance per contrast")
  expect_error(run(diag(3)), "finite K by K")
  negative <- matrix(c(1, 2, 2, 1), 2, dimnames = list(c("F1", "F2"), c("F1", "F2")))
  expect_error(run(negative), "positive semidefinite")
  Sigma <- matrix(c(1, 0.5, 0.5, 2), 2, dimnames = list(c("F2", "F1"), c("F2", "F1")))
  out <- run(Sigma)
  B <- view$effects[, c("F1", "F2")]
  expect_equal(out$effects$std.error, sqrt(rowSums((B %*% Sigma[c("F1", "F2"), c("F1", "F2")]) * B)), ignore_attr = TRUE)
})

test_that("backproject_contrast draw mode refuses ambiguous or partial draws", {
  view <- fixture_view()
  basis <- view$manifest
  draws <- matrix(c(2, 2.5, -1, -0.5), 2, dimnames = list(NULL, c("F1", "F2")))
  run <- function(x, ...) backproject_contrast(view, x, "draws", factor_basis = basis, ...)
  expect_error(run(draws, contrasts = c("a", "b")), "one selected contrast")
  expect_error(run(data.frame(draw = 1, factor = "F1")), "requires draw, factor and estimate")
  long <- data.frame(draw = rep(1:2, each = 2), factor = rep(c("F1", "F2"), 2), estimate = c(2, -1, 2.5, -0.5), term = rep(c("a", "b"), 2))
  expect_error(run(long), "Explicit term selection")
  expect_error(run(rbind(long[long$term == "a", ], long[1, ]), term = "a"), "unique complete")
  expect_error(run(draws[1, , drop = FALSE]), "at least two draws")
  expect_error(run(draws[, "F1", drop = FALSE]), "allow_partial is required")
  partial <- run(draws[, "F1", drop = FALSE], allow_partial = TRUE)
  expect_identical(partial$analysis_metadata$omitted_factors, "F2")
})
