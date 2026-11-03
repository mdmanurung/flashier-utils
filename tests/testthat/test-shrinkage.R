shrinkage_coefficients <- function(conditions = c("A", "B", "C")) {
  d <- association_fixture(K = 8)
  coefficients <- expand.grid(factor = colnames(d$view$activity), cohort = conditions, stringsAsFactors = FALSE)
  coefficients$estimate <- .with_seed(6, rnorm(nrow(coefficients)))
  coefficients$std.error <- rep(0.5, nrow(coefficients))
  coefficients$term <- "groupB"
  coefficients$basis_id <- d$view$manifest$basis_id
  attr(coefficients, "analysis_metadata") <- d$view$manifest
  coefficients
}

test_that("shrink_factor_effects refuses inputs it cannot shrink", {
  skip_if_not_installed("mashr")
  coefficients <- shrinkage_coefficients()
  run <- function(x, ...) shrink_factor_effects(x, "cohort", term = "groupB", independent_conditions = TRUE, ...)
  expect_error(shrink_factor_effects(coefficients, "cohort", method = "ash", term = "groupB"), "Choose mash")
  expect_error(shrink_factor_effects(coefficients, "cohort"), "Choose mash")
  expect_error(run(unclass(coefficients)), "coefficient table")
  expect_error(shrink_factor_effects(coefficients, "cohort", term = "other", independent_conditions = TRUE), "one comparable term")
  expect_error(run(shrinkage_coefficients("A")), "at least two conditions")
  bad <- coefficients
  bad$estimate[1] <- NA_real_
  expect_error(run(bad), "finite estimates")
  bad <- coefficients
  bad$std.error[1] <- 0
  expect_error(run(bad), "positive SEs")
  bad <- rbind(coefficients, coefficients[1, ])
  attr(bad, "analysis_metadata") <- attr(coefficients, "analysis_metadata")
  expect_error(run(bad), "Duplicate factor/condition")
  expect_error(run(coefficients, control = list(nonsense = 1)), "Unsupported native mash")
  expect_error(run(coefficients, control = list(outputlevel = 1L)), "outputlevel")
  expect_error(run(coefficients, covariances = list()), "covariances must be")
})

test_that("shrink_factor_effects validates and reorders a supplied V", {
  skip_if_not_installed("mashr")
  coefficients <- shrinkage_coefficients()
  run <- function(V) shrink_factor_effects(coefficients, "cohort", term = "groupB", V = V)
  expect_error(run(diag(2)), "finite named condition correlation")
  asymmetric <- diag(3)
  asymmetric[1, 2] <- 0.5
  dimnames(asymmetric) <- list(c("A", "B", "C"), c("A", "B", "C"))
  expect_error(run(asymmetric), "symmetric positive definite")
  V <- matrix(c(1, 0.3, 0, 0.3, 1, 0, 0, 0, 1), 3, dimnames = list(c("B", "A", "C"), c("B", "A", "C")))
  out <- run(V)
  expect_identical(dimnames(out$analysis_metadata$error_correlation), list(c("A", "B", "C"), c("A", "B", "C")))
  expect_equal(out$analysis_metadata$error_correlation["A", "B"], 0.3)
})

test_that("shrink_factor_effects accepts a result list and a caller covariance list", {
  skip_if_not_installed("mashr")
  coefficients <- shrinkage_coefficients()
  wrapped <- list(table = coefficients, analysis_metadata = attr(coefficients, "analysis_metadata"))
  attr(wrapped$table, "analysis_metadata") <- NULL
  U <- list(identity = diag(3))
  out <- shrink_factor_effects(wrapped, "cohort", term = "groupB", independent_conditions = TRUE, covariances = U)
  expect_identical(out$covariances, U)
  expect_identical(out$analysis_metadata$covariance_family, "caller_supplied")
  native <- mashr::mash(mashr::mash_set_data(out$Bhat, out$Shat), U, verbose = FALSE)
  expect_equal(out$table$estimate, as.vector(ashr::get_pm(native)))
})
