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
  expect_error(
    shrink_factor_effects(coefficients, "cohort", term = "other", independent_conditions = TRUE),
    "one comparable term"
  )
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

test_that("a df column is passed to mash and matches a manual mashr run", {
  skip_if_not_installed("mashr")
  coefficients <- shrinkage_coefficients()
  coefficients$df <- 4
  run <- function(x) shrink_factor_effects(x, "cohort", term = "groupB", independent_conditions = TRUE)
  with_df <- run(coefficients)
  coefficients$df <- NULL
  without_df <- run(coefficients)
  ids <- rownames(with_df$Bhat)
  manual_data <- mashr::mash_set_data(with_df$Bhat, with_df$Shat, df = 4)
  manual <- mashr::mash(manual_data, Ulist = mashr::cov_canonical(manual_data), verbose = FALSE)
  expect_equal(unname(as.vector(ashr::get_lfsr(manual))), with_df$table$lfsr)
  expect_identical(with_df$analysis_metadata$df_source, "table_df_column")
  expect_identical(without_df$analysis_metadata$df_source, "infinite_assumed")
  # few residual df make the same estimates less certain, so lfsr cannot be smaller on the whole
  expect_gt(median(with_df$table$lfsr / without_df$table$lfsr), 1)
})

test_that("an invalid df column is refused", {
  skip_if_not_installed("mashr")
  coefficients <- shrinkage_coefficients()
  coefficients$df <- 0
  expect_error(
    shrink_factor_effects(coefficients, "cohort", term = "groupB", independent_conditions = TRUE),
    "df column"
  )
})

test_that("one N(0, I) component with a unit grid gives the closed-form posterior", {
  skip_if_not_installed("mashr")
  m <- association_fixture(K = 8)$view$manifest
  set.seed(4)
  Bhat <- matrix(rnorm(24, sd = 2), 8, 3, dimnames = list(m$factor_ids, c("A", "B", "C")))
  Shat <- matrix(runif(24, 0.3, 1.5), 8, 3, dimnames = dimnames(Bhat))
  coefficients <- expand.grid(factor = rownames(Bhat), cohort = colnames(Bhat), stringsAsFactors = FALSE)
  coefficients$estimate <- as.vector(Bhat)
  coefficients$std.error <- as.vector(Shat)
  coefficients$term <- "groupB"
  coefficients$basis_id <- m$basis_id
  attr(coefficients, "analysis_metadata") <- m
  # one covariance and one grid value make mixsqp report its trivial solution; that notice is expected here
  closed_form <- function(x) {
    suppressWarnings(shrink_factor_effects(
      x,
      "cohort",
      term = "groupB",
      independent_conditions = TRUE,
      covariances = list(id = diag(3)),
      control = list(grid = 1, usepointmass = FALSE)
    ))$table
  }
  out <- closed_form(coefficients)
  posterior_mean <- Bhat / (1 + Shat^2)
  posterior_sd <- sqrt(Shat^2 / (1 + Shat^2))
  expect_equal(out$estimate, as.vector(posterior_mean))
  expect_equal(out$posterior.sd, as.vector(posterior_sd))
  expect_equal(out$lfsr, as.vector(pnorm(-abs(posterior_mean) / posterior_sd)))
  # row order of the input table does not change the answer
  shuffled <- coefficients[.with_seed(9, sample(nrow(coefficients))), ]
  attr(shuffled, "analysis_metadata") <- m
  again <- closed_form(shuffled)
  # conditions follow first appearance in the table, so compare by factor and condition
  key <- function(x) x[order(x$factor, x$condition), ]
  expect_equal(key(again)$estimate, key(out)$estimate)
  expect_equal(key(again)$lfsr, key(out)$lfsr)
})

test_that("a different df per cell is matched by factor and condition, not by row order", {
  skip_if_not_installed("mashr")
  coefficients <- shrinkage_coefficients()
  # df varies by cell and includes Inf; the table is shuffled so row order carries no meaning
  coefficients$df <- rep(c(3, 5, 8, 20, Inf, 12), length.out = nrow(coefficients))
  shuffled <- coefficients[.with_seed(11, sample(nrow(coefficients))), ]
  run <- function(x) shrink_factor_effects(x, "cohort", term = "groupB", independent_conditions = TRUE)
  out <- run(shuffled)
  ids <- rownames(out$Bhat)
  conditions <- colnames(out$Bhat)
  df_matrix <- matrix(NA_real_, length(ids), length(conditions), dimnames = list(ids, conditions))
  df_matrix[cbind(match(coefficients$factor, ids), match(coefficients$cohort, conditions))] <- coefficients$df
  expect_equal(out$df, df_matrix)
  manual_fit <- function(Dhat) {
    data <- mashr::mash_set_data(out$Bhat, out$Shat, df = Dhat)
    mashr::mash(data, Ulist = mashr::cov_canonical(data), verbose = FALSE)
  }
  manual <- manual_fit(df_matrix)
  expect_equal(out$table$estimate, unname(as.vector(ashr::get_pm(manual))), tolerance = 1e-10)
  expect_equal(out$table$lfsr, unname(as.vector(ashr::get_lfsr(manual))), tolerance = 1e-10)
  # cells with infinite df keep their input SE
  expect_equal(out$Shat_used[is.infinite(df_matrix)], out$Shat[is.infinite(df_matrix)])
  # a mis-indexed df must change the answer
  reversed <- shuffled
  reversed$df <- rev(reversed$df)
  expect_gt(max(abs(run(reversed)$table$lfsr - out$table$lfsr)), 1e-6)
})

test_that("very large t statistics with finite df run and give finite lfsr", {
  skip_if_not_installed("mashr")
  coefficients <- shrinkage_coefficients()
  coefficients$df <- 4
  # t = 1e8 still converts in mashr; t = 1e150 underflows pt() and made the old path fail
  coefficients$estimate[1] <- 1e8 * coefficients$std.error[1]
  coefficients$estimate[2] <- 1e150 * coefficients$std.error[2]
  # mash warns about non-finite mixture likelihoods at this extreme; the fit still completes
  out <- suppressWarnings(
    shrink_factor_effects(coefficients, "cohort", term = "groupB", independent_conditions = TRUE)
  )
  expect_true(all(is.finite(out$table$lfsr)))
  expect_true(all(is.finite(out$Shat_used)))
})

test_that("an exactly zero estimate with finite df stops and names the cell", {
  skip_if_not_installed("mashr")
  coefficients <- shrinkage_coefficients()
  coefficients$df <- 4
  coefficients$estimate[2] <- 0
  expect_error(
    shrink_factor_effects(coefficients, "cohort", term = "groupB", independent_conditions = TRUE),
    sprintf("exactly 0.*'%s'.*'%s'", coefficients$factor[2], coefficients$cohort[2])
  )
  # with infinite df the zero is harmless
  coefficients$df <- Inf
  expect_no_error(shrink_factor_effects(coefficients, "cohort", term = "groupB", independent_conditions = TRUE))
})
