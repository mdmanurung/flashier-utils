test_that("prediction R squared is negative even with perfect correlation", {
  observed <- data.frame(feature = c("g1", "g2", "g3"), estimate = c(1, 2, 3), contrast = "A - B")
  reconstructed <- observed
  reconstructed$estimate <- c(11, 12, 13)
  out <- validate_backprojection(
    reconstructed,
    observed,
    feature_scale = "arcsinh_intensity",
    top_n = 20L,
    plot = FALSE
  )
  expect_equal(out$metrics$value[out$metrics$metric == "prediction_r2"], -149)
  expect_equal(out$metrics$value[out$metrics$metric == "correlation_squared"], 1)
  expect_equal(out$top_n_overlap$effective_n, 3)
  expect_equal(out$residuals$residual, c(-10, -10, -10))
  attr(reconstructed, "analysis_metadata") <- list(feature_scale = "log1p")
  expect_error(validate_backprojection(reconstructed, observed, feature_scale = "arcsinh_intensity"), "scales")
  attr(reconstructed, "analysis_metadata") <- NULL
  reconstructed$contrast <- "B - A"
  expect_error(validate_backprojection(reconstructed, observed, feature_scale = "arcsinh_intensity"), "Contrast")
  observed$estimate <- 1
  out <- validate_backprojection(observed, observed, feature_scale = "arcsinh_intensity", plot = FALSE)
  expect_true(is.na(out$metrics$value[out$metrics$metric == "prediction_r2"]))
})

test_that("estimates constant up to rounding give NA, not a correlation read from noise", {
  features <- paste0("g", 1:6)
  observed <- data.frame(feature = features, estimate = c(0.1 + 0.2, rep(0.3, 5)), contrast = "A - B")
  reconstructed <- data.frame(feature = features, estimate = c(1, 2, 3, 4, 5, 6), contrast = "A - B")
  out <- validate_backprojection(reconstructed, observed, feature_scale = "arcsinh_intensity", plot = FALSE)
  value <- function(name) out$metrics$value[out$metrics$metric == name]
  expect_true(is.na(value("pearson_correlation")))
  expect_true(is.na(value("prediction_r2")))
  expect_true(is.na(value("slope")))
})

test_that("noise in the observed table attenuates slope but not calibration_slope", {
  features <- paste0("g", 1:60)
  truth <- seq(-3, 3, length.out = 60)
  set.seed(1)
  noise <- stats::resid(stats::lm(rnorm(60) ~ truth))
  observed <- data.frame(feature = features, estimate = truth + 0.5 * noise, contrast = "A - B")
  reconstructed <- data.frame(feature = features, estimate = truth, contrast = "A - B")
  out <- validate_backprojection(reconstructed, observed, feature_scale = "arcsinh_intensity", plot = FALSE)
  value <- function(name) out$metrics$value[out$metrics$metric == name]
  # with noise orthogonal to the truth: slope = var(truth) / (var(truth) + var(noise)); calibration slope = 1
  expect_equal(value("calibration_slope"), 1)
  expect_equal(value("slope"), stats::var(truth) / (stats::var(truth) + 0.25 * stats::var(noise)))
  expect_lt(value("slope"), 1)
})

test_that("metrics equal hand-computed values and do not depend on row order", {
  observed <- data.frame(feature = paste0("g", 1:4), estimate = c(1, 2, 3, 4), contrast = "A - B")
  reconstructed <- data.frame(feature = paste0("g", 1:4), estimate = c(1, 3, 2, 5), contrast = "A - B")
  shuffled <- reconstructed[c(3, 1, 4, 2), ]
  out <- validate_backprojection(shuffled, observed, feature_scale = "arcsinh_intensity", top_n = 2L, plot = FALSE)
  value <- function(name) out$metrics$value[out$metrics$metric == name]
  # sxx = 5, syy = 8.75, sxy = 5.5 (means 2.5 and 2.75)
  expect_equal(value("slope"), 5.5 / 5)
  expect_equal(value("intercept"), 2.75 - 1.1 * 2.5)
  expect_equal(value("calibration_slope"), 5.5 / 8.75)
  expect_equal(value("pearson_correlation"), 5.5 / sqrt(5 * 8.75))
  expect_equal(value("sign_concordance"), 1)
  # top 2 by absolute observed: g4, g3; by reconstructed: g4, g2
  expect_equal(out$top_n_overlap$overlap_n, 1)
  expect_equal(out$top_n_overlap$overlap_fraction, 0.5)
})
