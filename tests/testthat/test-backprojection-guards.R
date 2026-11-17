test_that("backproject_contrast refuses undeclared or malformed coefficient inputs", {
  view <- fixture_view()
  basis <- view$manifest
  beta <- c(F1 = 2, F2 = -1)
  run <- function(x, ...) backproject_contrast(view, x, "coefficients", factor_basis = basis, ...)
  expect_error(backproject_contrast(view, beta, factor_basis = basis), "input must explicitly declare")
  expect_error(run(beta, feature_block_size = 0), "feature_block_size")
  expect_error(run(beta, level = 1), "level must be between")
  expect_error(
    backproject_contrast(view, beta, "coefficients", factor_basis = list(fit_id = "x")),
    "Supply factor_basis"
  )
  inconsistent <- basis
  inconsistent$feature_scale <- "other_units"
  expect_error(backproject_contrast(view, beta, "coefficients", factor_basis = inconsistent), "inconsistent")
  expect_error(run(c(F1 = 2, F2 = NA)), "finite, complete")
  expect_error(run(c(F1 = 2, F2 = -1, F9 = 1)), "Unknown coefficient factors")
})

test_that("backproject_contrast selects rows of a coefficient coefficients explicitly", {
  view <- fixture_view()
  basis <- view$manifest
  coefficients <- data.frame(
    factor = rep(c("F1", "F2"), 2),
    contrast = rep(c("up", "down"), each = 2),
    estimate = c(2, -1, -2, 1)
  )
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
  run <- function(Sigma) {
    backproject_contrast(view, beta, "coefficients", factor_basis = basis, coefficient_covariance = Sigma)
  }
  expect_error(run(list(diag(2), diag(2))), "one full covariance per contrast")
  expect_error(run(diag(3)), "finite K by K")
  negative <- matrix(c(1, 2, 2, 1), 2, dimnames = list(c("F1", "F2"), c("F1", "F2")))
  expect_error(run(negative), "positive semidefinite")
  Sigma <- matrix(c(1, 0.5, 0.5, 2), 2, dimnames = list(c("F2", "F1"), c("F2", "F1")))
  out <- run(Sigma)
  B <- view$effects[, c("F1", "F2")]
  expect_equal(
    out$effects$std.error,
    sqrt(rowSums((B %*% Sigma[c("F1", "F2"), c("F1", "F2")]) * B)),
    ignore_attr = TRUE
  )
})

test_that("backproject_contrast draw mode refuses ambiguous or partial draws", {
  view <- fixture_view()
  basis <- view$manifest
  draws <- matrix(c(2, 2.5, -1, -0.5), 2, dimnames = list(NULL, c("F1", "F2")))
  run <- function(x, ...) backproject_contrast(view, x, "draws", factor_basis = basis, ...)
  expect_error(run(draws, contrasts = c("a", "b")), "one selected contrast")
  expect_error(run(data.frame(draw = 1, factor = "F1")), "requires draw, factor and estimate")
  long <- data.frame(
    draw = rep(1:2, each = 2),
    factor = rep(c("F1", "F2"), 2),
    estimate = c(2, -1, 2.5, -0.5),
    term = rep(c("a", "b"), 2)
  )
  expect_error(run(long), "Explicit term selection")
  expect_error(run(rbind(long[long$term == "a", ], long[1, ]), term = "a"), "unique complete")
  expect_error(run(draws[1, , drop = FALSE]), "at least two draws")
  expect_error(run(draws[, "F1", drop = FALSE]), "allow_partial is required")
  partial <- run(draws[, "F1", drop = FALSE], allow_partial = TRUE)
  expect_identical(partial$analysis_metadata$omitted_factors, "F2")
})

test_that("a known df gives a t interval, and the table df column is used when none is given", {
  view <- fixture_view()
  beta <- c(F1 = 2, F2 = -1)
  Sigma <- matrix(c(0.2, 0.05, 0.05, 0.1), 2, dimnames = list(names(beta), names(beta)))
  base <- backproject_contrast(view, beta, "coefficients", factor_basis = view$manifest, coefficient_covariance = Sigma)
  expect_equal(
    base$effects$upper - base$effects$estimate,
    stats::qnorm(0.975) * base$effects$std.error
  )
  expect_match(base$analysis_metadata$uncertainty_type, "^normal")
  with_df <- backproject_contrast(
    view,
    beta,
    "coefficients",
    factor_basis = view$manifest,
    coefficient_covariance = Sigma,
    df = 4
  )
  expect_equal(with_df$effects$upper - with_df$effects$estimate, stats::qt(0.975, 4) * with_df$effects$std.error)
  expect_match(with_df$analysis_metadata$uncertainty_type, "^t_")
  expect_equal(with_df$analysis_metadata$interval_df, 4)
  tab <- data.frame(factor = names(beta), estimate = unname(beta), df = c(10, 4))
  from_table <- backproject_contrast(
    view,
    tab,
    "coefficients",
    factor_basis = view$manifest,
    coefficient_covariance = Sigma
  )
  expect_equal(from_table$analysis_metadata$interval_df, 4)
  expect_error(
    backproject_contrast(view, beta, "coefficients", factor_basis = view$manifest, df = -1),
    "df must be one positive number"
  )
})

test_that("an intercept term is labelled a contribution and warns when the offset is nonzero", {
  view <- fixture_view()
  view <- .make_representation(
    view$activity,
    view$effects,
    view$manifest,
    offset = stats::setNames(rep(1, nrow(view$effects)), rownames(view$effects))
  )
  beta <- c(F1 = 2, F2 = -1)
  expect_warning(
    out <- backproject_contrast(
      view,
      beta,
      "coefficients",
      contrasts = "(Intercept)",
      factor_basis = view$manifest
    ),
    "not a contrast"
  )
  expect_identical(out$analysis_metadata$estimand, "factor_contribution_excluding_offset")
  contrast <- backproject_contrast(view, beta, "coefficients", contrasts = "groupB", factor_basis = view$manifest)
  expect_identical(contrast$analysis_metadata$estimand, "represented_feature_contrast")
})

test_that("a table whose term and contrast columns match needs one selection", {
  d <- association_fixture()
  K <- ncol(d$view$activity)
  set.seed(2)
  effects <- matrix(rnorm(6 * K), 6, K, dimnames = list(paste0("g", 1:6), colnames(d$view$activity)))
  view <- .make_representation(d$view$activity, effects, d$view$manifest)
  out <- test_factor_association(
    view,
    d$metadata,
    ~group,
    "lm",
    unit_col = "donor",
    independence = "independent",
    include_intercept = TRUE
  )
  by_term <- backproject_contrast(view, out$table, "coefficients", term = "groupB")
  by_contrast <- backproject_contrast(view, out$table, "coefficients", contrasts = "groupB")
  expect_equal(by_term$effects$estimate, by_contrast$effects$estimate)
  # oracle: the back-projected contrast equals a per-feature lm on the fitted matrix
  fitted <- view$activity %*% t(view$effects)
  direct <- apply(fitted, 2L, function(y) unname(stats::coef(stats::lm(y ~ group, data = d$metadata))[["groupB"]]))
  expect_equal(by_term$effects$estimate, unname(direct))
})

test_that("term and conditions are refused for a vector or matrix, and a vector takes one contrast", {
  view <- fixture_view()
  beta <- c(F1 = 2, F2 = -1)
  expect_error(
    backproject_contrast(view, beta, "coefficients", term = "x", factor_basis = view$manifest),
    "select rows of a coefficient table"
  )
  expect_error(
    backproject_contrast(view, beta, "coefficients", conditions = "x", factor_basis = view$manifest),
    "select rows of a coefficient table"
  )
  expect_error(
    backproject_contrast(view, beta, "coefficients", contrasts = c("a", "b"), factor_basis = view$manifest),
    "one contrast label"
  )
  draws <- matrix(rnorm(20 * 2), 20, dimnames = list(NULL, names(beta)))
  expect_error(
    backproject_contrast(view, draws, "draws", term = "x", factor_basis = view$manifest),
    "select rows of a coefficient table"
  )
})
