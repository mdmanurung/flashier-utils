test_that("selection preserves dimensions, order and matched basis", {
  fit <- fixture_fit()
  view <- standardize_factors(fit)
  A <- factor_activity(view, factors = "F2", samples = c("s3", "s1"))
  expect_equal(dim(A), c(2, 1))
  expect_equal(rownames(A), c("s3", "s1"))
  expect_equal(attr(A, "analysis_metadata")$basis_id, view$manifest$basis_id)
  B <- factor_effects(view, factors = "F2", features = "g3", format = "long")
  expect_equal(B$estimate, view$effects["g3", "F2"])
  expect_equal(B$basis_id, view$manifest$basis_id)
  expect_error(factor_activity(view, samples = "unknown"), "Unknown")
  expect_error(factor_activity(view, factors = c("F1", "F1")), "unique")
})

test_that("factor_pve names its denominator, and a feature intercept stays in it", {
  plain <- factor_pve(fixture_fit())
  expect_identical(unique(plain$denominator), "second_moment")
  with_intercept <- fit_ebmf(
    fixture_matrix(),
    "rows",
    ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace),
    max_factors = 2L,
    intercept = "feature",
    nullcheck = FALSE,
    backfit = TRUE,
    verbose = 0,
    seed = 12
  )
  expect_identical(unique(factor_pve(with_intercept)$denominator), "second_moment_including_intercept")
  expect_identical(
    unique(factor_pve(standardize_factors(with_intercept))$denominator),
    "second_moment_including_intercept"
  )
})
