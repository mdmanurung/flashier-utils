test_that("semi-nonnegative fitting matches native calls in both orientations", {
  A <- matrix(seq_len(12), 12, 1)
  B <- matrix(c(2, -1, 0, 1, 3, -2, 4, -1), 8, 1)
  X <- tcrossprod(A, B) + outer(1:12, 1:8, function(i, j) .1 * sin(i * j))
  dimnames(X) <- list(paste0("s", 1:12), paste0("g", 1:8))
  for (side in c("rows", "columns")) {
    Y <- if (side == "rows") X else t(X)
    priors <- list(ebnm::ebnm_point_exponential, ebnm::ebnm_point_laplace)
    if (side == "columns") priors <- rev(priors)
    set.seed(93)
    before <- .Random.seed
    fit <- fit_seminonnegative_ebmf(Y, side, max_factors = 2, seed = 4,
      backfit = TRUE, nullcheck = FALSE, verbose = 0, feature_scale = "intensity")
    expect_identical(.Random.seed, before)
    native <- .with_seed(4, flashier::flash(Y, ebnm_fn = priors, greedy_Kmax = 2,
      backfit = TRUE, nullcheck = FALSE, verbose = 0))
    expect_equal(fit$L_pm, native$L_pm, tolerance = 1e-8)
    expect_equal(fit$F_pm, native$F_pm, tolerance = 1e-8)
    expect_equal(fit$elbo, native$elbo, tolerance = 1e-8)
    expect_true(all(factor_activity(fit) >= 0))
    expect_true(any(factor_effects(fit) < 0))
    expect_identical(attr(fit, "flashbridge_metadata")$prior_support,
      c(activity = "nonnegative", effects = "signed"))
    expect_error(standardize_factors(fit, orientation = "anchor_feature"), "signed prior")
    expect_identical(attr(fit, "analysis_metadata")$engine, "flashier")
  }
})

test_that("semi-nonnegative priors reject invalid support and retain unknown custom support", {
  X <- fixture_matrix()
  expect_error(fit_seminonnegative_ebmf(X, "rows", activity_prior = ebnm::ebnm_normal), "nonnegative")
  expect_error(fit_seminonnegative_ebmf(X, "rows", effect_prior = ebnm::ebnm_point_exponential), "fit_nonnegative_ebmf")
  expect_error(fit_seminonnegative_ebmf(X, "rows", initialization = "NNLM"), "native_greedy")
  expect_error(fit_seminonnegative_ebmf(X, "rows", ebnm_fn = ebnm::ebnm_normal), "activity_prior")
  custom <- function(x, s, ...) ebnm::ebnm_point_laplace(x, s, ...)
  fit <- fit_seminonnegative_ebmf(X, "rows", effect_prior = custom,
    max_factors = 1, seed = 4, verbose = 0, feature_scale = "intensity")
  expect_identical(attr(fit, "flashbridge_metadata")$prior_support["effects"], c(effects = "unknown"))
  expect_identical(attr(fit, "flashbridge_metadata")$refit_policy, "callback_required")
  metadata <- data.frame(sample = rownames(X), donor = rownames(X))
  expect_error(factor_stability(X, fit, metadata, "donor", B = 1, seed = 7), "explicit fit_fn")
})
