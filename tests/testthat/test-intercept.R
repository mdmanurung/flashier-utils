test_that("feature intercepts are offsets while adaptive native maps remain correct", {
  X <- sweep(fixture_matrix(), 2, seq_len(8), "+")
  for (side in c("rows", "columns")) {
    Y <- if (side == "rows") X else t(X)
    priors <- list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)
    if (side == "columns") priors <- rev(priors)
    fit <- fit_ebmf(Y, side, ebnm_fn = priors, intercept = "feature", max_factors = 2,
      backfit = TRUE, nullcheck = FALSE, verbose = 0, seed = 12, feature_scale = "intensity")
    set.seed(12)
    native <- flashier::flash_init(Y)
    native <- flashier::flash_set_verbose(native, verbose = 0)
    native <- flashier::flash_add_intercept(native, rowwise = side == "rows")
    native <- flashier::flash_greedy(native, Kmax = 2, ebnm_fn = priors, verbose = 0)
    native <- flashier::flash_backfit(native, verbose = 0)
    expect_equal(fit$L_pm, native$L_pm, tolerance = 1e-10)
    expect_equal(fit$F_pm, native$F_pm, tolerance = 1e-10)
    expect_equal(fit$elbo, native$elbo, tolerance = 1e-10)
    expect_identical(attr(fit, "flashbridge_metadata")$native_factor_ids, c("Intercept", "F1", "F2"))
    for (d_location in c("activity", "effects")) for (scaling in c("raw", "ldf")) {
      view <- standardize_factors(fit, scaling = scaling, d_location = d_location)
      expect_identical(view$manifest$factor_ids, c("F1", "F2"))
      expect_identical(view$manifest$native_indices, 2:3)
      native_offset <- if (side == "rows") fit$F_pm[, 1] else fit$L_pm[, 1]
      expect_equal(view$offset, native_offset, tolerance = 1e-8)
      expected <- if (side == "rows") stats::fitted(fit) else t(stats::fitted(fit))
      expect_equal(unname(sweep(tcrossprod(view$activity, view$effects), 2, view$offset, "+")),
        unname(expected), tolerance = 1e-8)
      out <- factor_uncertainty(fit, "activity", representation = view, factors = "F2")
      native <- if (side == "rows") fit$L_psd[, 3] else fit$F_psd[, 3]
      expect_equal(out$posterior.sd, unname(native * abs(view$manifest$activity_multiplier[2])), tolerance = 1e-8)
      expect_equal(factor_pve(view)$pve, unname(fit$pve[2:3]), tolerance = 1e-8)
      beta <- c(F1 = 1, F2 = 2)
      projection <- backproject_contrast(view, beta, "coefficients", factor_basis = view$manifest)
      expect_identical(projection$factor_mapping$native_index, 2:3)
    }
    intervals <- factor_intervals(fit, "activity", representation = list(scaling = "raw"),
      factors = "F2", ids = "s2", nsamp = 20, seed = 39)
    set.seed(39)
    draws <- fit$sampler(20)
    values <- vapply(draws, function(d) d[[if (side == "rows") 1 else 2]][2, 3], numeric(1))
    expect_equal(intervals$mean, mean(values))
    expect_equal(intervals$lower, unname(quantile(values, .025)))
    reordered <- reorder_factors(fit, c("F2", "F1"))
    expect_equal(unname(stats::fitted(reordered)), unname(stats::fitted(fit)), tolerance = 1e-8)
    expect_equal(factor_effects(reordered)[, 1], factor_effects(fit)[, 2], ignore_attr = TRUE)
    fixed <- fix_factors(fit, "F1", "effects")
    fixed <- refit_factors(fixed, control = list(maxiter = 20, verbose = 0), seed = 12)
    expect_equal(factor_effects(fixed, representation = list(scaling = "raw"))[, 1],
      factor_effects(fit, representation = list(scaling = "raw"))[, 1], tolerance = 1e-8)
    removed <- remove_factors(fit, "F1")
    expect_identical(.factor_ids(removed), "F2")
    expect_identical(standardize_factors(removed)$manifest$native_indices, 2L)
    expect_identical(standardize_factors(removed)$manifest$native_factor_ids, c("Intercept", "F2"))
    expect_error(fix_factors(fit, "(Intercept)", "activity", action = "unfix"), "Unknown")
    expect_error(remove_factors(fit, "Intercept"), "Unknown")
    broken <- fit
    attr(broken, "flashbridge_metadata")$factor_roles[1] <- "factor"
    expect_error(standardize_factors(broken), "protected intercept role")
  }
})

test_that("column-oriented intercept-only fits preserve offsets and empty axes", {
  X <- fixture_matrix()
  for (Y in list(X, sweep(X, 2, colMeans(X), "-"))) {
    fit <- fit_ebmf(t(Y), "columns", max_factors = 0, intercept = "feature", verbose = 0)
    view <- standardize_factors(fit)
    expect_identical(dim(view$activity), c(12L, 0L))
    expect_identical(dim(view$effects), c(8L, 0L))
    expect_equal(view$offset, fit$L_pm[, 1])
    expect_equal(unname(t(stats::fitted(fit))), unname(matrix(view$offset, 12, 8, byrow = TRUE)))
    expect_equal(nrow(factor_intervals(fit, "effects", nsamp = 2)), 0L)
  }
})

test_that("intercept-only fits have semantic rank zero, including a zero offset", {
  X <- fixture_matrix()
  centered <- sweep(X, 2, colMeans(X), "-")
  for (Y in list(X, centered)) {
    fit <- fit_ebmf(Y, "rows", max_factors = 0, intercept = "feature", verbose = 0)
    view <- standardize_factors(fit)
    expect_identical(dim(view$activity), c(12L, 0L))
    expect_identical(view$manifest$native_indices, integer())
    expect_equal(nrow(factor_pve(fit)), 0L)
    expected <- matrix(view$offset, 12, 8, byrow = TRUE, dimnames = dimnames(X))
    expect_equal(unname(stats::fitted(fit)), unname(expected), tolerance = 1e-8)
    expect_equal(nrow(summarize_factors(fit)), 0L)
    expect_identical(view$manifest$offset_model, "feature")
    for (scaling in c("raw", "ldf")) {
      expect_identical(dim(standardize_factors(fit, scaling = scaling)$effects), c(8L, 0L))
    }
    uncertainty <- factor_uncertainty(fit, "activity")
    expect_equal(nrow(uncertainty), 0L)
    expect_type(uncertainty$status, "character")
    expect_type(uncertainty$reason, "character")
    intervals <- factor_intervals(fit, "activity", nsamp = 2, seed = 7)
    expect_equal(nrow(intervals), 0L)
    expect_type(intervals$mean, "double")
    expect_equal(factor_diagnostics(fit)$samples$programs, rep(0, 12))
    refit <- refit_factors(fit, control = list(verbose = 0))
    expect_equal(refit$n_factors, 1L)
  }
  expect_error(fit_ebmf(X, "rows", intercept = "feature",
    control = list(nullcheck = list(kset = 1))), "kset is protected")
})

test_that("adaptive nullchecking preserves a pre-zeroed protected factor", {
  fit <- fixture_fit()
  zeroed <- flashier::flash_factors_set_to_zero(fit, 1:2)
  result <- .nullcheck_adaptive(zeroed, 2L)
  expect_equal(result$fit$n_factors, 1L)
  expect_identical(result$retained, 1L)
  expect_error(.nullcheck_adaptive(zeroed, 2L, control = list(kset = 1L)), "protected")
  retained <- .nullcheck_adaptive(zeroed, 2L, control = list(remove = FALSE))
  expect_equal(retained$fit$n_factors, 2L)
})
