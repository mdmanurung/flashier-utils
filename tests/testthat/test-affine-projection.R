test_that("new measurement offsets agree with public native affine fits", {
  X <- fixture_matrix()
  source <- standardize_factors(fixture_fit(), scaling = "raw")
  source <- .make_representation(source$activity, source$effects,
    utils::modifyList(source$manifest, list(offset_model = "feature")))
  A <- source$activity[12:1, , drop = FALSE]
  Y <- sweep(X[12:1, 1:3], 2, c(-5, 2, 8), "+")
  design <- cbind(1, A)
  b0 <- t(qr.coef(qr(design), Y))
  for (side in c("rows", "columns")) {
    priors <- list(ebnm::ebnm_normal, ebnm::ebnm_point_exponential)
    if (side == "columns") priors <- rev(priors)
    result <- suppressWarnings(project_factor_features(source,
      if (side == "rows") Y else t(Y), side, "protein", ebnm_fn = priors))
    offset_priors <- priors
    offset_priors[[if (side == "rows") 2 else 1]] <- ebnm::ebnm_point_normal
    native <- flashier::flash_init(if (side == "rows") Y else t(Y),
      var_type = if (side == "rows") 2 else 1)
    pair <- list(design[, 1, drop = FALSE], b0[, 1, drop = FALSE])
    native <- flashier::flash_factors_init(native,
      if (side == "rows") pair else rev(pair), ebnm_fn = offset_priors)
    pair <- list(A, b0[, -1, drop = FALSE])
    native <- flashier::flash_factors_init(native,
      if (side == "rows") pair else rev(pair), ebnm_fn = priors)
    native <- flashier::flash_factors_fix(native, 1:3,
      which_dim = if (side == "rows") "loadings" else "factors")
    native <- suppressWarnings(flashier::flash_backfit(native,
      maxiter = 200, extrapolate = FALSE, verbose = 0))
    expected <- if (side == "rows") native$F_pm else native$L_pm
    view <- result$representation
    expect_identical(view$activity, A)
    expect_equal(view$offset, expected[, 1], tolerance = 1e-8)
    expect_equal(unname(view$effects), unname(expected[, -1]), tolerance = 1e-8)
    expect_true(any(view$offset < 0))
    expect_true(all(view$effects >= 0))
    expect_identical(view$manifest$native_indices, 2:3)
    restored <- standardize_factors(result$fit, scaling = "raw")
    expect_equal(restored$activity, view$activity)
    expect_equal(restored$effects, view$effects)
    expect_equal(restored$offset, view$offset)
    expect_equal(factor_pve(view)$pve, unname(native$pve[2:3]))
    out <- factor_uncertainty(result$fit, "effects", representation = view)
    psd <- if (side == "rows") native$F_psd else native$L_psd
    expect_equal(out$posterior.sd, as.vector(psd[, -1]))
  }
  expect_error(project_factor_features(source, Y, "rows", "protein",
    ebnm_fn = ebnm::ebnm_normal, dense_budget_bytes = 1), "allocations")
  constant <- source$activity
  constant[, 1] <- 1
  bad <- .make_representation(constant, source$effects, source$manifest)
  expect_error(project_factor_features(bad, Y, "rows", "protein",
    ebnm_fn = ebnm::ebnm_normal), "Rank-deficient")
})

test_that("intercept-only sources learn offsets for new measurements", {
  X <- fixture_matrix()
  fit <- fit_ebmf(sweep(X, 2, colMeans(X), "-"), "rows",
    max_factors = 0, intercept = "feature", verbose = 0)
  Y <- sweep(X[, 1:3], 2, c(-5, 2, 8), "+")
  for (side in c("rows", "columns")) {
    result <- project_factor_features(fit,
      if (side == "rows") Y else t(Y), side, "protein")
    view <- result$representation
    expect_identical(dim(view$activity), c(12L, 0L))
    expect_identical(dim(view$effects), c(3L, 0L))
    expect_true(any(view$offset != 0))
    expect_identical(view$manifest$native_indices, integer())
    expect_equal(standardize_factors(result$fit, scaling = "raw")$offset, view$offset)
    expected <- if (side == "rows") stats::fitted(result$fit) else t(stats::fitted(result$fit))
    expect_equal(unname(matrix(view$offset, 12, 3, byrow = TRUE)), unname(expected))
  }
})

test_that("activity transfer retains external offsets without a native intercept", {
  X <- sweep(fixture_matrix(), 2, seq_len(8), "+")
  fit <- fit_ebmf(X, "rows", intercept = "feature", max_factors = 2,
    backfit = TRUE, nullcheck = FALSE, seed = 12, verbose = 0, feature_scale = "intensity")
  view <- standardize_factors(fit, scaling = "raw")
  result <- suppressWarnings(project_factor_activity(view, X, "rows", "intensity",
    list(scope = "identity")))
  expect_null(result$representation$manifest$intercept_index)
  expect_equal(result$representation$manifest$source_coordinates$intercept_index, 1L)
  expect_equal(result$fit$n_factors, 2L)
  expect_identical(result$representation$manifest$native_indices, 1:2)
  restored <- standardize_factors(result$fit, scaling = "raw")
  expect_equal(restored$offset, view$offset)
  expect_equal(restored$activity, result$representation$activity)
})
