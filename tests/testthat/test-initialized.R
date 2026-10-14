initialized_fixture <- function() {
  X <- fixture_matrix()
  A <- cbind(F1 = seq(-1, 1, length.out = 12), anchor = sin(1:12))
  B <- matrix(c(2,-1,0,1,3,-2,4,0,1,2,-3,0,-2,1,0,4), 8, 2,
    dimnames = list(colnames(X), colnames(A)))
  rownames(A) <- rownames(X)
  list(X = X, A = A, B = B)
}

test_that("initialized fitting matches public calls with exact IDs in both orientations", {
  x <- initialized_fixture()
  controls <- list(backfit = list(maxiter = 50, extrapolate = FALSE))
  for (side in c("rows", "columns")) for (fixed in c("none", "activity", "effects")) {
    Y <- if (side == "rows") x$X else t(x$X)
    set.seed(31)
    previous <- .Random.seed
    fit <- suppressWarnings(fit_initialized_ebmf(Y, side, x$A[12:1, ], x$B[8:1, ],
      activity_prior = ebnm::ebnm_normal, fixed_side = fixed, seed = 17,
      control = controls, verbose = 0))
    expect_identical(.Random.seed, previous)
    priors <- list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)
    pair <- list(x$A, x$B)
    if (side == "columns") {
      pair <- rev(pair)
      priors <- rev(priors)
    }
    set.seed(17)
    native <- flashier::flash_init(Y)
    native <- flashier::flash_set_verbose(native, 0)
    native <- flashier::flash_factors_init(native, pair, ebnm_fn = priors)
    if (fixed != "none") native <- flashier::flash_factors_fix(native, 1:2,
      which_dim = if ((fixed == "activity") == (side == "rows")) "loadings" else "factors")
    native <- suppressWarnings(flashier::flash_backfit(native, maxiter = 50,
      extrapolate = FALSE, verbose = 0))
    expect_equal(fit$L_pm, native$L_pm, tolerance = 1e-10)
    expect_equal(fit$F_pm, native$F_pm, tolerance = 1e-10)
    expect_equal(fit$elbo, native$elbo, tolerance = 1e-10)
    expect_identical(.factor_ids(fit), c("F1", "anchor"))
    view <- standardize_factors(fit, scaling = "raw")
    expect_identical(colnames(view$activity), colnames(x$A))
    expect_identical(unname(view$manifest$factor_roles), rep("initialized", 2))
    if (fixed != "none") {
      expect_equal(view[[fixed]], if (fixed == "activity") x$A else x$B)
      other <- if (fixed == "activity") "effects" else "activity"
      expect_error(fix_factors(fit, "F1", other), "single mode")
    }
    metadata <- attr(fit, "flashbridge_metadata")
    expect_identical(metadata$refit_policy, "callback_required")
    expect_identical(metadata$initialization$factor_ids, c("F1", "anchor"))
    expect_false(any(c("activity_init", "effect_init") %in% names(attr(fit, "analysis_metadata")$parameters)))
    md <- data.frame(sample = rownames(x$X), donor = rep(1:4, each = 3))
    expect_error(factor_stability(Y, fit, md, "donor", B = 1, seed = 9), "explicit fit_fn")
  }
})

test_that("new factors preserve initialized IDs through selected-only nullchecking", {
  x <- initialized_fixture()
  fit <- suppressWarnings(fit_initialized_ebmf(x$X, "rows", x$A, x$B,
    activity_prior = ebnm::ebnm_normal, fixed_side = "activity", add_factors = 2,
    nullcheck_new = FALSE, seed = 17, verbose = 0,
    control = list(backfit = list(maxiter = 50), greedy = list(maxiter = 50))))
  ids <- .factor_ids(fit)
  expect_gt(length(ids), 2L)
  expect_identical(ids[1:2], c("F1", "anchor"))
  expect_identical(ids[-(1:2)], paste0("F", seq_len(length(ids) - 2) + 1L))
  zeroed <- flashier::flash_factors_set_to_zero(fit, seq_len(fit$n_factors))
  checked <- .nullcheck_adaptive(zeroed, seq_len(fit$n_factors)[-(1:2)])
  expect_identical(checked$retained, 1:2)
  expect_equal(checked$fit$n_factors, 2L)
  updated <- .native_changed(fit, checked$fit, ids[checked$retained], "nullcheck_new")
  expect_identical(.factor_ids(updated), c("F1", "anchor"))
  expect_identical(unname(standardize_factors(updated)$manifest$factor_roles), rep("initialized", 2))
  reordered <- reorder_factors(fit, rev(ids))
  expect_identical(.factor_ids(reordered), rev(ids))
  expect_equal(unname(stats::fitted(reordered)), unname(stats::fitted(fit)))
})

test_that("initialized fitting validates axes, support and protected controls", {
  x <- initialized_fixture()
  run <- function(A = x$A, B = x$B, ...) fit_initialized_ebmf(x$X, "rows", A, B, verbose = 0, ...)
  bad <- x$B[, 2:1]
  expect_error(run(B = bad), "names and order")
  bad <- x$A
  rownames(bad)[1] <- "absent"
  expect_error(run(A = bad), "exactly match")
  bad <- x$A
  bad[1, 1] <- NA
  expect_error(run(A = bad), "Missing")
  bad[1, 1] <- Inf
  expect_error(run(A = bad), "Inf")
  bad <- x$A
  bad[, 1] <- 0
  expect_error(run(A = bad), "nonzero")
  expect_error(run(activity_prior = ebnm::ebnm_point_exponential), "nonnegative")
  expect_error(run(add_factors = -1), "add_factors")
  expect_error(run(fixed_side = "both"), "arg")
  expect_error(run(control = list(backfit = list(kset = 2))), "protected")
  expect_error(run(control = list(nullcheck = list(kset = 1))), "protected")
  expect_error(run(control = list(greedy = list(Kmax = 5))), "protected")
  expect_error(run(control = list(backfit = 1)), "Pipeline control")
  set.seed(42)
  previous <- .Random.seed
  failing_prior <- function(x, s, g_init = NULL, fix_g = FALSE, output = NULL) {
    if (!identical(unname(x), c(-10, 0, 10))) stop("deliberate_native_failure")
    ebnm::ebnm_normal(x, s, g_init = g_init, fix_g = fix_g, output = output)
  }
  failed <- tryCatch(run(seed = 17, activity_prior = failing_prior), error = identity)
  expect_s3_class(failed, "error")
  expect_match(conditionMessage(failed), "deliberate_native_failure")
  expect_identical(.Random.seed, previous)
  colnames(x$A)[1] <- colnames(x$B)[1] <- "Intercept"
  expect_error(run(), "reserved")
})

test_that("added-factor pipeline matches independent selected-only native removal", {
  x <- initialized_fixture()
  control <- list(backfit = list(maxiter = 50), greedy = list(maxiter = 50))
  for (side in c("rows", "columns")) {
    Y <- if (side == "rows") x$X else t(x$X)
    fit <- suppressWarnings(fit_initialized_ebmf(Y, side, x$A, x$B,
      activity_prior = ebnm::ebnm_normal, fixed_side = "activity", add_factors = 2,
      seed = 17, control = control, verbose = 0))
    priors <- list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)
    pair <- list(x$A, x$B)
    if (side == "columns") {
      priors <- rev(priors)
      pair <- rev(pair)
    }
    set.seed(17)
    native <- flashier::flash_set_verbose(flashier::flash_init(Y), 0)
    native <- flashier::flash_factors_init(native, pair, ebnm_fn = priors)
    native <- flashier::flash_factors_fix(native, 1:2,
      which_dim = if (side == "rows") "loadings" else "factors")
    native <- suppressWarnings(flashier::flash_backfit(native, maxiter = 50, verbose = 0))
    native <- suppressWarnings(flashier::flash_greedy(native, Kmax = 2, ebnm_fn = priors,
      maxiter = 50, verbose = 0))
    native <- suppressWarnings(flashier::flash_backfit(native, maxiter = 50, verbose = 0))
    new_indices <- setdiff(seq_len(native$n_factors), 1:2)
    expected_ids <- c("F1", "anchor", paste0("F", seq_along(new_indices) + 1L))
    native <- flashier::flash_nullcheck(native, kset = new_indices, remove = FALSE, verbose = 0)
    moments <- lapply(1:2, function(n) flashier::flash_fit_get_p2m(flashier::flash_fit(native), n))
    for (k in rev(new_indices)) if (all(moments[[1]][, k] == 0) && all(moments[[2]][, k] == 0)) {
      native <- flashier::flash_factors_remove(native, k)
      expected_ids <- expected_ids[-k]
    }
    expect_equal(fit$L_pm, native$L_pm, tolerance = 1e-10)
    expect_equal(fit$F_pm, native$F_pm, tolerance = 1e-10)
    expect_equal(fit$elbo, native$elbo, tolerance = 1e-10)
    expect_identical(.factor_ids(fit), expected_ids)
    expect_equal(standardize_factors(fit, scaling = "raw")$activity[, 1:2], x$A)
  }
})
