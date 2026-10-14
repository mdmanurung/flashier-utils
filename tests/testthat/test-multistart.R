test_that("multistart candidates reproduce direct native seeded fits", {
  X <- fixture_matrix()
  for (side in c("rows", "columns")) {
    Y <- if (side == "rows") X else t(X)
    priors <- list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)
    if (side == "columns") priors <- rev(priors)
    set.seed(41)
    rng <- .Random.seed
    result <- fit_multistart_ebmf(Y, side, ebnm_fn = priors, n_starts = 2,
      init_seeds = c(31, 72), max_factors = 2, backfit = TRUE,
      nullcheck = FALSE, verbose = 0, feature_scale = "intensity")
    expect_identical(.Random.seed, rng)
    expect_identical(result$candidates$status, c("supported", "supported"))
    expect_true(all(result$candidates$initializer_calls > 0L))
    for (i in 1:2) {
      native <- .with_seed(result$candidates$init_seed[i], {
        init <- function(f) flashier::flash_greedy_init_default(f,
          sign_constraints = c(0, 0), seed = result$candidates$init_seed[i])
        fit <- flashier::flash_init(Y)
        fit <- flashier::flash_greedy(fit, Kmax = 2, ebnm_fn = priors,
          init_fn = init, verbose = 0)
        flashier::flash_backfit(fit, verbose = 0)
      })
      expect_equal(result$fits[[i]]$L_pm, native$L_pm, tolerance = 1e-10)
      expect_equal(result$fits[[i]]$F_pm, native$F_pm, tolerance = 1e-10)
      expect_equal(result$candidates$elbo[i], flashier::flash_fit_get_elbo(flashier::flash_fit(native)), tolerance = 1e-10)
      expect_identical(attr(result$fits[[i]], "flashbridge_metadata")$refit_policy, "callback_required")
      expect_identical(attr(result$fits[[i]], "analysis_metadata")[["function"]], "fit_multistart_ebmf")
    }
    expect_identical(result$selected_start, which.max(result$candidates$elbo))
    metadata <- data.frame(sample = rownames(X), donor = rownames(X))
    expect_error(factor_stability(Y, result$best_fit, metadata, "donor", B = 1, seed = 3), "explicit fit_fn")
    old <- result$best_fit
    attr(old, "flashbridge_metadata")$refit_policy <- NULL
    expect_error(factor_stability(Y, old, metadata, "donor", B = 1, seed = 3), "explicit fit_fn")
  }
})

test_that("actual native random initializers receive seeds and sign constraints", {
  original <- flashier::flash_greedy_init_default
  calls <- list()
  testthat::local_mocked_bindings(flash_greedy_init_default = function(flash,
      sign_constraints = NULL, seed = 666, ...) {
    values <- original(flash, sign_constraints = sign_constraints, seed = seed, ...)
    calls[[length(calls) + 1L]] <<- list(seed = seed, constraints = sign_constraints, values = values)
    values
  }, .package = "flashier")
  X <- fixture_matrix()
  for (side in c("rows", "columns")) {
    priors <- list(ebnm::ebnm_point_exponential, ebnm::ebnm_point_laplace)
    if (side == "columns") priors <- rev(priors)
    calls <- list()
    result <- fit_multistart_ebmf(if (side == "rows") X else t(X), side,
      ebnm_fn = priors, n_starts = 2, init_seeds = c(11, 12),
      max_factors = 1, nullcheck = FALSE, verbose = 0)
    expect_setequal(vapply(calls, `[[`, integer(1), "seed"), c(11L, 12L))
    constrained <- if (side == "rows") 1L else 2L
    for (call in calls) {
      expect_equal(call$constraints, if (side == "rows") c(1, 0) else c(0, 1))
      expect_true(all(call$values[[constrained]] >= 0))
    }
    expect_true(all(factor_activity(result$best_fit, representation = list(scaling = "raw")) >= 0))
  }
})

test_that("candidate failures and all-fail reasons survive scoped RNG", {
  X <- fixture_matrix()
  original <- flashier::flash_greedy_init_default
  testthat::local_mocked_bindings(flash_greedy_init_default = function(flash,
      sign_constraints = NULL, seed = 666, ...) {
    if (seed == 11L) stop("deliberate initializer failure")
    original(flash, sign_constraints = sign_constraints, seed = seed, ...)
  }, .package = "flashier")
  set.seed(3)
  rng <- .Random.seed
  result <- fit_multistart_ebmf(X, "rows", n_starts = 2, init_seeds = c(11, 12),
    max_factors = 1, nullcheck = FALSE, verbose = 0)
  expect_identical(.Random.seed, rng)
  expect_identical(result$candidates$status, c("failed", "supported"))
  expect_true(nzchar(result$candidates$reason[1]))
  expect_identical(result$selected_start, 2L)
  expect_null(result$fits[[1]])
  failed <- tryCatch(fit_multistart_ebmf(X, "rows", n_starts = 2,
    init_seeds = c(11, 11), max_factors = 1, verbose = 0), error = identity)
  expect_match(conditionMessage(failed), "unique")
})

test_that("ties, finite ELBO eligibility and no-selection are explicit", {
  X <- fixture_matrix()
  fixed <- fixture_fit()
  testthat::local_mocked_bindings(fit_ebmf = function(...) fixed, .package = "flashier.utils")
  result <- fit_multistart_ebmf(X, "rows", n_starts = 2, init_seeds = c(12, 11),
    max_factors = 1, verbose = 0, keep_fits = FALSE)
  expect_equal(result$candidates$elbo[1], result$candidates$elbo[2])
  expect_identical(result$selected_start, 1L)
  expect_null(result$fits)
  expect_s3_class(result$best_fit, "flash")
  expect_identical(attr(fixed, "flashbridge_metadata")$refit_policy, "ordinary")
  unselected <- fit_multistart_ebmf(X, "rows", n_starts = 2,
    init_seeds = c(11, 12), max_factors = 1, verbose = 0, select = "none")
  expect_null(unselected$best_fit)
  expect_true(is.na(unselected$selected_start))
  expect_false(any(unselected$candidates$selected))
  testthat::local_mocked_bindings(flash_fit_get_elbo = function(f) Inf, .package = "flashier")
  failed <- tryCatch(fit_multistart_ebmf(X, "rows", n_starts = 2,
    init_seeds = c(11, 12), max_factors = 1, verbose = 0), error = identity)
  expect_s3_class(failed, "error")
  expect_match(conditionMessage(failed), "All multistart candidates failed")
  expect_identical(failed$candidates$status, c("failed", "failed"))
  expect_true(all(grepl("ELBO is not finite", failed$candidates$reason)))
})

test_that("all initializer failures retain reasons and restore caller RNG", {
  X <- fixture_matrix()
  testthat::local_mocked_bindings(flash_greedy_init_default = function(...) {
    stop("deliberate initializer failure")
  }, .package = "flashier")
  set.seed(19)
  rng <- .Random.seed
  failed <- tryCatch(fit_multistart_ebmf(X, "rows", n_starts = 2,
    init_seeds = c(11, 12), max_factors = 1, verbose = 0), error = identity)
  expect_identical(.Random.seed, rng)
  expect_match(conditionMessage(failed), "All multistart candidates failed")
  expect_identical(failed$candidates$init_seed, c(11L, 12L))
  expect_identical(failed$candidates$status, c("failed", "failed"))
  expect_true(all(nzchar(failed$candidates$reason)))
})

test_that("fully nonnegative starts keep both native dimensions constrained", {
  result <- fit_multistart_ebmf(abs(fixture_matrix()), "rows", n_starts = 2,
    init_seeds = c(11, 12), max_factors = 1, nullcheck = FALSE,
    ebnm_fn = ebnm::ebnm_point_exponential, verbose = 0)
  expect_equal(result$analysis_metadata$native_sign_constraints, c(1, 1))
  expect_true(all(factor_activity(result$best_fit, representation = list(scaling = "raw")) >= 0))
  expect_true(all(factor_effects(result$best_fit, representation = list(scaling = "raw")) >= 0))
})

test_that("master-seeded numerical summaries are reproducible and compact", {
  X <- fixture_matrix()
  run <- function(...) fit_multistart_ebmf(X, "rows", n_starts = 2,
    seed = 3, max_factors = 1, nullcheck = FALSE, verbose = 0, ...)
  first <- run(keep_fits = FALSE)
  second <- run(keep_fits = FALSE)
  fields <- c("start", "init_seed", "status", "reason", "rank", "elbo", "candidate_digest", "initialization_digest", "initializer_calls", "selected")
  expect_identical(first$candidates[fields], second$candidates[fields])
  expect_false(anyDuplicated(first$candidates$init_seed) > 0)
  expect_null(first$fits)
  anonymous <- X
  dimnames(anonymous) <- NULL
  generated <- fit_multistart_ebmf(anonymous, "rows", n_starts = 2,
    seed = 3, max_factors = 1, nullcheck = FALSE, id_policy = "generate", verbose = 0)
  expect_identical(attr(generated$best_fit, "flashbridge_metadata")$generated_ids,
    c(rows = TRUE, columns = TRUE))
})

test_that("initializer callbacks do not serialize the caller's data and candidates", {
  make <- function(method, n) {
    retained_input <- numeric(n)
    retained_fits <- list(retained_input)
    .multistart_initializer(method, 11, c(0, 0))$callback
  }
  for (method in c("native_random", "softImpute", "irlba")) {
    expect_equal(length(serialize(make(method, 100000), NULL)),
      length(serialize(make(method, 1), NULL)))
  }
  result <- fit_multistart_ebmf(fixture_matrix(), "rows", n_starts = 2,
    init_seeds = c(11, 12), max_factors = 1, nullcheck = FALSE,
    verbose = 0, keep_fits = FALSE)
  path <- tempfile()
  on.exit(unlink(path))
  saveRDS(result$best_fit, path)
  restored <- readRDS(path)
  expect_identical(standardize_factors(restored)$manifest,
    standardize_factors(result$best_fit)$manifest)
  expect_equal(restored$L_pm, result$best_fit$L_pm)
  expect_equal(restored$F_pm, result$best_fit$F_pm)
})

test_that("signed spectral initializers use their public native helpers", {
  X <- fixture_matrix()
  for (method in c("softImpute", "irlba")) {
    skip_if_not_installed(method)
    result <- fit_multistart_ebmf(X, "rows", n_starts = 2, init_seeds = c(11, 12),
      init_method = method, max_factors = 1, nullcheck = FALSE,
      ebnm_fn = ebnm::ebnm_normal, verbose = 0)
    expect_identical(result$candidates$status, c("supported", "supported"))
    expect_identical(result$analysis_metadata$initializer_engine, method)
    expect_identical(result$analysis_metadata$initializer_version, as.character(packageVersion(method)))
    helper <- getExportedValue("flashier", paste0("flash_greedy_init_", method))
    init <- function(f) helper(f, seed = 11)
    native <- .with_seed(11, flashier::flash_greedy(flashier::flash_init(X),
      Kmax = 1, ebnm_fn = ebnm::ebnm_normal, init_fn = init, verbose = 0))
    expect_equal(result$fits[[1]]$L_pm, native$L_pm, tolerance = 1e-10)
    expect_equal(result$fits[[1]]$F_pm, native$F_pm, tolerance = 1e-10)
  }
})

test_that("multistart refuses ambiguous seeds and incompatible initialization", {
  X <- fixture_matrix()
  expect_error(fit_multistart_ebmf(X), "sample_side")
  expect_error(fit_multistart_ebmf(X, "rows", n_starts = 1, seed = 3), "n_starts")
  expect_error(fit_multistart_ebmf(X, "rows", max_factors = 0, seed = 3), "max_factors")
  expect_error(fit_multistart_ebmf(X, "rows"), "Supply seed")
  for (seeds in list(c(1, 1), c(1, NA), c(1, 2.5), c(1, -2), 1)) {
    expect_error(fit_multistart_ebmf(X, "rows", n_starts = 2, init_seeds = seeds), "init_seeds")
  }
  expect_error(fit_multistart_ebmf(X, "rows", seed = 3,
    ebnm_fn = function(...) ebnm::ebnm_normal(...)), "recognized")
  expect_error(fit_multistart_ebmf(X, "rows", seed = 3,
    control = list(greedy = list(init_fn = function(...) NULL))), "owns greedy")
  for (method in c("softImpute", "irlba")) {
    expect_error(fit_multistart_ebmf(X, "rows", seed = 3, init_method = method,
      ebnm_fn = ebnm::ebnm_point_exponential), "signed priors")
  }
  X[1, 1] <- NA
  expect_error(fit_multistart_ebmf(X, "rows", seed = 3, init_method = "irlba"), "dense complete")
  expect_error(fit_multistart_ebmf(Matrix::Matrix(fixture_matrix(), sparse = TRUE),
    "rows", seed = 3, init_method = "irlba"), "dense complete")
})
