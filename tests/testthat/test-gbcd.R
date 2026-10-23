gbcd_preprocessing <- function() list(scope = "synthetic", description = "natural-log simulated expression", log_base = "e")
gbcd_fixture <- function() {
  .with_seed(713, {
    n <- 48L; p <- 30L
    A <- cbind(seq(.3, 1.3, length.out = n), rep(c(0, 1), each = n / 2), rep(c(0, 1), n / 2))
    B <- matrix(rnorm(p * 3, 1, .7), p, 3)
    X <- tcrossprod(A, B) + matrix(rnorm(n * p, sd = .15), n)
    dimnames(X) <- list(paste0("s", seq_len(n)), paste0("g", seq_len(p)))
    X
  })
}
gbcd_fake <- function(Y) {
  A <- cbind(Baseline = seq(.1, 1, length.out = nrow(Y)), GEP1 = rep(c(0, 1), length.out = nrow(Y)))
  rownames(A) <- rownames(Y)
  B <- matrix(seq_len(ncol(Y) * 2) / 10, ncol(Y), 2,
    dimnames = list(colnames(Y), colnames(A)))
  list(L = A, F = list(lfc = B, z_score = B, lfsr = B * 0 + .1))
}

test_that("GBCD preflight mirrors the pinned positive-entry branch and guards workspaces", {
  Y <- matrix(0, 8, 4)
  Y[1] <- -1
  expect_equal(.gbcd_preflight(Y, 2, 1e9)$strategy, "lowrank_plus_sparse")
  Y[] <- 1
  expect_equal(.gbcd_preflight(Y, 2, 1e9)$strategy, "dense_covariance")
  expect_equal(.gbcd_preflight(Matrix::Matrix(Y, sparse = TRUE), 2, 1e9),
    .gbcd_preflight(Y, 2, 1e9))
  expect_error(.gbcd_preflight(Y, 2, 1), "workspace")
  expect_error(.gbcd_preflight(matrix(0, 8, 1000), 100, 1e6), "workspace")
})

test_that("GBCD conversion keeps effects, labels and conditional statistics separate", {
  Y <- gbcd_fixture(); raw <- gbcd_fake(Y)
  raw$F$z_score[1, 1] <- Inf
  out <- .gbcd_output(raw, Y, "rows", "natural_log", gbcd_preprocessing(), list(), TRUE)
  view <- out$representation
  expect_equal(view$activity, raw$L)
  expect_equal(view$effects, log(2) * raw$F$lfc)
  expect_equal(tcrossprod(view$activity, view$effects), tcrossprod(raw$L, log(2) * raw$F$lfc))
  expect_true(all(view$offset == 0))
  expect_true(all(is.na(view$manifest$native_indices)))
  expect_identical(view$manifest$mapping_scope, "gbcd_components")
  expect_identical(view$manifest$constraints, c(activity = "unknown", effects = "signed"))
  expect_true(all(.factor_support(view$manifest)$activity_support == "unknown"))
  expect_true(is.na(out$signatures$z_score[1]))
  expect_equal(out$signatures$lfsr[1], .1)
  expect_false(out$signatures$statistics_available[1])
  expect_identical(out$native, raw)
  expect_true(all(c("gbcd", "ashr") %in% attr(provenance(result = view), "dependencies")$package))
  expect_equal(nrow(factor_activity(view)), nrow(Y))
  expect_equal(nrow(factor_effects(view)), ncol(Y))
  expect_s3_class(plot_factor_features(view, "GEP1"), "ggplot")
  expect_error(factor_pve(view), "unavailable")
  expect_error(factor_uncertainty(view), "native")
  expect_error(project_factor_activity(view, Y, "rows", "natural_log", gbcd_preprocessing()), "prior")
  contrast <- setNames(c(1, -1), colnames(raw$L))
  back <- backproject_contrast(view, contrast, "coefficients", factor_basis = view$manifest)
  expect_true(all(is.na(back$factor_mapping$native_index)))
  expect_equal(back$effects$estimate, as.vector(view$effects %*% contrast))
  broken <- raw; broken$L <- broken$L[, 1, drop = FALSE]
  expect_error(.gbcd_output(broken, Y, "rows", "natural_log", gbcd_preprocessing(), list(), FALSE), "two finite")
  broken <- raw; broken$F$lfc[1, 1] <- Inf
  expect_error(.gbcd_output(broken, Y, "rows", "natural_log", gbcd_preprocessing(), list(), FALSE), "finite")
  broken <- raw; colnames(broken$F$lfsr) <- rev(colnames(broken$L))
  expect_error(.gbcd_output(broken, Y, "rows", "natural_log", gbcd_preprocessing(), list(), FALSE), "axes")
  broken <- raw; broken$F$lfsr[1] <- 2
  expect_error(.gbcd_output(broken, Y, "rows", "natural_log", gbcd_preprocessing(), list(), FALSE), "lfsr")
  broken <- view; broken$manifest$native_indices <- 1:2
  expect_error(.validate_view(broken), "GBCD component mapping")
})

test_that("GBCD validates caller declarations before engine execution", {
  Y <- gbcd_fixture(); prep <- gbcd_preprocessing()
  call <- function(...) fit_gbcd_programs(Y, "rows", 3, "natural_log", prep, ...)
  expect_error(fit_gbcd_programs(Y, "rows", 1, "natural_log", prep), "kmax")
  expect_error(fit_gbcd_programs(Y, "rows", 3, "unknown", prep), "units")
  expect_error(fit_gbcd_programs(Y, "rows", 3, "natural_log", list()), "preprocessing")
  prep$log_base <- "2"
  expect_error(fit_gbcd_programs(Y, "rows", 3, "natural_log", prep), "log_base")
  prep <- gbcd_preprocessing()
  expect_error(call(prior = ebnm::ebnm_normal), "nonnegative")
  expect_error(call(Y = Y), "protected")
  expect_error(call(maxiter1 = 2^32), "maxiter1")
  expect_error(call(maxiter1 = 0), "maxiter1")
  expect_error(call(control = list(corr_thres = 2)), "corr_thres")
  expect_error(call(control = list(warmstart = NULL)), "control")
  expect_error(call(control = NULL), "control")
  expect_error(call(verbose = NULL), "verbose")
  expect_error(call(dense_budget_bytes = 1), "workspace")
  Y[1] <- NA_real_
  expect_error(fit_gbcd_programs(Y, "rows", 3, "natural_log", gbcd_preprocessing()), "Missing")
})

test_that("installed pinned GBCD agrees exactly with direct calls in both orientations", {
  skip_if_not_installed("gbcd")
  expect_identical(utils::packageDescription("gbcd")$RemoteSha, .gbcd_pin)
  Y <- gbcd_fixture()
  settings <- list(maxiter1 = 20L, maxiter2 = 20L, maxiter3 = 20L,
    control = list(corr_thres = 0), verbose = 0L)
  before <- serialize(Y, NULL)
  # Bounded parity checks compare identical iteration limits; convergence is separate.
  direct <- suppressWarnings(.with_seed(12, do.call(gbcd::fit_gbcd, c(list(Y = Y, Kmax = 3), settings))))
  for (side in c("rows", "columns")) {
    X <- if (side == "rows") Y else t(Y)
    wrapped <- suppressWarnings(.with_seed(12, do.call(fit_gbcd_programs, c(list(X = X, sample_side = side,
      kmax = 3, feature_scale = "natural_log", preprocessing = gbcd_preprocessing(), keep_native = TRUE), settings))))
    expect_identical(wrapped$native, direct)
    expect_equal(wrapped$representation$activity, direct$L)
    expect_equal(wrapped$representation$effects, log(2) * direct$F$lfc)
    expect_identical(wrapped$analysis_metadata$engine, "gbcd")
    deps <- attr(provenance(result = wrapped), "dependencies")
    expect_true(all(c("gbcd", "ashr") %in% deps$package))
    expect_identical(deps$runtime_source_commit[deps$package == "gbcd"], .gbcd_pin)
  }
  expect_identical(serialize(Y, NULL), before)
})


test_that("GBCD refuses an unqualified installed engine before fitting", {
  original <- utils::packageDescription
  testthat::local_mocked_bindings(packageDescription = function(pkg, ...) {
    if (identical(pkg, "gbcd")) return(list(Version = "0.2-10", RemoteSha = "wrong"))
    original(pkg, ...)
  }, .package = "utils")
  expected <- if (requireNamespace("gbcd", quietly = TRUE)) "audited source pin" else "Install optional"
  expect_error(fit_gbcd_programs(gbcd_fixture(), "rows", 3, "natural_log", gbcd_preprocessing()), expected)
})

test_that("missing GBCD fails before engine execution", {
  original <- base::requireNamespace
  testthat::local_mocked_bindings(requireNamespace = function(package, ...) {
    if (identical(package, "gbcd")) return(FALSE)
    original(package, ...)
  }, .package = "base")
  expect_error(fit_gbcd_programs(gbcd_fixture(), "rows", 3, "natural_log",
    gbcd_preprocessing()), "Install optional gbcd")
})

test_that("installed GBCD sparse low-rank branch agrees with direct engine output", {
  skip_if_not_installed("gbcd")
  Y <- .with_seed(714, {
    X <- matrix(0, 60, 30, dimnames = list(paste0("s", 1:60), paste0("g", 1:30)))
    for (k in 1:3) X[((k-1)*20+1):(k*20), ((k-1)*10+1):(k*10)] <- matrix(rexp(200, .5), 20)
    Matrix::Matrix(X, sparse = TRUE)
  })
  settings <- list(maxiter1 = 20L, maxiter2 = 20L, maxiter3 = 20L,
    control = list(corr_thres = 0), verbose = 0L)
  direct <- suppressWarnings(.with_seed(14, do.call(gbcd::fit_gbcd,
    c(list(Y = Y, Kmax = 3), settings))))
  wrapped <- suppressWarnings(.with_seed(14, do.call(fit_gbcd_programs,
    c(list(X = Y, sample_side = "rows", kmax = 3, feature_scale = "natural_log",
      preprocessing = gbcd_preprocessing(), keep_native = TRUE), settings))))
  expect_identical(wrapped$diagnostics$strategy, "lowrank_plus_sparse")
  expect_identical(wrapped$native, direct)
})
