# FU-0: backend capability preflight — fail fast on missing native symbols
# instead of failing late inside an expensive fit/nullcheck path.

test_that("flashier preflight passes on the pinned runtime", {
  expect_true(flashier.utils:::.flashier_preflight())
  expect_true(flashier.utils:::.flashier_preflight(need_nullcheck_moments = TRUE))
})

test_that("preflight produces a targeted diagnostic for missing native capability", {
  ## simulate an incompatible flashier: one late-needed symbol absent
  local_mocked_bindings(
    .ns_has_function = function(pkg, sym) !identical(sym, "flash_fit_get_p2m"),
    .package = "flashier.utils"
  )
  ## core-only check still passes (the missing symbol is late-path only)
  expect_true(flashier.utils:::.flashier_preflight())
  ## nullcheck path fails fast with an actionable message
  expect_error(
    flashier.utils:::.flashier_preflight(need_nullcheck_moments = TRUE),
    "flash_fit_get_p2m"
  )
  expect_error(
    flashier.utils:::.flashier_preflight(need_nullcheck_moments = TRUE),
    "fd44811"  # message names the pinned upstream revision
  )
})

test_that("fit_ebmf fails fast (before fitting) when a late capability is absent", {
  local_mocked_bindings(
    .ns_has_function = function(pkg, sym) !identical(sym, "flash_fit_get_p2m"),
    .package = "flashier.utils"
  )
  X <- matrix(rpois(40, 5), 8, 5, dimnames = list(paste0("s", 1:8), paste0("g", 1:5)))
  t0 <- Sys.time()
  expect_error(
    fit_ebmf(X, sample_side = "rows", max_factors = 2, nullcheck = TRUE, verbose = 0),
    "lacks native capabilities"
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  expect_lt(elapsed, 5)  # failed at entry, not after an expensive fit
  ## disabling nullcheck avoids the late-path requirement and still fits
  expect_silent(
    fit_ebmf(X, sample_side = "rows", max_factors = 2, nullcheck = FALSE,
             backfit = FALSE, verbose = 0, seed = 1)
  )
})

test_that("pinned runtime fits fresh-native zero/affine/fixed constraints (FU-0 oracle)", {
  skip_if_not_installed("flashier")
  X <- matrix(rpois(60, 3), 10, 6, dimnames = list(paste0("s", 1:10), paste0("g", 1:6)))
  ## zero-factor request is a typed empty result, not an error
  f0 <- fit_ebmf(X, sample_side = "rows", max_factors = 0, nullcheck = FALSE,
                 backfit = FALSE, verbose = 0, seed = 1)
  expect_true(inherits(f0, "flash") || inherits(f0, "list"))
  ## affine intercept path runs on the pinned engine
  f1 <- fit_ebmf(X, sample_side = "rows", max_factors = 2, intercept = "feature",
                 nullcheck = FALSE, backfit = FALSE, verbose = 0, seed = 1)
  expect_true(!is.null(f1))
})
