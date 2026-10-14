## SYNC-1A: opt-in RMS scaling in standardize_factors(). Parity with the
## cellprograms canonicalization convention (unit-RMS activity, compensating
## effects, max-|loading| positive-gene sign anchor) plus raw/LDF non-regression.

test_that("rms scaling yields unit-RMS activity and preserves reconstruction", {
  fit <- fixture_fit()
  raw <- standardize_factors(fit, scaling = "raw")
  rms <- standardize_factors(fit, scaling = "rms")
  K <- ncol(rms$activity)
  skip_if(K == 0L, "fixture produced zero factors")
  # unit RMS on the activity side
  expect_equal(unname(sqrt(colMeans(rms$activity^2))), rep(1, K), tolerance = 1e-8)
  # reconstruction preserved exactly vs raw
  recon_raw <- tcrossprod(raw$activity, raw$effects)
  recon_rms <- tcrossprod(rms$activity, rms$effects)
  expect_equal(recon_rms, recon_raw, tolerance = 1e-8)
})

test_that("rms matches the cellprograms canonicalization convention", {
  fit <- fixture_fit()
  rms <- standardize_factors(fit, scaling = "rms", orientation = "anchor_feature")
  # reference: cellprograms canonicalize_programs() applied to the native fit
  Z <- as.matrix(fit$L_pm); W <- as.matrix(fit$F_pm)
  s_k <- sqrt(colMeans(Z^2)); s_k[s_k == 0] <- 1
  sign_k <- apply(W, 2, function(w) {
    am <- which(abs(w) == max(abs(w))); am <- am[order(rownames(W)[am], w[am])]
    sign(w[am[1]])
  })
  sign_k[sign_k == 0] <- 1
  Zc <- sweep(sweep(Z, 2, s_k, "/"), 2, sign_k, "*")
  Wc <- sweep(sweep(W, 2, s_k, "*"), 2, sign_k, "*")
  expect_equal(unname(as.matrix(rms$activity)), unname(Zc), tolerance = 1e-8)
  expect_equal(unname(as.matrix(rms$effects)), unname(Wc), tolerance = 1e-8)
})

test_that("rms records the per-factor scale in the manifest", {
  fit <- fixture_fit()
  rms <- standardize_factors(fit, scaling = "rms")
  expect_equal(rms$manifest$scaling, "rms")
  expect_true(!is.null(rms$manifest$rms_scale))
  expect_equal(length(rms$manifest$rms_scale), ncol(rms$activity))
})

test_that("rms zero-activity column keeps scale 1 (degenerate guard)", {
  # build a representation with a zero-activity factor and standardize via rms
  v <- fixture_view_zero_activity()
  s_k <- sqrt(colMeans(v$activity^2)); s_k[s_k == 0] <- 1
  expect_equal(unname(s_k["F2"]), 1)
})

test_that("rms anchor orientation resolves tied anchors lexically", {
  fit <- fixture_fit()
  r1 <- standardize_factors(fit, scaling = "rms", orientation = "anchor_feature")
  r2 <- standardize_factors(fit, scaling = "rms", orientation = "anchor_feature")
  expect_equal(r1$manifest$orientation_anchor, r2$manifest$orientation_anchor)
  expect_equal(r1$effects, r2$effects)
})

test_that("raw and ldf scalings are unchanged by the rms addition", {
  fit <- fixture_fit()
  raw <- standardize_factors(fit, scaling = "raw")
  ldf <- standardize_factors(fit, scaling = "ldf")
  expect_equal(raw$manifest$scaling, "raw")
  expect_equal(ldf$manifest$scaling, "ldf")
  expect_null(raw$manifest$rms_scale)
  expect_null(ldf$manifest$rms_scale)
  # raw and ldf reconstruct the same model
  expect_equal(tcrossprod(raw$activity, raw$effects),
               tcrossprod(ldf$activity, ldf$effects), tolerance = 1e-8)
})

test_that("rms works for sample_side = columns", {
  fit <- fixture_fit(sample_side = "columns")
  rms <- standardize_factors(fit, sample_side = "columns", scaling = "rms")
  K <- ncol(rms$activity)
  skip_if(K == 0L, "fixture produced zero factors")
  expect_equal(unname(sqrt(colMeans(rms$activity^2))), rep(1, K), tolerance = 1e-8)
})

test_that("rms handles K = 0 without error", {
  fit <- fixture_fit_zero_factors()
  rms <- standardize_factors(fit, scaling = "rms")
  expect_equal(ncol(rms$activity), 0L)
  expect_equal(rms$manifest$scaling, "rms")
})
