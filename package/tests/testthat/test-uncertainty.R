test_that("posterior SD uses absolute fixed multipliers and native PVE", {
  fit <- fixture_fit()
  for (location in c("effects","activity")) for (side in c("activity","effects")) {
    view <- standardize_factors(fit,d_location=location,orientation="anchor_feature")
    out <- factor_uncertainty(fit,side,representation=view)
    raw <- if(side=="activity") fit$L_psd else fit$F_psd
    expected <- sweep(raw,2,abs(view$manifest[[paste0(side,"_multiplier")]]),"*")
    expect_equal(out$posterior.sd,as.vector(expected),tolerance=1e-10)
    expect_false("std.error" %in% names(out))
  }
  expect_equal(factor_pve(fit)$pve,fit$pve)
  missing <- fit; missing$L_psd <- NULL
  out <- factor_uncertainty(missing,"activity")
  expect_true(all(is.na(out$posterior.sd)))
  expect_error(factor_uncertainty(missing,"activity",require_uncertainty=TRUE),"unavailable")
})
