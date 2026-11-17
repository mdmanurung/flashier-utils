test_that("lm contrasts, SE and joint covariance match native engines", {
  d <- association_fixture()
  out <- test_factor_association(
    d$view,
    d$metadata,
    ~group,
    "lm",
    unit_col = "donor",
    independence = "independent",
    control = list(joint_coefficient_covariance = TRUE)
  )
  native <- lm(d$view$activity[, 1] ~ group, data = d$metadata)
  expect_equal(out$table$estimate[1], unname(coef(native)["groupB"]))
  expect_equal(out$table$std.error[1], unname(summary(native)$coefficients["groupB", "Std. Error"]))
  expect_equal(out$table$adj.p.value, p.adjust(out$table$p.value, "BH"))
  joint <- lm(d$view$activity ~ group, data = d$metadata)
  V <- vcov(joint)
  index <- match(paste0(colnames(d$view$activity), ":groupB"), rownames(V))
  expect_equal(unname(out$coefficient_covariance$groupB), unname(V[index, index]))
  reverse <- matrix(c(0, -1), 2, dimnames = list(c("(Intercept)", "groupB"), "AminusB"))
  reversed <- test_factor_association(
    d$view,
    d$metadata,
    ~group,
    "lm",
    unit_col = "donor",
    independence = "independent",
    contrasts = reverse
  )
  expect_equal(reversed$table$estimate, -out$table$estimate)
  expect_equal(attr(out$table, "analysis_metadata")$basis_id, d$view$manifest$basis_id)
})

test_that("a non-estimable lm row has no p-value and stays out of the BH family", {
  d <- association_fixture()
  # a zero standard error on F2 only: statistic is Inf and the old p-value was exactly 0
  real_vcov <- stats::vcov
  zero_for_f2 <- function(object, ...) {
    V <- real_vcov(object, ...)
    if (isTRUE(all.equal(unname(object$model[[1]]), unname(d$view$activity[, "F2"])))) {
      V[] <- 0
    }
    V
  }
  testthat::local_mocked_bindings(vcov = zero_for_f2, .package = "stats")
  out <- test_factor_association(
    d$view,
    d$metadata,
    ~group,
    "lm",
    unit_col = "donor",
    independence = "independent"
  )
  bad <- out$table$status == "non_estimable"
  expect_identical(out$table$factor[bad], "F2")
  expect_true(all(is.na(out$table$p.value[bad])))
  expect_true(all(is.na(out$table$adj.p.value[bad])))
  expect_true(all(nzchar(out$table$reason[bad])))
  ok <- !bad
  expect_equal(out$table$adj.p.value[ok], p.adjust(out$table$p.value[ok], "BH"))
})
