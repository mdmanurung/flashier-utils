test_that("limma moderation and transposition are native-equivalent", {
  skip_if_not_installed("limma")
  d <- association_fixture()
  out <- test_factor_association(d$view, d$metadata, ~group, "limma", unit_col = "donor", independence = "independent")
  native <- limma::lmFit(t(d$view$activity), model.matrix(~group, d$metadata))
  C <- matrix(c(0, 1), 2, dimnames = list(colnames(native$coefficients), "groupB"))
  native <- limma::eBayes(limma::contrasts.fit(native, C))
  expect_equal(out$table$estimate, unname(drop(native$coefficients)))
  expect_equal(out$table$std.error, unname(drop(native$stdev.unscaled) * sqrt(native$s2.post)))
  expect_equal(out$table$p.value, unname(drop(native$p.value)))
  expect_equal(out$table$adj.p.value, p.adjust(unname(drop(native$p.value)), "BH"))
})

test_that("a single factor is tested, and equals lm at K = 1", {
  skip_if_not_installed("limma")
  d <- association_fixture()
  args <- list(d$view, d$metadata, ~group, unit_col = "donor", independence = "independent", factors = "F1")
  limma_out <- do.call(test_factor_association, c(args, engine = "limma"))
  lm_out <- do.call(test_factor_association, c(args, engine = "lm"))
  expect_false(anyNA(limma_out$table$p.value))
  expect_equal(limma_out$table$estimate, lm_out$table$estimate)
  expect_equal(limma_out$table$std.error, lm_out$table$std.error)
  expect_equal(limma_out$table$p.value, lm_out$table$p.value)
})

test_that("factors with very different residual variance scales trigger the moderation warning", {
  skip_if_not_installed("limma")
  d <- association_fixture()
  activity <- d$view$activity
  activity[, "F2"] <- activity[, "F2"] * 20
  view <- .make_representation(activity, d$view$effects, d$view$manifest)
  expect_warning(
    test_factor_association(view, d$metadata, ~group, "limma", unit_col = "donor", independence = "independent"),
    "differ more than tenfold"
  )
})
