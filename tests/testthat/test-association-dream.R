test_that("dream retains donor random effects and native test extraction", {
  skip_if_not_installed("variancePartition")
  d <- association_fixture(TRUE)
  rownames(d$metadata) <- d$metadata$sample
  formula <- ~ group + visit + (1 | donor)
  C <- matrix(c(0, 1, -1), 3, dimnames = list(colnames(model.matrix(~ group + visit, d$metadata)), "groupMinusVisit"))
  out <- test_factor_association(
    d$view,
    d$metadata,
    formula,
    "dream",
    unit_col = "donor",
    independence = "repeated",
    contrasts = C,
    control = list(dream = list(ddf = "Satterthwaite"))
  )
  fitted <- variancePartition::dream(t(d$view$activity), formula, d$metadata, L = C, ddf = "Satterthwaite")
  native <- variancePartition::eBayes(fitted)
  # default: dream's own Satterthwaite t, df and p, with no moderation
  expect_equal(out$table$estimate, unname(drop(fitted$coefficients[, "groupMinusVisit"])), tolerance = 1e-8)
  expect_equal(out$table$p.value, unname(drop(fitted$p.value[, "groupMinusVisit"])), tolerance = 1e-8)
  expect_equal(out$table$df, unname(drop(fitted$df.residual[, "groupMinusVisit"])), tolerance = 1e-8)
  expect_identical(out$analysis_metadata$moderation, "none")
  expect_true(all(out$table$n_units == 12))
  expect_equal(out$analysis_metadata$formula, formula)
  # moderate = TRUE restores the eBayes path
  moderated <- test_factor_association(
    d$view,
    d$metadata,
    formula,
    "dream",
    unit_col = "donor",
    independence = "repeated",
    contrasts = C,
    control = list(dream = list(ddf = "Satterthwaite")),
    moderate = TRUE
  )
  expect_equal(moderated$table$p.value, unname(drop(native$p.value[, "groupMinusVisit"])), tolerance = 1e-8)
  expect_identical(moderated$analysis_metadata$moderation, "ebayes")
})

test_that("eBayes pools residual df across factors, which inflates the df of a between-donor term", {
  skip_if_not_installed("variancePartition")
  d <- association_fixture(TRUE)
  formula <- ~ group + (1 | donor)
  run <- function(moderate) {
    test_factor_association(
      d$view,
      d$metadata,
      formula,
      "dream",
      unit_col = "donor",
      independence = "repeated",
      moderate = moderate
    )$table
  }
  unmoderated <- run(FALSE)
  moderated <- run(TRUE)
  between <- unmoderated$term == "groupB"
  # 12 donors, 2 per group level pair: the between-donor df is about 10, not hundreds
  expect_true(all(unmoderated$df[between] < 15))
  expect_true(all(moderated$df[moderated$term == "groupB"] > 100))
})

test_that("a dream random slope without a unit intercept is rejected", {
  skip_if_not_installed("variancePartition")
  d <- association_fixture(TRUE)
  expect_error(
    test_factor_association(
      d$view,
      d$metadata,
      ~ group + (0 + visit | donor),
      "dream",
      unit_col = "donor",
      independence = "repeated"
    ),
    "random intercept for the declared biological unit"
  )
})

test_that("control$ebayes needs moderate = TRUE for dream", {
  skip_if_not_installed("variancePartition")
  d <- association_fixture(TRUE)
  expect_error(
    test_factor_association(
      d$view,
      d$metadata,
      ~ group + (1 | donor),
      "dream",
      unit_col = "donor",
      independence = "repeated",
      control = list(ebayes = list(robust = TRUE))
    ),
    "moderate = TRUE"
  )
})
