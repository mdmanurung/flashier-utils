test_that("design choices cannot silently introduce pseudoreplication", {
  d <- association_fixture()
  expect_error(test_factor_association(d$view, d$metadata, ~group), "explicitly")
  repeated <- association_fixture(TRUE)
  expect_error(
    test_factor_association(
      repeated$view,
      repeated$metadata,
      ~group,
      "lm",
      unit_col = "donor",
      independence = "independent"
    ),
    "repeated records"
  )
  expect_error(
    test_factor_association(
      repeated$view,
      repeated$metadata,
      ~ group + (1 | donor),
      "lm",
      unit_col = "donor",
      independence = "independent"
    ),
    "Random effects"
  )
  md <- d$metadata
  md$duplicate <- md$group
  expect_error(
    test_factor_association(d$view, md, ~ group + duplicate, "lm", unit_col = "donor", independence = "independent"),
    "rank deficient"
  )
  md$group[1] <- NA
  expect_error(
    test_factor_association(d$view, md, ~group, "lm", unit_col = "donor", independence = "independent"),
    "Missing metadata"
  )
  out <- test_factor_association(
    d$view,
    md,
    ~group,
    "lm",
    unit_col = "donor",
    independence = "independent",
    na_action = "omit"
  )
  expect_true(all(out$table$n_samples == 11))
  expect_equal(out$analysis_metadata$exclusions$sample, "s1")
})

test_that("using the sample ID as the unit is announced", {
  d <- association_fixture()
  expect_message(
    test_factor_association(d$view, d$metadata, ~group, "lm", unit_col = "sample", independence = "independent"),
    "every sample is treated as its own biological unit"
  )
})

test_that("intercept tests do not change the adjusted p-values of the other terms", {
  d <- association_fixture()
  args <- list(d$view, d$metadata, ~group, "lm", unit_col = "donor", independence = "independent")
  without <- do.call(test_factor_association, args)$table
  with_intercept <- do.call(test_factor_association, c(args, include_intercept = TRUE))$table
  keep <- with_intercept$term != "(Intercept)"
  expect_equal(with_intercept$adj.p.value[keep], without$adj.p.value)
})
