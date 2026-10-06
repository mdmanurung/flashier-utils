test_that("design choices cannot silently introduce pseudoreplication", {
  d <- association_fixture()
  expect_error(test_factors(d$view,d$metadata,~group),"explicitly")
  repeated <- association_fixture(TRUE)
  expect_error(test_factors(repeated$view,repeated$metadata,~group,"lm",unit_col="donor",independence="independent"),"repeated records")
  expect_error(test_factors(repeated$view,repeated$metadata,~group+(1|donor),"lm",unit_col="donor",independence="independent"),"Random effects")
  md <- d$metadata; md$duplicate <- md$group
  expect_error(test_factors(d$view,md,~group+duplicate,"lm",unit_col="donor",independence="independent"),"rank deficient")
  md$group[1] <- NA
  expect_error(test_factors(d$view,md,~group,"lm",unit_col="donor",independence="independent"),"Missing metadata")
  out <- test_factors(d$view,md,~group,"lm",unit_col="donor",independence="independent",na_action="omit")
  expect_true(all(out$table$n_samples==11))
  expect_equal(out$analysis_metadata$exclusions$sample,"s1")
})
