test_that("dream retains donor random effects and native test extraction", {
  skip_if_not_installed("variancePartition")
  d <- association_fixture(TRUE)
  rownames(d$metadata) <- d$metadata$sample
  formula <- ~group+visit+(1|donor)
  C <- matrix(c(0,1,-1),3,dimnames=list(colnames(model.matrix(~group+visit,d$metadata)),"groupMinusVisit"))
  out <- test_factors(d$view,d$metadata,formula,"dream",unit_col="donor",independence="repeated",contrasts=C,control=list(dream=list(ddf="Satterthwaite")))
  native <- variancePartition::eBayes(variancePartition::dream(t(d$view$activity),formula,d$metadata,L=C,ddf="Satterthwaite"))
  expect_equal(out$table$estimate,unname(drop(native$coefficients[,"groupMinusVisit"])),tolerance=1e-8)
  expect_equal(out$table$p.value,unname(drop(native$p.value[,"groupMinusVisit"])),tolerance=1e-8)
  expect_true(all(out$table$n_units==12))
  expect_equal(out$analysis_metadata$formula,formula)
})
