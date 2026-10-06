test_that("base correlations match native output and retain constants without p-values", {
  d<-association_fixture()
  out<-factor_cor(d$view,"pearson")
  expect_equal(out$correlation,cor(d$view$activity))
  expect_equal(factor_cor(d$view,"spearman")$correlation,cor(d$view$activity,method="spearman"))
  expect_true(all(out$pair_counts==nrow(d$view$activity)))
  d$view$activity[,1]<-1;d$view<-.make_representation(d$view$activity,d$view$effects,d$view$manifest)
  out<-factor_cor(d$view,"pearson")
  expect_true(all(is.na(out$correlation[1,])))
  expect_equal(out$constant_factors,"F1")
  expect_false("p.value" %in% names(out))
})

test_that("CorShrink delegates only the audited public method", {
  if(!nzchar(system.file(package="CorShrink"))) skip("CorShrink unavailable")
  # The pinned package's CVXR/stats sd import warning is recorded in the audit.
  # Only this known namespace notice is muffled; fitting warnings stay visible.
  withCallingHandlers(requireNamespace("CorShrink"),warning=function(w) {
    if(grepl("replacing previous import.*CVXR::sd",conditionMessage(w))) invokeRestart("muffleWarning")
  })
  d<-association_fixture(K=6)
  previous_warn<-getOption("warn")
  out<-factor_cor(d$view,"corshrink")
  expect_identical(getOption("warn"),previous_warn)
  expect_error(factor_cor(d$view,"corshrink",control=list(dosym=TRUE)),"native branch")
  native<-CorShrink::CorShrinkData(d$view$activity,image="null",sd_boot=FALSE,type="cor")
  options(warn=previous_warn)
  expect_equal(unname(out$correlation),unname(native$cor),tolerance=1e-8)
  expect_true(all(out$pair_counts==nrow(d$view$activity)))
  expect_true(is.finite(out$analysis_metadata$matrix_properties$minimum_eigenvalue_of_symmetric_part))
})
