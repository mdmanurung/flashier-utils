test_that("ECDF display retains signs, reference ties and excludes inference", {
  view<-fixture_view()
  reference<-c("s1","s2")
  out<-factor_activity_percentile(view,reference,factors="F2")
  expect_equal(as.numeric(out),100*ecdf(view$activity[reference,"F2"])(view$activity[,"F2"]))
  expect_equal(attr(out,"analysis_metadata")$reference_sample_ids,reference)
  expect_true(attr(out,"analysis_metadata")$display_only)
  expect_error(factor_cor(out,"pearson"),"Display-only")
  expect_error(backproject_contrast(out,c(F1=1,F2=0),input="coefficients",factor_basis=view$manifest),"Display-only")
  expect_error(factor_activity_percentile(view,"absent"),"reference sample")
  tied<-association_fixture()$view
  tied$activity[1:2,1]<-0
  tied<-.make_representation(tied$activity,tied$effects,tied$manifest)
  result<-factor_activity_percentile(tied,c("s1","s2"),"F1")
  expect_equal(as.numeric(result[1:2,]),c(100,100))
})
