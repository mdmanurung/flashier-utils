test_that("mash uses every common-basis row and native summary generics", {
  skip_if_not_installed("mashr")
  d<-association_fixture(K=8)
  set.seed(6)
  table<-expand.grid(factor=colnames(d$view$activity),cohort=c("A","B","C"),stringsAsFactors=FALSE)
  table$estimate<-rnorm(nrow(table));table$std.error<-rep(0.5,nrow(table));table$term<-"groupB";table$basis_id<-d$view$manifest$basis_id
  attr(table,"analysis_metadata")<-d$view$manifest
  expect_error(shrink_factor_effects(table,"cohort",term="groupB"),"Supply V")
  out<-shrink_factor_effects(table,"cohort",term="groupB",independent_conditions=TRUE)
  data<-mashr::mash_set_data(out$Bhat,out$Shat)
  native<-mashr::mash(data,mashr::cov_canonical(data),verbose=FALSE)
  expect_equal(out$table$estimate,as.vector(ashr::get_pm(native)))
  expect_equal(out$table$posterior.sd,as.vector(ashr::get_psd(native)))
  expect_equal(out$table$lfsr,as.vector(ashr::get_lfsr(native)))
  expect_error(shrink_factor_effects(table[-1,],"cohort",term="groupB",independent_conditions=TRUE),"every eligible")
  bad<-table;bad$basis_id[1]<-"wrong"
  expect_error(shrink_factor_effects(bad,"cohort",term="groupB",independent_conditions=TRUE),"incompatible")
})
