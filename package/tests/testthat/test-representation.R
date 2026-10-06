test_that("both orientations and D placements reconstruct native fitted means", {
  for (side in c("rows","columns")) {
    fit <- fixture_fit(side)
    original <- serialize(fit,NULL)
    expected <- stats::fitted(fit)
    if(side=="columns") expected <- t(expected)
    for (scaling in c("raw","ldf")) for(type in c("f","o","i")) for(location in c("activity","effects")) {
      view <- standardize_factors(fit,scaling=scaling,type=type,d_location=location)
      expect_equal(unname(tcrossprod(view$activity,view$effects)),unname(expected),tolerance=1e-10)
      expect_identical(rownames(view$activity),paste0("s",1:12))
    }
    anchor <- standardize_factors(fit,orientation="anchor_feature")
    for(k in seq_len(ncol(anchor$effects))) expect_gte(anchor$effects[anchor$manifest$orientation_anchor[k],k],0)
    expect_identical(serialize(fit,NULL),original)
    expect_error(standardize_factors(fit,orientation="anchor_feature",prior_support=c(activity="nonnegative",effects="signed")),"signed prior")
  }
  for(K in c(0L,1L)) {
    view <- standardize_factors(fixture_fit(K=K))
    expect_equal(dim(view$activity),c(12,K))
    expect_equal(dim(view$effects),c(8,K))
    expect_equal(nrow(view$factors),K)
  }
  bare <- flashier::flash_init(fixture_matrix())
  expect_error(standardize_factors(bare),"sample_side")
  expect_equal(dim(standardize_factors(bare,"rows")$effects),c(8,0))
})
