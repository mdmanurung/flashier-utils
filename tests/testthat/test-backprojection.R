test_that("named point reconstruction agrees with hand arithmetic and rejects wrong bases", {
  view <- fixture_view(); beta <- c(F1=2,F2=-1)
  out <- backproject_contrast(view,beta,"coefficients",factor_basis=view$manifest)
  expect_equal(out$effects$estimate,c(1,0,-4,2))
  expect_equal(out$analysis_metadata$feature_scale,"arcsinh_intensity")
  expect_false("p.value" %in% names(out$effects))
  expect_error(backproject_contrast(view,unname(beta),"coefficients",factor_basis=view$manifest),"IDs")
  expect_error(backproject_contrast(view,beta,"coefficients"),"factor_basis")
  other <- view$manifest; other$fit_id <- "different"
  expect_error(backproject_contrast(view,beta,"coefficients",factor_basis=other),"different fit")
  partial <- backproject_contrast(view,beta[1],"coefficients",factor_basis=view$manifest,allow_partial=TRUE)
  expect_equal(partial$analysis_metadata$omitted_factors,"F2")
  expect_equal(partial$effects$estimate,c(4,-2,0,2))
  expect_error(backproject_contrast(view,beta[1],"coefficients",factor_basis=view$manifest),"Missing factors")
})

test_that("matching scale sign and permutation transforms preserve feature contrasts", {
  view <- fixture_view(); beta <- c(F1=2,F2=-1)
  for(q in list(c(3,0.25),c(-1,1))) {
    changed <- .make_representation(sweep(view$activity,2,q,"*"),sweep(view$effects,2,q,"/"),view$manifest)
    out <- backproject_contrast(changed,beta*q,"coefficients",factor_basis=changed$manifest)
    expect_equal(out$effects$estimate,c(1,0,-4,2))
    expect_error(backproject_contrast(changed,beta,"coefficients",factor_basis=view$manifest),"Incompatible")
  }
  changed <- .make_representation(view$activity[,2:1],view$effects[,2:1],view$manifest)
  out <- backproject_contrast(changed,beta[2:1],"coefficients",factor_basis=changed$manifest)
  expect_equal(out$effects$estimate,c(1,0,-4,2))
})
