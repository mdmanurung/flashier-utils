test_that("native operations preserve identity and numerical native equivalence", {
  fit <- fixture_fit()
  reordered <- reorder_factors(fit,c("F2","F1"))
  native <- flashier::flash_factors_reorder(fit,2:1)
  expect_equal(reordered$L_pm,native$L_pm)
  expect_equal(standardize_factors(reordered)$manifest$factor_ids,c("F2","F1"))
  expect_equal(unname(stats::fitted(reordered)),unname(stats::fitted(fit)))
  removed <- remove_factors(fit,"F1")
  expect_equal(standardize_factors(removed)$manifest$factor_ids,"F2")
  expect_equal(removed$F_pm,flashier::flash_factors_remove(fit,1)$F_pm)
  expect_error(reorder_factors(fit,"F1"),"permutation")
  fixed <- fix_factors(fit,"F1","effects")
  refit <- refit_factors(fixed,control=list(maxiter=30L,verbose=0L))
  expect_equal(refit$F_pm[,1],fit$F_pm[,1],tolerance=1e-12)
  expect_false(identical(attr(refit,"flashbridge_metadata")$fit_id,attr(fit,"flashbridge_metadata")$fit_id))
  expect_s3_class(fix_factors(fixed,"F1","effects",action="unfix"),"flash")
  expect_error(fix_factors(fixed,"F1","effects",ids="g1",action="unfix"),"partial")
  for(side in c("rows","columns")) {
    f <- fixture_fit(side)
    fixed <- fix_factors(f,"F1","effects")
    native <- flashier::flash_factors_fix(f,1,which_dim=if(side=="rows") "factors" else "loadings")
    expect_equal(fixed$L_pm,native$L_pm)
    expect_equal(fixed$F_pm,native$F_pm)
  }
})

test_that("partial ids fix exactly the named entries in both orientations and axes", {
  for(side in c("rows","columns")) {
    f <- fixture_fit(side)
    view <- standardize_factors(f,scaling="raw")
    for(axis in c("activity","effects")) {
      all_ids <- rownames(view[[axis]])
      for(ids in list(all_ids[c(2,5)],all_ids[3])) {
        fixed <- fix_factors(f,"F1",axis,ids=ids)
        mode <- if((axis=="activity")==(side=="rows")) "loadings" else "factors"
        native <- flashier::flash_factors_fix(f,1,which_dim=mode,fixed_idx=all_ids %in% ids)
        expect_equal(flashier:::get.fix.idx(flashier:::get.fit(fixed))[[1]],which(all_ids %in% ids))
        expect_equal(fixed$L_pm,native$L_pm)
        expect_equal(fixed$F_pm,native$F_pm)
      }
    }
  }
  f <- fixture_fit()
  expect_error(fix_factors(f,"F1","effects",ids=character()),"at least one")
  expect_error(fix_factors(f,"F1","effects",ids="nope"),"Unknown")
  expect_error(fix_factors(f,NULL,"effects"),"factors is required")
  expect_error(remove_factors(f,NULL),"factors is required")
})

test_that("reorder gives a new identity so a stale view is refused", {
  fit <- fixture_fit()
  view <- standardize_factors(fit)
  reordered <- reorder_factors(fit,c("F2","F1"))
  expect_false(identical(attr(reordered,"flashbridge_metadata")$fit_id,attr(fit,"flashbridge_metadata")$fit_id))
  expect_error(factor_uncertainty(reordered,"effects",representation=view),"different fit")
  expect_s3_class(factor_uncertainty(reordered,"effects",representation=standardize_factors(reordered)),"data.frame")
})

test_that("provenance engine fields are single strings", {
  meta <- attr(fixture_fit(),"analysis_metadata")
  expect_type(meta[["engine"]],"character")
  expect_type(meta[["engine_version"]],"character")
  expect_length(meta[["engine_version"]],1L)
})

test_that("nonnegative wrapper agrees with bounded native point-exponential fit", {
  X <- abs(fixture_matrix())
  wrapped <- fit_nonnegative_ebmf(X,"rows",max_factors=2L,verbose=0L,seed=2)
  set.seed(2)
  native <- flashier::flash(X,ebnm_fn=ebnm::ebnm_point_exponential,greedy_Kmax=2L,verbose=0L)
  expect_equal(wrapped$L_pm,native$L_pm,tolerance=1e-10)
  expect_true(all(wrapped$L_pm >= 0)); expect_true(all(wrapped$F_pm >= 0))
  expect_lte(wrapped$n_factors,2)
  expect_error(fit_nonnegative_ebmf(X,"rows",initialization="NNLM"),"unsupported|unavailable")
})
