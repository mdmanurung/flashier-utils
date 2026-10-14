test_that("donor occurrences preserve visits, auxiliary maps and serial/parallel identities", {
  md <- data.frame(sample=paste0("s",1:8),donor=rep(c("d1","d2","d3"),c(2,3,3)),visit=c(1,2,1,2,3,1,2,3),cohort=c(rep("A",5),rep("B",3)))
  out <- .bootstrap_units(md,"donor","cohort",seed=7)
  expect_false(anyDuplicated(out$metadata$sample)>0)
  for(unit in unique(out$mapping$unit_occurrence)) {
    rows<-out$mapping[out$mapping$unit_occurrence==unit,]
    expect_setequal(rows$original_sample,md$sample[md$donor==rows$original_unit[1]])
  }
  X <- matrix(seq_len(16),8,dimnames=list(md$sample,c("g1","g2")))
  aux <- .resample_aux(list(S=setNames(seq_len(8),md$sample)),out$mapping,X,"rows")
  expect_equal(unname(aux$S),match(out$mapping$original_sample,md$sample))
  expect_equal(names(aux$S),out$mapping$sample)
  fit <- fixture_fit(); md <- data.frame(sample=paste0("s",1:12),donor=rep(paste0("d",1:6),each=2))
  callback <- function(X,metadata,sample_side,seed,settings,observation_aux=list()) fit_ebmf(X,sample_side,max_factors=1L,backfit=TRUE,seed=seed,verbose=0L,feature_scale="unknown")
  # Stability cannot compare unknown feature units: use the same declared units.
  callback <- function(X,metadata,sample_side,seed,settings,observation_aux=list()) fit_ebmf(X,sample_side,max_factors=1L,backfit=TRUE,seed=seed,verbose=0L,feature_scale="intensity")
  attr(fit,"flashbridge_metadata")$feature_scale <- "intensity"
  dir <- tempfile();on.exit(unlink(dir,recursive=TRUE))
  first <- factor_stability(fixture_matrix(),fit,md,"donor",B=2L,fit_fn=callback,seed=33,min_congruence=0.5,checkpoint_dir=dir,control=list(callback_configuration=list(prior="point_normal",K=1L)))
  second <- factor_stability(fixture_matrix(),fit,md,"donor",B=2L,fit_fn=callback,seed=33,min_congruence=0.5,checkpoint_dir=dir,control=list(callback_configuration=list(prior="point_normal",K=1L)))
  expect_equal(first$resamples,second$resamples)
  expect_equal(first$summary,second$summary)
  expect_equal(first$replicates,second$replicates)
  if(.Platform$OS.type=="unix") {
    parallel <- factor_stability(fixture_matrix(),fit,md,"donor",B=2L,fit_fn=callback,seed=33,min_congruence=0.5,workers=2L,control=list(callback_configuration=list(prior="point_normal",K=1L)))
    expect_equal(first$resamples,parallel$resamples)
    expect_equal(first$matches,parallel$matches)
  }
  expect_error(factor_stability(fixture_matrix(),fit,md,"donor",B=2L,fit_fn=callback,seed=34,checkpoint_dir=dir,control=list(callback_configuration=list(prior="point_normal",K=1L))),"mismatch")
  failed <- factor_stability(fixture_matrix(),fit,md,"donor",B=2L,fit_fn=function(...) stop("deliberate fitting failure"),seed=3,min_congruence=0.5)
  zero_callback <- function(X,metadata,sample_side,seed,settings,observation_aux=list()) fit_ebmf(X,sample_side,max_factors=0L,seed=seed,verbose=0L,feature_scale="intensity")
  zero <- factor_stability(fixture_matrix(),fit,md,"donor",B=1L,fit_fn=zero_callback,seed=3,min_congruence=0.5)
  expect_equal(zero$replicates$rank,0)
  expect_true(all(zero$summary$recovery_rate_all_requested==0))
  expect_equal(nrow(failed$failures),2)
  expect_true(all(failed$summary$recovery_rate_all_requested==0))
})

test_that("donor strata cannot hide different metadata tuples", {
  metadata <- data.frame(sample = c("s1", "s2"), donor = "d1",
    first = c("x.y", "x"), second = c("z", "y.z"))
  expect_error(.bootstrap_units(metadata, "donor", c("first", "second"), seed = 7), "crosses strata")
})

test_that("automatic resampling preserves ordinary recipes and guards edited fits", {
  fit <- fixture_fit()
  attr(fit, "flashbridge_metadata")$feature_scale <- "intensity"
  metadata <- data.frame(sample = paste0("s", 1:12), donor = rep(paste0("d", 1:6), each = 2))
  run <- function(reference, ...) factor_stability(fixture_matrix(), reference, metadata,
    "donor", B = 1L, seed = 7, ...)
  reordered <- reorder_factors(fit, c("F2", "F1"))
  expect_identical(attr(reordered, "flashbridge_metadata")$refit_policy, "ordinary")
  expect_no_error(run(reordered))
  removed <- remove_factors(fit, "F2")
  expect_error(run(removed), "explicit fit_fn")
  fixed <- fix_factors(fit, "F1", "effects")
  expect_error(run(fixed), "explicit fit_fn")
  reordered_edit <- reorder_factors(fixed, c("F2", "F1"))
  expect_error(run(reordered_edit), "explicit fit_fn")
  projected <- suppressWarnings(project_factor_activity(fit, fixture_matrix(), "rows", "intensity",
    list(scope = "identity"), maxiter = 200))
  expect_error(run(projected$fit), "explicit fit_fn")
  custom <- run(removed, fit_fn = function(...) stop("deliberate custom failure"))
  expect_equal(nrow(custom$failures), 1L)
})
