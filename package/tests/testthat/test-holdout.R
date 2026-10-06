holdout_callback <- function(X,metadata,sample_side,seed,settings,observation_aux=list()) {
  priors<-list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace)
  if(sample_side=="columns") priors<-rev(priors)
  fit_ebmf(X,sample_side,ebnm_fn=priors,max_factors=2L,backfit=TRUE,nullcheck=FALSE,seed=seed,verbose=0,feature_scale="intensity")
}

test_that("held-out fixed programs stay in their raw training basis in both orientations", {
  X<-fixture_matrix();md<-data.frame(sample=rownames(X),donor=rep(paste0("d",1:6),each=2),cohort=rep(c("A","B"),each=6))
  folds<-data.frame(sample=md$sample,fold=c(rep("train",8),rep("test",4)))
  for(side in c("rows","columns")) {
    data<-if(side=="rows") X else t(X)
    out<-factor_holdout(data,md,"donor",folds,holdout_callback,seed=12,control=list(sample_side=side,preprocessing_manifest=list(scope="identity"),evaluated_folds="test",cohort_col="cohort"))
    expect_equal(nrow(out$failures),0)
    expect_equal(rownames(out$activities$test),md$sample[9:12])
    expect_setequal(out$fold_results$test$training_ids,md$sample[1:8])
    expect_equal(attr(out$activities$test,"analysis_metadata")$training_basis_id,out$fold_results$test$train_basis$basis_id)
    expect_true(all(out$metrics$evaluation_type=="projection"))
    expect_true("macro_unit" %in% out$metrics$weighting)
    expect_true("training_mean" %in% out$metrics$model)
  }
  bad<-folds;bad$fold[1]<-"test"
  expect_error(factor_holdout(X,md,"donor",bad,holdout_callback,seed=1,control=list(sample_side="rows",preprocessing_manifest=list(scope="identity"))),"cross folds")
})

test_that("altering only withheld truth leaves masked activities unchanged", {
  X<-fixture_matrix();md<-data.frame(sample=rownames(X),donor=rep(paste0("d",1:6),each=2))
  folds<-data.frame(sample=md$sample,fold=c(rep("train",8),rep("test",4)))
  masked<-matrix(FALSE,12,8,dimnames=dimnames(X));masked[9:12,c(1,3)]<-TRUE
  settings<-list(sample_side="rows",preprocessing_manifest=list(scope="identity"),evaluated_folds="test")
  before<-factor_holdout(X,md,"donor",folds,holdout_callback,evaluation="masked",mask=masked,seed=4,control=settings)
  absurd<-X;absurd[masked]<-1e8
  after<-factor_holdout(absurd,md,"donor",folds,holdout_callback,evaluation="masked",mask=masked,seed=4,control=settings)
  expect_equal(nrow(before$failures),0)
  expect_equal(unname(before$activities$test),unname(after$activities$test),ignore_attr=TRUE,tolerance=1e-12)
  expect_equal(before$fold_results$test$effects,after$fold_results$test$effects)
  expect_true(all(before$metrics$evaluation_type=="masked"))
  expect_error(factor_holdout(X,md,"donor",folds,holdout_callback,seed=1,freeze_noise=TRUE),"unsupported")
})

 test_that("training preprocessing receives no held-out entries and masking precedes transform", {
  X<-fixture_matrix();md<-data.frame(sample=rownames(X),donor=rep(paste0("d",1:6),each=2))
  folds<-data.frame(sample=md$sample,fold=c(rep("train",8),rep("test",4)))
  preprocessing<-list(fit=function(X,metadata,training_ids,seed,settings) {
    stopifnot(setequal(rownames(X),training_ids),all(training_ids %in% paste0("s",1:8)))
    mean<-colMeans(X)
    list(X=sweep(X,2,mean,"-"),metadata=metadata,manifest=list(scope="training_only",center=mean),state=mean)
  },transform=function(X,metadata,state,seed,settings) list(X=sweep(X,2,state,"-"),metadata=metadata))
  settings<-list(sample_side="rows",evaluated_folds="test")
  before<-factor_holdout(X,md,"donor",folds,holdout_callback,preprocess_fn=preprocessing,evaluation="masked",activity_features=paste0("g",2:8),scoring_features="g1",seed=4,control=settings)
  absurd<-X;absurd[9:12,1]<-1e8
  after<-factor_holdout(absurd,md,"donor",folds,holdout_callback,preprocess_fn=preprocessing,evaluation="masked",activity_features=paste0("g",2:8),scoring_features="g1",seed=4,control=settings)
  expect_equal(nrow(before$failures),0)
  expect_equal(unname(before$activities$test),unname(after$activities$test),ignore_attr=TRUE,tolerance=1e-12)
})

test_that("named observation errors follow dense folds and sparse limits are explicit", {
  X<-fixture_matrix();md<-data.frame(sample=rownames(X),donor=rep(paste0("d",1:6),each=2))
  folds<-data.frame(sample=md$sample,fold=c(rep("train",8),rep("test",4)))
  S<-setNames(seq(0.2,0.4,length.out=12),rownames(X))
  for(side in c("rows","columns")) {
    data<-if(side=="rows")X else t(X)
    callback<-function(X,metadata,sample_side,seed,settings,observation_aux=list()) {
      samples<-if(sample_side=="rows")rownames(X) else colnames(X)
      stopifnot(identical(names(observation_aux$S),samples),all(observation_aux$S==S[samples]))
      priors<-list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace)
      if(sample_side=="columns")priors<-rev(priors)
      fit_ebmf(X,sample_side,S=observation_aux$S,max_factors=2,ebnm_fn=priors,seed=seed,backfit=TRUE,nullcheck=FALSE,verbose=0,control=list(S_dim=if(sample_side=="rows")1L else 2L),feature_scale="intensity")
    }
    control<-list(sample_side=side,preprocessing_manifest=list(scope="identity"),evaluated_folds="test",observation_aux=list(S=S))
    out<-factor_holdout(data,md,"donor",folds,callback,evaluation="masked",activity_features=paste0("g",1:6),scoring_features=c("g7","g8"),seed=3,control=control)
    expect_equal(nrow(out$failures),0L,info=paste(out$failures$reason,collapse=";"))
    expect_true(nrow(out$metrics)>0)
    expect_equal(rownames(out$activities$test),md$sample[9:12])
    sparse<-factor_holdout(Matrix::Matrix(data,sparse=TRUE),md,"donor",folds,callback,evaluation="masked",activity_features=paste0("g",1:6),scoring_features=c("g7","g8"),seed=3,control=control)
    expect_match(sparse$failures$reason,"Sparse observations with explicit S")
    supported<-factor_holdout(Matrix::Matrix(data,sparse=TRUE),md,"donor",folds,holdout_callback,evaluation="masked",activity_features=paste0("g",1:6),scoring_features=c("g7","g8"),seed=3,control=control[names(control)!="observation_aux"])
    expect_equal(nrow(supported$failures),0L,info=paste(supported$failures$reason,collapse=";"))
  }
})
