fixture_matrix <- function() {
  set.seed(55)
  A <- matrix(rnorm(24),12,2)
  B <- matrix(c(2,-1,0,1,3,-2,4,0,1,2,-3,0,-2,1,0,4),8,2)
  X <- tcrossprod(A,B) + matrix(rnorm(96,sd=0.2),12,8)
  dimnames(X) <- list(paste0("s",1:12),paste0("g",1:8))
  X
}
fixture_fit <- function(sample_side="rows", K=2L) {
  X <- fixture_matrix()
  prior <- list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace)
  if (sample_side=="columns") { X <- t(X); prior <- rev(prior) }
  fit_ebmf(X, sample_side, ebnm_fn=prior, max_factors=K, nullcheck=FALSE,
           backfit=TRUE, verbose=0L, seed=12)
}
fixture_view <- function() {
  A <- matrix(c(1,2,4,2,0,-1),3,dimnames=list(c("s1","s2","s3"),c("F1","F2")))
  B <- matrix(c(2,-1,0,1,3,-2,4,0),4,dimnames=list(paste0("g",1:4),colnames(A)))
  .make_representation(A,B,list(fit_id="exact_fixture",schema_version="1.0.0",
       sample_side="rows", native_indices=1:2, scaling="raw", orientation="as_fit",
       constraints=c(activity="signed",effects="signed"),feature_scale="arcsinh_intensity"))
}

## SYNC-0 contract fixtures: signed, tied-anchor, zero-activity, zero-factor.

## Signed loadings on both sides (positive and negative entries), so sign
## anchoring and signed-prior constraints are exercised.
fixture_view_signed <- function() {
  A <- matrix(c(1,-2,4,-2,0,-1),3,dimnames=list(c("s1","s2","s3"),c("F1","F2")))
  B <- matrix(c(-2,-1,0,1,3,-2,-4,0),4,dimnames=list(paste0("g",1:4),colnames(A)))
  .make_representation(A,B,list(fit_id="signed_fixture",schema_version="1.0.0",
       sample_side="rows", native_indices=1:2, scaling="raw", orientation="as_fit",
       constraints=c(activity="signed",effects="signed"),feature_scale="arcsinh_intensity"))
}

## Two genes tie for the maximal absolute loading on F1 (g1 and g3 both -2/2),
## so the deterministic lexical tie-break in anchor selection is exercised.
fixture_view_tied_anchor <- function() {
  A <- matrix(c(1,2,4,2,0,-1),3,dimnames=list(c("s1","s2","s3"),c("F1","F2")))
  B <- matrix(c(2,-1,0,1,-2,-2,4,0),4,dimnames=list(paste0("g",1:4),colnames(A)))
  .make_representation(A,B,list(fit_id="tied_fixture",schema_version="1.0.0",
       sample_side="rows", native_indices=1:2, scaling="raw", orientation="as_fit",
       constraints=c(activity="signed",effects="signed"),feature_scale="arcsinh_intensity"))
}

## F2 has zero activity on the sample side, exercising the RMS fallback (s_k = 1)
## and zero-activity handling in standardization. A is 3x2 column-major, so
## F1 = (1,2,4) and F2 = (0,0,0).
fixture_view_zero_activity <- function() {
  A <- matrix(c(1,2,4,0,0,0),3,dimnames=list(c("s1","s2","s3"),c("F1","F2")))
  B <- matrix(c(2,-1,0,1,3,-2,4,0),4,dimnames=list(paste0("g",1:4),colnames(A)))
  .make_representation(A,B,list(fit_id="zeroact_fixture",schema_version="1.0.0",
       sample_side="rows", native_indices=1:2, scaling="raw", orientation="as_fit",
       constraints=c(activity="signed",effects="signed"),feature_scale="arcsinh_intensity"))
}

## A fit with K = 0 factors (max_factors = 0), exercising the zero-factor path.
fixture_fit_zero_factors <- function() {
  X <- fixture_matrix()
  fit_ebmf(X, "rows", ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
           max_factors=0L, nullcheck=FALSE, backfit=FALSE, verbose=0L, seed=12)
}
