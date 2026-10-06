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
