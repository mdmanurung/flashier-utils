association_fixture <- function(repeated=FALSE, K=4L) {
  set.seed(82)
  donors <- paste0("d",1:12)
  units <- if(repeated) rep(donors,each=2) else donors
  n <- length(units)
  md <- data.frame(sample=paste0("s",seq_len(n)),donor=units,
                   group=rep(c("A","B"),each=n/2),visit=if(repeated) rep(0:1,12) else rep(0,n))
  random <- matrix(rnorm(12*K,sd=2),12,K)
  A <- random[match(md$donor,donors),,drop=FALSE] + matrix(rnorm(n*K),n,K)
  A[,1] <- A[,1] + as.numeric(md$group=="B")*2 + md$visit
  dimnames(A) <- list(md$sample,paste0("F",seq_len(K)))
  B <- diag(K); dimnames(B) <- list(paste0("g",seq_len(K)),colnames(A))
  view <- .make_representation(A,B,list(fit_id="association_fixture",schema_version="1.0.0",sample_side="rows",native_indices=seq_len(K),scaling="raw",orientation="as_fit",constraints=c(activity="signed",effects="signed"),feature_scale="intensity"))
  list(view=view,metadata=md)
}
