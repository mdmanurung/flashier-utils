# Exact pinned helpers, MIT, Peter Carbonetto and Matthew Stephens 2025.
# commit df62f8d169fcac07cd0cb294a3c15051f573ccd5, R/rank_effects.R,
# R/compute_le_effects.R. Notice: ../../inst/licenses/singlecelljamboreeR.txt.
upstream_rank_effects <- function (effects_matrix) {
  if (!(is.matrix(effects_matrix) & is.numeric(effects_matrix)))
    stop("Input \"effects_matrix\" should be a numeric matrix")
  return(apply(-effects_matrix,2,rank))
}

upstream_compute_le_effects <- function (effects_matrix) {
  if (!(is.matrix(effects_matrix) & is.numeric(effects_matrix)))
    stop("Input \"effects_matrix\" should be a numeric matrix")
  k <- ncol(effects_matrix)
  if (k <= 1)
    return(effects_matrix)
  out <- effects_matrix

  # Repeat for each column of the effects matrix.
  for (j in 1:k) {
    x <- effects_matrix[,j]
    i <- which(x > 0)
    if (length(i) > 0) {
      y <- cbind(0,effects_matrix[i,-j])
      out[i,j] <- pmax(0,x[i] - apply(y,1,max))
    }
    i <- which(x < 0)
    if (length(i) > 0) {
      y <- cbind(0,effects_matrix[i,-j])
      out[i,j] <- pmin(0,x[i] - apply(y,1,min))
    }
  }

  return(out)
}
