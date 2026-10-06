# Modified from singlecelljamboreeR, Peter Carbonetto and Matthew Stephens.
# commit df62f8d169fcac07cd0cb294a3c15051f573ccd5.
# Sources: R/rank_effects.R and R/compute_le_effects.R, MIT.
# Full notice: inst/licenses/singlecelljamboreeR.txt.
# Changes: private names and explicit dimension preservation (one-feature case).
.rank_effects <- function (effects_matrix) {
  if (!(is.matrix(effects_matrix) & is.numeric(effects_matrix)))
    stop("Input \"effects_matrix\" should be a numeric matrix")
  return(matrix(apply(-effects_matrix,2,rank), nrow(effects_matrix), ncol(effects_matrix),
                dimnames = dimnames(effects_matrix)))
}

.compute_le_effects <- function (effects_matrix) {
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
      y <- cbind(0,effects_matrix[i,-j,drop=FALSE])
      out[i,j] <- pmax(0,x[i] - apply(y,1,max))
    }
    i <- which(x < 0)
    if (length(i) > 0) {
      y <- cbind(0,effects_matrix[i,-j,drop=FALSE])
      out[i,j] <- pmin(0,x[i] - apply(y,1,min))
    }
  }

  return(out)
}
