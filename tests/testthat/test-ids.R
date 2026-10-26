test_that("IDs and missingness are validated without changing sparse zeros", {
  X <- matrix(c(0,1,2,0),2,dimnames=list(c("s1","s2"),c("g1","g2")))
  sparse <- methods::as(Matrix::Matrix(X,sparse=TRUE),"dgCMatrix")
  expect_identical(.validate_matrix(sparse), sparse)
  expect_error(.validate_matrix(unname(X)), "IDs")
  expect_equal(rownames(.validate_matrix(unname(X),"generate")), c("row1","row2"))
  bad <- X; rownames(bad) <- c("s1","s1")
  expect_error(.validate_matrix(bad), "unique")
  bad <- X; bad[1] <- NaN
  expect_error(.validate_matrix(bad), "NaN")
  metadata <- data.frame(sample=c("s2","s1"), donor=c("d2","d1"))
  joined <- .match_ids(rownames(X),metadata)
  expect_equal(joined$metadata$donor,c("d1","d2"))
  expect_error(.match_ids(rownames(X),metadata[1,,drop=FALSE]),"Missing metadata")
  extra <- rbind(metadata,data.frame(sample="s3",donor="d3"))
  expect_error(.match_ids(rownames(X),extra,extra="fail"),"Extra")
  expect_equal(.match_ids(rownames(X),extra)$exclusions$sample,"s3")
  metadata$donor[1] <- NA
  expect_error(.match_ids(rownames(X),metadata),"Missing metadata values")
  expect_equal(.match_ids(rownames(X),metadata,na_action="omit")$matching$sample,"s1")
})

test_that("selections take character or factor IDs, never positions", {
  v <- fixture_view()
  expect_equal(factor_activity(v, samples = factor(c("s2", "s1"))), factor_activity(v, samples = c("s2", "s1")), ignore_attr = TRUE)
  expect_error(factor_activity(v, factors = 2), "character IDs, not positions")
  expect_error(factor_activity(v, factors = c(TRUE, FALSE)), "character IDs, not positions")
})
