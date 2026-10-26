gom_fixture <- function() {
  set.seed(1); n <- 40; m <- 20
  L <- matrix(rgamma(n * 2, 1), n, 2); L <- L / rowSums(L)
  dimnames(L) <- list(paste0("s", 1:n), c("T1", "T2"))
  F <- matrix(rgamma(m * 2, 1), m, 2)
  X <- matrix(rpois(n * m, 20 * tcrossprod(L, F)), n, m, dimnames = list(rownames(L), paste0("g", 1:m)))
  list(X = X, L = L)
}

test_that("GoM DE matches native fastTopics under the same seed", {
  skip_if_not_installed("fastTopics")
  d <- gom_fixture()
  for (shrink in c("ash", "none")) {
    out <- test_gom_de(d$X, d$L, "rows", seed = 3, shrink_method = shrink, control = list(ns = 100))
    set.seed(3)
    native <- suppressWarnings(fastTopics::de_analysis(d$L, d$X, shrink.method = shrink,
      control = list(ns = 100, nc = 1), verbose = FALSE))
    expect_equal(out$statistics$log2_fc, as.vector(native$postmean))
    expect_equal(out$statistics$lower, as.vector(native$lower))
    expect_identical(out$statistics$feature[1:20], rownames(native$postmean))
    expect_identical(unique(out$statistics$factor), c("T1", "T2"))
    if (shrink == "ash") {
      expect_equal(out$statistics$lfsr, as.vector(native$lfsr))
      expect_true(all(is.na(out$statistics$neg_log10_p)))
    } else {
      expect_equal(out$statistics$neg_log10_p, as.vector(native$lpval))
      expect_true(all(is.na(out$statistics$lfsr)))
    }
  }
  sparse <- test_gom_de(Matrix::Matrix(t(d$X) * 1, sparse = TRUE), d$L[40:1, ], "columns", seed = 3,
    control = list(ns = 100))
  expect_equal(sparse$statistics$log2_fc, out <- test_gom_de(d$X, d$L, "rows", seed = 3,
    control = list(ns = 100))$statistics$log2_fc)
})

test_that("GoM DE restores RNG state and records its settings", {
  skip_if_not_installed("fastTopics")
  d <- gom_fixture()
  set.seed(9); before <- .Random.seed
  out <- test_gom_de(d$X, d$L, "rows", seed = 3, control = list(ns = 50))
  expect_identical(.Random.seed, before)
  meta <- out$analysis_metadata
  expect_identical(meta$inference, "conditional_on_fixed_memberships")
  expect_identical(meta$control$conf.level, 0.68)
  expect_identical(meta$control$nc, 1)
  expect_type(meta$warnings, "character")
  expect_null(out$native)
})

test_that("GoM DE refuses invalid counts, memberships and settings", {
  d <- gom_fixture()
  run <- function(X = d$X, L = d$L, ...) test_gom_de(X, L, "rows", seed = 1, ...)
  bad <- d$X; bad[1, 1] <- 0.5; expect_error(run(bad), "nonnegative integers")
  bad <- d$X; bad[1, 1] <- -1; expect_error(run(bad), "nonnegative integers")
  bad <- d$X; bad[1, ] <- 0; expect_error(run(bad), "positive library")
  for (scale in c(0.9, 1.1)) { bad <- d$L; bad[1, ] <- bad[1, ] * scale; expect_error(run(L = bad), "sum to one") }
  expect_error(run(L = d$L[, 1, drop = FALSE]), "at least two")
  bad <- d$L; rownames(bad)[1] <- "other"; expect_error(run(L = bad), "same sample IDs")
  bad <- d$L; bad[1, 1] <- -0.1; bad[1, 2] <- 1.1; expect_error(run(L = bad), "nonnegative")
  expect_error(test_gom_de(d$X, d$L, "rows"), "seed")
  skip_if_not_installed("fastTopics")
  expect_error(run(control = list(unknown = 1)), "Unknown fastTopics")
})

test_that("missing fastTopics fails before engine execution", {
  original <- base::requireNamespace
  testthat::local_mocked_bindings(requireNamespace = function(package, ...) {
    if (identical(package, "fastTopics")) return(FALSE)
    original(package, ...)
  }, .package = "base")
  d <- gom_fixture()
  expect_error(test_gom_de(d$X, d$L, "rows", seed = 1), "unavailable")
})
