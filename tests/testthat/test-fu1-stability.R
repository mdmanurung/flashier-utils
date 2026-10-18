# FU-1: stability/uncertainty semantics — edge-case regression tests.
# Scoped after auditing test-resampling.R, which already covers: deliberate
# fit_fn failure ledgering, zero-factor refits, checkpoint resume identity,
# serial/parallel identity, and donor occurrences preserving visits. The
# tests below pin only the remaining unaudited edges.

.stab_md <- data.frame(sample = paste0("s", 1:12),
                       donor = rep(paste0("donor", 1:6), each = 2),
                       stringsAsFactors = FALSE)

# like the existing resampling tests, alignment requires declared units
.stab_fit <- function() {
  fit <- fixture_fit()
  attr(fit, "flashbridge_metadata")$feature_scale <- "intensity"
  fit
}
.stab_cb <- function(X, metadata, sample_side, seed, settings, observation_aux = list()) {
  fit_ebmf(X, sample_side, max_factors = 2L, backfit = TRUE, seed = seed,
           verbose = 0L, feature_scale = "intensity")
}

test_that("zero-factor REFERENCE is refused with a typed message (FU-1)", {
  X <- fixture_matrix()
  f0 <- fit_ebmf(X, "rows", ebnm_fn = ebnm::ebnm_point_normal, max_factors = 0,
                 nullcheck = FALSE, backfit = FALSE, verbose = 0, seed = 1)
  expect_error(
    factor_stability(X, f0, .stab_md, unit_col = "donor", B = 2, seed = 7),
    "factor|biological|zero|no ",
    ignore.case = TRUE
  )
})

test_that("sample-axis reversal is rejected by the ID contract (FU-1)", {
  fit <- .stab_fit()
  expect_error(
    factor_stability(t(fixture_matrix()), fit, .stab_md, unit_col = "donor",
                     B = 2, seed = 7, fit_fn = .stab_cb),
    "sample/feature IDs differ"
  )
})

test_that("preprocessing-stage failures are ledgered with their reason (FU-1)", {
  fit <- .stab_fit()
  calls <- 0L
  failing_pre <- function(Y, metadata, training_ids, seed, settings) {
    calls <<- calls + 1L
    if (calls %% 2L == 0L) stop("injected preprocessing failure")
    list(X = Y, metadata = metadata, manifest = list(scope = "test_passthrough"),
         observation_aux = list())
  }
  res <- factor_stability(fixture_matrix(), fit, .stab_md, unit_col = "donor",
                          B = 4, seed = 7, workers = 1,
                          fit_fn = .stab_cb, preprocess_fn = failing_pre)
  expect_identical(nrow(res$failures), 2L)
  expect_true(all(grepl("injected preprocessing failure", res$failures$reason)))
  expect_identical(unique(res$summary$n_requested), 4)
  expect_identical(unique(as.integer(res$summary$n_failed)), 2L)
  expect_identical(unique(as.integer(res$summary$n_successful)), 2L)
})

test_that("rotated rank-2 subspace is recoverable via existing alignment (FU-1/U07)", {
  ## A reconstruction-preserving rotation of a 2-factor view spans the same
  ## subspace; congruence alignment should match both rotated factors back.
  v <- fixture_view()
  theta <- 0.4
  R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
  A_rot <- v$activity %*% R
  B_rot <- v$effects %*% solve(R)
  colnames(A_rot) <- colnames(B_rot) <- c("F1", "F2")
  v_rot <- flashier.utils:::.make_representation(
    A_rot, B_rot, modifyList(v$manifest, list(fit_id = "rotated")))
  al <- flashier.utils:::.align_programs(v, v_rot, min_congruence = 0.5)
  expect_identical(nrow(al$matches), 2L)
  expect_true(all(al$matches$similarity > 0.5))
})
