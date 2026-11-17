# FU-1: stability/uncertainty semantics — edge-case regression tests.
# Scoped after auditing test-resampling.R, which already covers: deliberate
# fit_fn failure ledgering, zero-factor refits, checkpoint resume identity,
# serial/parallel identity, and donor occurrences preserving visits. The
# tests below pin only the remaining unaudited edges.

.stab_md <- data.frame(
  sample = paste0("s", 1:12),
  donor = rep(paste0("donor", 1:6), each = 2),
  stringsAsFactors = FALSE
)

# like the existing resampling tests, alignment requires declared units
.stab_fit <- function() {
  fit <- fixture_fit()
  attr(fit, "flashbridge_metadata")$feature_scale <- "intensity"
  fit
}
.stab_cb <- function(X, metadata, sample_side, seed, settings, observation_aux = list()) {
  fit_ebmf(X, sample_side, max_factors = 2L, backfit = TRUE, seed = seed, verbose = 0L, feature_scale = "intensity")
}

test_that("zero-factor REFERENCE is refused with a typed message (FU-1)", {
  X <- fixture_matrix()
  f0 <- fit_ebmf(
    X,
    "rows",
    ebnm_fn = ebnm::ebnm_point_normal,
    max_factors = 0,
    nullcheck = FALSE,
    backfit = FALSE,
    verbose = 0,
    seed = 1
  )
  expect_error(
    factor_stability(X, f0, .stab_md, unit_col = "donor", B = 2, seed = 7),
    "no biological factors available",
    fixed = TRUE
  )
})

test_that("sample-axis reversal is rejected by the ID contract (FU-1)", {
  fit <- .stab_fit()
  expect_error(
    factor_stability(t(fixture_matrix()), fit, .stab_md, unit_col = "donor", B = 2, seed = 7, fit_fn = .stab_cb),
    "sample/feature IDs differ"
  )
})

test_that("preprocessing-stage failures are ledgered with their reason (FU-1)", {
  fit <- .stab_fit()
  calls <- 0L
  failing_pre <- function(Y, metadata, training_ids, seed, settings) {
    calls <<- calls + 1L
    if (calls %% 2L == 0L) {
      stop("injected preprocessing failure")
    }
    list(
      X = Y,
      metadata = metadata,
      manifest = list(scope = "test_passthrough"),
      observation_aux = list()
    )
  }
  res <- factor_stability(
    fixture_matrix(),
    fit,
    .stab_md,
    unit_col = "donor",
    B = 4,
    seed = 7,
    workers = 1,
    fit_fn = .stab_cb,
    preprocess_fn = failing_pre
  )
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
  B_rot <- v$effects %*% t(solve(R))
  expect_equal(tcrossprod(A_rot, B_rot), tcrossprod(v$activity, v$effects), tolerance = 1e-12)
  colnames(A_rot) <- colnames(B_rot) <- c("F1", "F2")
  v_rot <- flashier.utils:::.make_representation(
    A_rot,
    B_rot,
    modifyList(v$manifest, list(fit_id = "rotated"))
  )
  al <- flashier.utils:::.align_programs(v, v_rot, min_congruence = 0.5)
  expect_identical(nrow(al$matches), 2L)
  expect_true(all(al$matches$similarity > 0.5))
})

test_that("n_matched counts refits that assigned the factor, and the quantiles describe only those (review)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("clue")
  X <- fixture_matrix()
  fit <- .stab_fit()
  # every second seed refits one factor only, so the other reference factor is lost in those refits
  losing_cb <- function(X, metadata, sample_side, seed, settings, observation_aux = list()) {
    fit_ebmf(
      X,
      sample_side,
      max_factors = if (seed %% 2L == 0L) 1L else 2L,
      backfit = TRUE,
      seed = seed,
      verbose = 0L,
      feature_scale = "intensity"
    )
  }
  out <- allow_maxiter(factor_stability(X, fit, .stab_md, unit_col = "donor", B = 6, seed = 7, fit_fn = losing_cb))
  finite <- tapply(is.finite(out$matches$similarity), out$matches$reference_factor, sum)
  expect_equal(out$summary$n_matched, as.integer(finite[out$summary$factor]), ignore_attr = TRUE)
  expect_true(all(out$summary$n_matched <= out$summary$n_eligible))
  expect_true(any(out$summary$n_matched < out$summary$n_successful))
})

test_that("a subsample stratum that draws all its units warns and the plan is recorded", {
  skip_if_not_installed("psych")
  skip_if_not_installed("clue")
  X <- fixture_matrix()
  fit <- .stab_fit()
  expect_warning(
    allow_maxiter(factor_stability(
      X,
      fit,
      .stab_md,
      unit_col = "donor",
      B = 1,
      seed = 7,
      fit_fn = .stab_cb,
      scheme = "subsample",
      fraction = 0.9
    )),
    "draws all its units"
  )
  out <- allow_maxiter(factor_stability(
    X,
    fit,
    .stab_md,
    unit_col = "donor",
    B = 1,
    seed = 7,
    fit_fn = .stab_cb,
    scheme = "subsample",
    fraction = 0.5
  ))
  plan <- out$analysis_metadata$resampling_plan
  expect_equal(plan$n_units, 6L)
  expect_equal(plan$n_drawn, 3)
  expect_equal(plan$realised_fraction, 0.5)
})

test_that("best_similarity exposes a near-perfect candidate that one-to-one assignment gave away", {
  skip_if_not_installed("psych")
  skip_if_not_installed("clue")
  view <- function(B) {
    dimnames(B) <- list(paste0("g", 1:3), paste0("F", seq_len(ncol(B))))
    .make_representation(
      matrix(1, 3, ncol(B), dimnames = list(paste0("s", 1:3), colnames(B))),
      B,
      list(
        fit_id = "stub",
        schema_version = "1.0.0",
        sample_side = "rows",
        native_indices = seq_len(ncol(B)),
        scaling = "raw",
        orientation = "as_fit",
        constraints = c(activity = "signed", effects = "signed"),
        feature_scale = "intensity"
      )
    )
  }
  reference <- view(cbind(c(1, 0, 0), c(1, 0.2, 0)))
  target <- view(cbind(c(1, 0.15, 0), c(0, 0, 1)))
  aligned <- .align_programs(reference, target)
  f1 <- aligned$matches[aligned$matches$reference_factor %in% "F1", ]
  # F1 is closest to T1 (0.989) but F2 (0.999) takes T1, leaving F1 with no assigned match
  expect_true(is.na(f1$similarity))
  expect_equal(f1$best_similarity, aligned$similarity_matrix["F1", "F1"])
  expect_gt(f1$best_similarity, 0.98)
})
