# Fresh-native verification oracles for FND-001/FND-002 (2026-10-16).
#
# The original delivery receipts for the five defect-closure commits
# (91c9fc8..0e0bdc0) were session-local and are lost (see HANDOFF reality
# check). These oracles re-establish the closures from public flashier calls
# only: every expected value is rebuilt from scratch beside the package result,
# so inherited source metadata or positional/native-index confusion cannot pass.

test_that("FND-001 oracle: activity projection equals a fresh native pipeline", {
  for (orientation in c("rows", "columns")) {
    fit <- fixture_fit(orientation)
    view <- standardize_factors(fit, scaling = "ldf", d_location = "activity")
    # fixture_fit deliberately leaves units unknown; declare them explicitly.
    view <- .make_representation(view$activity, view$effects,
      utils::modifyList(view$manifest, list(feature_scale = "simulated")))
    # A query that is NOT the source data, so a copied source PVE could never
    # masquerade as the projection fit's own PVE.
    Y <- 2 * fixture_matrix()
    query <- if (orientation == "rows") Y else t(Y)
    priors <- list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)
    result <- suppressWarnings(project_factor_activity(view, query, orientation,
      "simulated", list(scope = "identity"), ebnm_fn = priors, maxiter = 200, var_type = 0L))
    # Fresh pipeline: public flashier calls only, source effects fixed exactly.
    genes <- rownames(view$effects)
    Ysem <- if (orientation == "rows") query[, genes] else t(query)[, genes]
    B <- view$effects[genes, , drop = FALSE]
    A0 <- matrix(0.1, nrow(Ysem), ncol(B), dimnames = list(rownames(Ysem), colnames(B)))
    native_data <- if (orientation == "rows") Ysem else t(Ysem)
    object <- flashier::flash_init(native_data, S = NULL, var_type = 0L)
    pair <- if (orientation == "rows") list(A0, B) else list(B, A0)
    object <- flashier::flash_factors_init(object, pair, ebnm_fn = priors)
    object <- flashier::flash_factors_fix(object, seq_len(ncol(B)),
      which_dim = if (orientation == "rows") "factors" else "loadings")
    object <- flashier::flash_backfit(object, maxiter = 200, extrapolate = FALSE, verbose = 0L)
    # Mean oracle: posterior means come from the projection fit, not the source.
    expect_equal(unname(result$fit$L_pm), unname(object$L_pm), tolerance = 1e-12)
    expect_equal(unname(result$fit$F_pm), unname(object$F_pm), tolerance = 1e-12)
    # SD oracle: native posterior SDs of the projection fit, both sides.
    expect_equal(unname(result$fit$L_psd), unname(object$L_psd), tolerance = 1e-12)
    expect_equal(unname(result$fit$F_psd), unname(object$F_psd), tolerance = 1e-12)
    # PVE oracle: the projection fit's own native PVE, never the source's.
    expect_equal(result$fit$pve, object$pve, tolerance = 1e-12)
    expect_false(isTRUE(all.equal(result$fit$pve, view$manifest$native_pve)))
    # Draw oracle: posterior draws reproduce under the fresh pipeline and seed.
    set.seed(39); d1 <- result$fit$sampler(20)
    set.seed(39); d2 <- object$sampler(20)
    expect_equal(d1, d2, tolerance = 1e-12)
    # Manifest resets and preserved source provenance.
    projected <- result$representation
    expect_identical(projected$manifest$scaling, "raw")
    expect_null(projected$manifest$ldf_type)
    expect_null(projected$manifest$d_location)
    expect_equal(projected$manifest$activity_multiplier, rep(1, ncol(B)))
    expect_equal(projected$manifest$effects_multiplier, rep(1, ncol(B)))
    expect_identical(projected$manifest$orientation, "as_fit")
    expect_equal(projected$manifest$native_pve, object$pve[projected$manifest$native_indices],
      tolerance = 1e-12)
    expect_equal(projected$manifest$source_coordinates$native_pve, view$manifest$native_pve,
      tolerance = 1e-12)
    expect_identical(projected$manifest$source_coordinates$scaling, "ldf")
    expect_equal(projected$manifest$source_basis_id, view$manifest$basis_id)
    # Semantic reconstruction identity on the shared genes.
    expect_true(result$reconstruction[["fit_r2"]] > 0.9)
  }
})

test_that("FND-001 oracle: feature projection moments and draws equal a fresh native pipeline", {
  M <- cbind(P1 = c(s1 = 1, s2 = 2, s3 = 5), P2 = c(3, 0, -2))
  priors <- list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)
  for (orientation in c("rows", "columns")) {
    fit <- fixture_fit(orientation)
    view <- standardize_factors(fit, scaling = "raw")
    p <- if (orientation == "rows") priors else rev(priors)
    result <- suppressWarnings(project_factor_features(view,
      if (orientation == "rows") M else t(M), orientation, "protein", ebnm_fn = p, maxiter = 200))
    a <- view$activity[rownames(M), ]
    b <- t(qr.coef(qr(a), as.matrix(M)))
    object <- flashier::flash_init(if (orientation == "rows") M else t(M),
      var_type = if (orientation == "rows") 2L else 1L)
    object <- flashier::flash_factors_init(object,
      if (orientation == "rows") list(a, b) else list(b, a), ebnm_fn = p)
    object <- flashier::flash_factors_fix(object, seq_len(ncol(a)),
      which_dim = if (orientation == "rows") "loadings" else "factors")
    object <- flashier::flash_backfit(object, maxiter = 200, extrapolate = FALSE, verbose = 0L)
    expect_equal(unname(result$fit$L_pm), unname(object$L_pm), tolerance = 1e-12)
    expect_equal(unname(result$fit$F_pm), unname(object$F_pm), tolerance = 1e-12)
    expect_equal(unname(result$fit$L_psd), unname(object$L_psd), tolerance = 1e-12)
    expect_equal(unname(result$fit$F_psd), unname(object$F_psd), tolerance = 1e-12)
    expect_equal(result$fit$pve, object$pve, tolerance = 1e-12)
    set.seed(39); d1 <- result$fit$sampler(20)
    set.seed(39); d2 <- object$sampler(20)
    expect_equal(d1, d2, tolerance = 1e-12)
    expect_equal(result$representation$activity, a, tolerance = 1e-12)
  }
})

test_that("FND-002 oracle: partial fixing pins named entries through backfit in both orientations", {
  for (side in c("rows", "columns")) {
    f <- fixture_fit(side)
    view <- standardize_factors(f, scaling = "raw")
    for (axis in c("activity", "effects")) {
      axis_ids <- rownames(view[[axis]])
      ids <- axis_ids[c(2, 5)]
      fixed <- fix_factors(f, "F1", axis, ids = ids)
      refit <- refit_factors(fixed, control = list(maxiter = 30L, verbose = 0L))
      mode <- if ((axis == "activity") == (side == "rows")) "L" else "F"
      # The named entries keep their exact pinned values through the refit.
      expect_equal(refit[[paste0(mode, "_pm")]][match(ids, axis_ids), 1],
        f[[paste0(mode, "_pm")]][match(ids, axis_ids), 1], tolerance = 1e-12)
      # The refit equals a fresh native mask-fix plus backfit of every factor.
      native <- flashier::flash_factors_fix(f, 1,
        which_dim = if ((axis == "activity") == (side == "rows")) "loadings" else "factors",
        fixed_idx = rownames(view[[axis]]) %in% ids)
      native <- flashier::flash_backfit(native, kset = seq_len(f$n_factors),
        maxiter = 30L, verbose = 0L)
      expect_equal(unname(refit[[paste0(mode, "_pm")]]),
        unname(native[[paste0(mode, "_pm")]]), tolerance = 1e-12)
    }
  }
})

# October 8 raw-helper counterexample, reconstructed as a public-API regression:
# the lost benchmark seed-consensus helper was input-order dependent because it
# did not sort stable member IDs. The public multistart primitive must be
# invariant to the input order of init_seeds: candidates are matched by seed,
# and the selected solution may not depend on start numbering. Exact-ELBO ties
# resolve to the smallest start ID, which is order-dependent numbering of the
# same supplied seed set; identical solutions make that relabeling harmless.
test_that("multistart selection is invariant to init_seeds input order", {
  X <- fixture_matrix()
  seeds <- c(101, 202, 303)
  args <- list(X = X, sample_side = "rows", ebnm_fn = ebnm::ebnm_point_laplace,
    n_starts = 3L, max_factors = 2L, backfit = TRUE, nullcheck = FALSE,
    verbose = 0L, feature_scale = "simulated_centered_intensity")
  r1 <- do.call(fit_multistart_ebmf, c(args, list(init_seeds = seeds)))
  r2 <- do.call(fit_multistart_ebmf, c(args, list(init_seeds = rev(seeds))))
  expect_setequal(r1$candidates$init_seed, r2$candidates$init_seed)
  m <- match(r2$candidates$init_seed, r1$candidates$init_seed)
  expect_equal(r2$candidates$elbo, r1$candidates$elbo[m])
  expect_equal(r2$candidates$candidate_digest, r1$candidates$candidate_digest[m])
  expect_equal(r2$candidates$rank, r1$candidates$rank[m])
  expect_equal(sum(r1$candidates$selected), 1L)
  expect_equal(r1$candidates$elbo[r1$candidates$selected],
    r2$candidates$elbo[r2$candidates$selected])
  expect_equal(r1$best_fit$L_pm, r2$best_fit$L_pm, tolerance = 1e-12)
  expect_equal(r1$best_fit$F_pm, r2$best_fit$F_pm, tolerance = 1e-12)
})
