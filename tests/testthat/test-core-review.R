test_that("bare native fits declare units without private metadata", {
  X <- fixture_matrix()
  for (side in c("rows", "columns")) {
    Y <- if (side == "rows") X else t(X)
    set.seed(12)
    fit <- flashier::flash(Y, greedy_Kmax = 2, nullcheck = FALSE,
      ebnm_fn = ebnm::ebnm_normal, backfit = TRUE, verbose = 0)
    expect_null(attr(fit, "flashbridge_metadata"))
    declaration <- list(scope = "identity", transformation = "none")
    v <- standardize_factors(fit, sample_side = side, scaling = "raw",
      feature_scale = "intensity", preprocessing = declaration,
      independent_unit = "sample", prior_support = c(activity = "signed", effects = "signed"))
    expect_identical(v$manifest$feature_scale, "intensity")
    expect_identical(v$manifest$preprocessing, declaration)
    expected <- stats::fitted(fit)
    if (side == "columns") expected <- t(expected)
    expect_equal(unname(tcrossprod(v$activity, v$effects)), unname(expected), tolerance = 1e-10)
    expect_equal(sum(factor_replicability(v, v, min_congruence = .99)$matches$matched), 2)
    expect_true(all(is.finite(factor_uncertainty(fit, sample_side = side,
      representation = v)$posterior.sd)))
    projected <- project_factor_activity(v, X, sample_side = "rows", feature_scale = "intensity",
      preprocessing = declaration, ebnm_fn = ebnm::ebnm_normal, maxiter = 1000)
    expect_equal(projected$representation$effects, v$effects)
    reloaded <- unserialize(serialize(fit, NULL))
    expect_equal(standardize_factors(reloaded, sample_side = side)$activity,
      standardize_factors(fit, sample_side = side)$activity)
    for (axis in c("activity", "effects")) {
      which_dim <- if ((axis == "activity") == (side == "rows")) "loadings" else "factors"
      native <- flashier::flash_factors_fix(fit, 1L, which_dim = which_dim)
      fixed <- fix_factors(fit, "F1", axis, sample_side = side)
      expect_equal(fixed$L_pm, native$L_pm)
      expect_equal(fixed$F_pm, native$F_pm)
      expect_equal(flashier::flash_fit(fixed), flashier::flash_fit(native))
      expect_identical(attr(fixed, "flashbridge_metadata")$sample_side, side)
    }
  }
  recorded <- fixture_fit()
  expect_error(standardize_factors(recorded, feature_scale = ""), "units")
  recorded <- fit_ebmf(X, "rows", max_factors = 1, verbose = 0, feature_scale = "intensity")
  expect_error(standardize_factors(recorded, feature_scale = "log_counts"), "conflicts")
  expect_error(standardize_factors(recorded, preprocessing = list(1)), "named list")
})

test_that("zero-factor references and targets compare without fabricated matches", {
  v <- fixture_view()
  z <- .make_representation(v$activity[, FALSE, drop = FALSE], v$effects[, FALSE, drop = FALSE],
    modifyList(v$manifest, list(native_indices = integer())))
  expect_equal(factor_replicability(z, v)$unmatched_target, c("F1", "F2"))
  expect_equal(nrow(factor_replicability(z, z)$matches), 0L)
  expect_error(factor_features_global(v, character()), "at least one")
})

test_that("draw selectors and declared feature units cannot be silently ignored", {
  v <- fixture_view()
  draws <- expand.grid(draw = 1:4, factor = c("F1", "F2"))
  draws$estimate <- seq_len(nrow(draws)) / 10
  for (selection in c("term", "contrasts", "conditions")) {
    args <- list(fit = v, factor_effects = draws, input = "draws", factor_basis = v$manifest)
    args[[selection]] <- "absent"
    expect_error(do.call(backproject_contrast, args), "Draw table has no")
  }
  observed <- data.frame(feature = c("g1", "g2"), estimate = c(1, 2))
  for (scale in list(NULL, NA_character_, "", "unknown", 1, c("a", "b"))) {
    expect_error(validate_backprojection(observed, observed, feature_scale = scale), "units")
  }
})

test_that("group labels cannot overwrite reserved summary columns", {
  v <- fixture_view()
  for (column in c("factor", "basis_id", "mean", "n", "centered")) {
    md <- data.frame(sample = rownames(v$activity))
    md[[column]] <- c("a", "a", "b")
    expect_error(summarize_factor_activity(v, md, column), "reserved output")
    expect_error(center_factor_activity(v, md, column), "reserved output")
  }
})

test_that("non-lm engines refuse unsupported formula offsets", {
  f <- association_fixture()
  f$metadata$z <- 5 * (f$metadata$group == "B")
  expect_error(test_factors(f$view, f$metadata, ~group + offset(z), "limma",
    unit_col = "donor", independence = "independent"), "offsets.*only")
})

test_that("holdout validates scoring axes as well as inference axes", {
  X <- fixture_matrix()
  md <- data.frame(sample = rownames(X), donor = rep(paste0("d", 1:6), each = 2))
  folds <- data.frame(sample = md$sample, fold = c(rep("train", 8), rep("test", 4)))
  prep <- list(fit = function(X, metadata, training_ids, seed, settings) {
    list(X = X, metadata = metadata, manifest = list(scope = "identity"), state = NULL)
  }, transform = function(X, metadata, state, seed, settings) {
    list(X = if (anyNA(X)) X else X[, rev(seq_len(ncol(X))), drop = FALSE], metadata = metadata)
  })
  callback <- function(X, metadata, sample_side, seed, settings, observation_aux = list()) {
    flashier::flash(X, greedy_Kmax = 2, ebnm_fn = ebnm::ebnm_normal,
      nullcheck = FALSE, backfit = TRUE, verbose = 0)
  }
  out <- factor_holdout(X, md, "donor", folds, callback, preprocess_fn = prep,
    evaluation = "masked", activity_features = colnames(X)[1:6],
    scoring_features = colnames(X)[7:8], seed = 4,
    control = list(sample_side = "rows", evaluated_folds = "test",
      projection_ebnm_fn = ebnm::ebnm_normal))
  expect_match(out$failures$reason, "Scoring transform changed")
  expect_null(out$metrics)
})

test_that("descriptive ranking accepts literal nonsyntactic metadata names", {
  v <- fixture_view()
  md <- data.frame(sample = rownames(v$activity), age = c(1, 3, 2))
  ordinary <- rank_factors_by_metadata(v, md, "age")
  names(md)[2] <- "age years"
  literal <- rank_factors_by_metadata(v, md, "age years")
  expect_equal(literal$value, ordinary$value)
})

test_that("stability accepts declared views from native-only refit callbacks", {
  X <- fixture_matrix()
  md <- data.frame(sample = rownames(X), donor = rep(paste0("d", 1:6), each = 2))
  callback <- function(X, metadata, sample_side, seed, settings, observation_aux = list()) {
    set.seed(seed)
    native <- flashier::flash(X, greedy_Kmax = 2, ebnm_fn = ebnm::ebnm_normal,
      nullcheck = FALSE, backfit = TRUE, verbose = 0)
    standardize_factors(native, sample_side = sample_side, feature_scale = "intensity",
      prior_support = c(activity = "signed", effects = "signed"))
  }
  for (side in c("rows", "columns")) {
    Y <- if (side == "rows") X else t(X)
    reference <- callback(Y, md, side, 12, NULL)
    out <- factor_stability(Y, reference, md, unit_col = "donor", B = 2,
      fit_fn = callback, seed = 5, min_congruence = .2)
    expect_equal(nrow(out$failures), 0L)
    expect_true(all(out$replicates$status == "supported"))
    expect_equal(nrow(out$replicates), 2L)
  }
})
