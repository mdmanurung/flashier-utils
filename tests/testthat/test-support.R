test_that("ordinary and serialized legacy views keep their original manifests", {
  fit <- fixture_fit()
  view <- standardize_factors(fit)
  expect_null(view$manifest$factor_support)
  expect_identical(view$manifest$schema_version, "1.0.0")
  expected <- data.frame(factor = c("F1", "F2"), native_index = 1:2,
    role = "factor", activity_support = "signed", effect_support = "signed")
  expect_identical(.factor_support(view$manifest), expected)
  path <- tempfile()
  on.exit(unlink(path))
  saveRDS(view, path)
  expect_identical(.validate_view(readRDS(path)), view)
  expect_identical(.make_representation(view$activity, view$effects, view$manifest)$manifest,
    view$manifest)
})

test_that("per-factor support and roles survive native edits and serialization", {
  X <- fixture_matrix()
  fit <- fit_ebmf(X, "rows", intercept = "feature", max_factors = 2,
    ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace),
    backfit = TRUE, nullcheck = FALSE, seed = 12, verbose = 0)
  view <- standardize_factors(fit)
  expected <- data.frame(factor = c("F1", "F2"), native_index = 2:3,
    role = "factor", activity_support = "signed", effect_support = "signed")
  expect_identical(view$manifest$factor_support, expected)
  expect_identical(view$manifest$schema_version, "1.1.0")
  reordered <- reorder_factors(fit, c("F2", "F1"))
  records <- standardize_factors(reordered)$manifest$factor_support
  expect_identical(records$factor, c("F2", "F1"))
  expect_identical(records$native_index, 2:3)
  remaining <- remove_factors(reordered, "F1")
  records <- standardize_factors(remaining)$manifest$factor_support
  expect_identical(records$factor, "F2")
  expect_identical(records$native_index, 2L)
  path <- tempfile()
  on.exit(unlink(path))
  saveRDS(remaining, path)
  expect_identical(standardize_factors(readRDS(path))$manifest, standardize_factors(remaining)$manifest)
  raw <- standardize_factors(fixture_fit(), scaling = "raw")
  initialized <- fit_initialized_ebmf(X, "rows", raw$activity, raw$effects,
    activity_prior = ebnm::ebnm_normal, effect_prior = ebnm::ebnm_point_laplace,
    fixed_side = "effects", verbose = 0, seed = 12)
  expect_identical(standardize_factors(initialized)$manifest$factor_support$role,
    rep("initialized", 2L))
})

test_that("support records are hashed, validated and cannot be overridden", {
  view <- fixture_view()
  records <- .factor_support(view$manifest)
  records$activity_support[2] <- "nonnegative"
  manifest <- view$manifest
  manifest$support_version <- "1.0.0"
  manifest$schema_version <- manifest$representation_version <- "1.1.0"
  manifest$factor_support <- records
  manifest$constraints["activity"] <- "unknown"
  mixed <- .make_representation(view$activity, view$effects, manifest)
  expect_identical(mixed$factors$constraints, c("signed/signed", "nonnegative/signed"))
  expect_false(identical(mixed$manifest$basis_id, view$manifest$basis_id))
  changed <- mixed
  changed$manifest$factor_support$activity_support[1] <- "unknown"
  expect_error(.validate_view(changed), "modified")
  for (key in c("factor", "native_index", "role", "activity_support", "effect_support")) {
    broken <- manifest
    broken$factor_support[[key]][1] <- switch(key,
      factor = "missing", native_index = 8L, role = "offset", "invalid")
    expect_error(.make_representation(view$activity, view$effects, broken), "factor_support")
  }
  broken <- manifest
  broken$support_version <- "9.0.0"
  expect_error(.make_representation(view$activity, view$effects, broken), "factor_support")
  broken <- manifest
  broken$constraints["activity"] <- "signed"
  expect_error(.make_representation(view$activity, view$effects, broken), "conservatively")
  fit <- fit_ebmf(fixture_matrix(), "rows", max_factors = 1,
    intercept = "feature", verbose = 0)
  expect_error(standardize_factors(fit, prior_support = c(activity = "nonnegative", effects = "signed")), "override")
  raw <- standardize_factors(fit, scaling = "raw")
  forged <- raw$manifest
  forged$factor_support$activity_support[] <- "nonnegative"
  forged$constraints["activity"] <- "nonnegative"
  forged <- .make_representation(raw$activity, raw$effects, forged, raw$offset)
  expect_error(factor_uncertainty(fit, "activity", representation = forged), "support/roles")
})

test_that("projection keeps fixed-side support and records the estimated prior", {
  view <- fixture_view()
  manifest <- view$manifest
  records <- .factor_support(manifest)
  records$effect_support[2] <- "unknown"
  manifest$constraints["effects"] <- "unknown"
  manifest$factor_support <- records
  manifest$support_version <- "1.0.0"
  source <- .make_representation(view$activity, view$effects, manifest)
  query <- tcrossprod(source$activity, source$effects)
  projected <- suppressWarnings(project_factor_activity(source, query, "rows",
    "arcsinh_intensity", list(scope = "identity"), factors = "F2",
    ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)))
  out <- projected$representation$manifest$factor_support
  expect_identical(out$factor, "F2")
  expect_identical(out$native_index, 1L)
  expect_identical(out$activity_support, "signed")
  expect_identical(out$effect_support, "unknown")
  expect_identical(standardize_factors(projected$fit)$manifest$factor_support, out)
  expect_error(standardize_factors(projected$fit, orientation = "anchor_feature"), "signed prior")
  transposed <- suppressWarnings(project_factor_features(source, t(query), "columns",
    "arcsinh_intensity", factors = "F1",
    ebnm_fn = list(ebnm::ebnm_point_laplace, ebnm::ebnm_normal)))
  expect_identical(transposed$representation$manifest$factor_support$activity_support, "signed")
  expect_identical(transposed$representation$manifest$factor_support$effect_support, "signed")
})
