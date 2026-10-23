test_that("the public API contains exactly the approved core", {
  core <- c(
    "backproject_contrast", "compare_factor_groups", "factor_activity",
    "factor_effects", "factor_features", "factor_holdout", "factor_intervals",
    "factor_pve", "factor_replicability", "factor_stability", "factor_uncertainty",
    "filter_cells_by_total_membership", "fit_ebmf", "fit_gbcd_programs", "plot_factor_activity", "plot_factor_features",
    "project_factor_activity", "project_factor_features", "provenance",
    "standardize_factors", "summarize_factor_activity", "test_factors",
    "validate_backprojection"
  )
  expect_setequal(getNamespaceExports("flashier.utils"), core)
  records <- provenance(include_dependencies = FALSE)$`function`
  internal <- setdiff(records, core)
  expect_length(internal, 34L)
  expect_true(all(vapply(internal, function(name) {
    is.function(get(name, asNamespace("flashier.utils"), inherits = FALSE))
  }, logical(1))))
  expect_error(getExportedValue("flashier.utils", "fit_multistart_ebmf"),
    "not an exported object")
})

test_that("approved positional calls bind their declared scientific arguments", {
  cases <- list(
    factor_features = quote(factor_features(fit, factor, representation, pool, selection, direction, n)),
    factor_holdout = quote(factor_holdout(X, metadata, unit, folds, callback, seed)),
    factor_intervals = quote(factor_intervals(fit, side, orientation, representation)),
    factor_stability = quote(factor_stability(X, reference, metadata, unit, seed)),
    factor_uncertainty = quote(factor_uncertainty(fit, side, orientation, representation)),
    plot_factor_activity = quote(plot_factor_activity(fit, metadata, group, factors, style)),
    plot_factor_features = quote(plot_factor_features(fit, factors, representation, selection)),
    test_factors = quote(test_factors(fit, metadata, formula, engine, unit, independence, id)),
    validate_backprojection = quote(validate_backprojection(reconstructed, observed, scale, contrast))
  )
  expected <- list(
    factor_features = list(representation = quote(representation), compare_factors = quote(pool),
      select = quote(selection), direction = quote(direction), n = quote(n)),
    factor_holdout = list(seed = quote(seed)),
    factor_intervals = list(sample_side = quote(orientation), representation = quote(representation)),
    factor_stability = list(seed = quote(seed)),
    factor_uncertainty = list(sample_side = quote(orientation), representation = quote(representation)),
    plot_factor_activity = list(factors = quote(factors), style = quote(style)),
    plot_factor_features = list(representation = quote(representation), select = quote(selection)),
    test_factors = list(unit_col = quote(unit), independence = quote(independence), id_col = quote(id)),
    validate_backprojection = list(feature_scale = quote(scale), contrast = quote(contrast))
  )
  for (name in names(cases)) {
    matched <- as.list(match.call(getExportedValue("flashier.utils", name), cases[[name]]))
    expect_identical(matched[names(expected[[name]])], expected[[name]], info = name)
  }
})

test_that("feature selection uses the new positional order with unchanged values", {
  view <- fixture_view()
  positional <- factor_features(view, "F1", NULL, c("F1", "F2"), "largest", "absolute", 2L)
  named <- factor_features(view, factor = "F1", representation = NULL,
    compare_factors = c("F1", "F2"), select = "largest", direction = "absolute", n = 2L)
  expect_identical(lapply(positional, identity), lapply(named, identity))
  expect_identical(positional$feature, c("g1", "g2"))
})
