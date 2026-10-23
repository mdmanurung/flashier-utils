membership_view <- function() {
  A <- rbind(broad = c(10, 10, 0), next_cell = c(4, 4, .2),
             a = c(5, 0, 0), b = c(0, 5, 0), c = c(0, 0, 1))
  colnames(A) <- paste0("F", 1:3)
  B <- matrix(c(1, -2, 3, 4, -1, 2), 2,
              dimnames = list(c("g1", "g2"), colnames(A)))
  .make_representation(A, B, list(fit_id = "membership_fixture",
    schema_version = "1.0.0", sample_side = "rows", native_indices = 1:3,
    scaling = "raw", orientation = "as_fit",
    constraints = c(activity = "nonnegative", effects = "signed"),
    feature_scale = "log_expression"), c(g1 = 2, g2 = -3))
}

test_that("ImmGenT filtering repeats normalization and preserves paired reconstruction", {
  view <- membership_view()
  before <- serialize(view, NULL)
  out <- filter_cells_by_total_membership(view, max_val = 1.5)
  expect_identical(out$kept_samples, c("a", "b", "c"))
  expect_identical(out$removed_samples, c("broad", "next_cell"))
  expect_equal(out$history$removed, c(1, 1, rep(0, 10)))
  expect_equal(unname(out$activity), diag(3))
  expect_equal(tcrossprod(out$activity, out$effects),
    tcrossprod(view$activity[out$kept_samples, , drop = FALSE], view$effects))
  expect_identical(out$offset, view$offset)
  expect_true(out$converged)
  expect_identical(serialize(view, NULL), before)
  expect_false(attr(out, "analysis_metadata")$native_fit_refitted)
  expect_equal(provenance("filter_cells_by_total_membership")$code_reuse, "none")
  expect_warning(partial <- filter_cells_by_total_membership(view,
    max_val = 1.5, numiter = 1), "reached numiter")
  expect_false(partial$converged)
  expect_identical(partial$kept_samples, c("next_cell", "a", "b", "c"))
})

test_that("filtering handles inclusive thresholds, zero factors and singleton axes", {
  view <- fixture_view_zero_activity()
  out <- filter_cells_by_total_membership(view, max_val = 1)
  expect_identical(out$kept_samples, rownames(view$activity))
  expect_true(all(out$activity[, "F2"] == 0))
  expect_identical(out$effects[, "F2"], view$effects[, "F2"])
  expect_identical(attr(out, "analysis_metadata")$zero_max_factors, "F2")
  single <- .make_representation(matrix(2, 1, 1, dimnames = list("s", "F1")),
    matrix(-3, 1, 1, dimnames = list("g", "F1")),
    list(fit_id = "single", schema_version = "1.0.0", sample_side = "rows",
      native_indices = 1L, scaling = "raw", orientation = "as_fit",
      constraints = c(activity = "nonnegative", effects = "signed")))
  one <- filter_cells_by_total_membership(single, max_val = 1)
  expect_identical(dim(one$activity), c(1L, 1L))
  expect_identical(dim(one$effects), c(1L, 1L))
  expect_equal(one$activity[1, 1], 1)
  expect_equal(one$effects[1, 1], -6)
  expect_length(one$removed_samples, 0)
  expect_error(filter_cells_by_total_membership(single, max_val = 0), "all samples")
})

test_that("filtering validates activities, orientation and controls", {
  expect_error(filter_cells_by_total_membership(fixture_view()), "nonnegative")
  expect_error(filter_cells_by_total_membership(membership_view(), "columns"), "conflicts")
  for (bad in list(-1, NA_real_, Inf, "10", c(1, 2))) {
    expect_error(filter_cells_by_total_membership(membership_view(), max_val = bad), "max_val")
  }
  for (bad in list(0, 1.5, NA_real_, Inf)) {
    expect_error(filter_cells_by_total_membership(membership_view(), numiter = bad), "numiter")
  }
  view <- membership_view()
  view$activity[1, 1] <- Inf
  expect_error(filter_cells_by_total_membership(view), "nonfinite")
})

test_that("native fits work in either sample orientation without mutation", {
  A <- membership_view()$activity
  B <- membership_view()$effects
  for (side in c("rows", "columns")) {
    pair <- if (side == "rows") list(A, B) else list(B, A)
    priors <- list(ebnm::ebnm_point_exponential, ebnm::ebnm_normal)
    if (side == "columns") priors <- rev(priors)
    fit <- flashier::flash_factors_init(flashier::flash_init(tcrossprod(pair[[1]], pair[[2]])),
      pair, ebnm_fn = priors)
    before <- serialize(fit, NULL)
    out <- filter_cells_by_total_membership(fit, side, max_val = 1.5)
    view <- standardize_factors(fit, sample_side = side, scaling = "raw", orientation = "as_fit")
    expect_equal(tcrossprod(out$activity, out$effects),
      tcrossprod(view$activity[out$kept_samples, , drop = FALSE], view$effects))
    expect_identical(serialize(fit, NULL), before)
  }
})
