test_that("signed displays retain estimates and bound long expansion", {
  view <- fixture_view()
  p <- plot_factor_features(view,c("F1","F2"),n=4)
  expect_s3_class(p,"ggplot")
  expect_true(any(p$data$estimate < 0))
  metadata <- data.frame(sample=c("s3","s1","s2"),cohort=c("b","a","a"))
  p <- plot_factor_activity(view,metadata,group="cohort",max_samples=2L,seed=3)
  expect_equal(length(unique(p$data$sample)),2)
  expect_equal(nrow(p$data),4)
  expect_setequal(p$data$group,c("a","b"))
  expect_error(plot_factor_activity(view,style="structure",display_normalization="row_sum"),"non-negative")
  expect_error(plot_factor_activity(view,style="structure"),"explicit")
  positive <- .make_representation(abs(view$activity),view$effects,view$manifest)
  p <- plot_factor_activity(positive,style="structure",display_normalization="row_sum")
  expect_equal(as.numeric(tapply(p$data$estimate,p$data$sample,sum)),c(1,1,1))
  expect_equal(positive$activity,abs(view$activity))
  positive$activity[1,] <- 0
  positive <- .make_representation(positive$activity,positive$effects,positive$manifest)
  p <- plot_factor_activity(positive,style="structure",display_normalization="row_sum",palette=c(F1="blue",F2="red"),limits=c(0,1))
  expect_equal(attr(p,"analysis_metadata")$zero_activity_samples,"s1")
})

## GB-005 / M8: STRUCTURE-style within-group 1D MDS ordering (display-only).
test_that("structure plots order displayed samples within groups by 1D MDS", {
  A <- matrix(c(0.9, 0.1, 0.1, 0.9, 0.5, 0.5, 0.3, 0.7, 1, 0, 0.5, 0.5, 0.5, 0.5),
    7, 2, byrow = TRUE, dimnames = list(paste0("s", 1:7), c("F1", "F2")))
  view <- .make_representation(A, matrix(rnorm(8), 4, 2, dimnames = list(paste0("g", 1:4), c("F1", "F2"))),
    list(fit_id = "ord", schema_version = "1.0.0",
    sample_side = "rows", native_indices = 1:2, scaling = "raw", orientation = "as_fit",
    constraints = c(activity = "nonnegative", effects = "signed"), feature_scale = "arcsinh_intensity"))
  metadata <- data.frame(sample = paste0("s", 1:7),
    cohort = c("a", "a", "a", "a", "b", "c", "c"))
  p <- plot_factor_activity(view, metadata, group = "cohort", style = "structure",
    display_normalization = "row_sum", order_within_group = "mds1")
  meta <- attr(p, "analysis_metadata")
  # Independent re-derivation of the documented algorithm for group a.
  d <- stats::dist(A[1:4, , drop = FALSE], method = "euclidean")
  mds <- stats::cmdscale(d, k = 1)
  reference <- paste0("s", 1:4)[order(mds[, 1], 1:4)]
  expect_equal(meta$displayed_sample_order[1:4], reference)
  # Group order is first appearance; singleton and zero-distance groups keep input order.
  expect_equal(meta$group_order, c("a", "b", "c"))
  expect_equal(meta$displayed_sample_order[5], "s5")
  expect_equal(meta$displayed_sample_order[6:7], c("s6", "s7"))
  expect_equal(meta$ordering_distance, "euclidean")
  expect_equal(meta$ordering_embedding, "cmdscale_k1")
  expect_equal(levels(p$data$sample), meta$displayed_sample_order)
  # Default keeps the current input order and records the choice.
  plain <- plot_factor_activity(view, metadata, group = "cohort", style = "structure",
    display_normalization = "row_sum")
  plain_meta <- attr(plain, "analysis_metadata")
  expect_equal(plain_meta$order_within_group, "input")
  expect_null(plain_meta$group_order)
  expect_equal(levels(plain$data$sample), paste0("s", 1:7))
  # Ordering applies to the displayed subsample only.
  subsampled <- plot_factor_activity(view, metadata, group = "cohort", style = "structure",
    display_normalization = "row_sum", order_within_group = "mds1", max_samples = 4L, seed = 4)
  sub_meta <- attr(subsampled, "analysis_metadata")
  expect_equal(length(sub_meta$displayed_sample_order), 4L)
  expect_true(all(sub_meta$displayed_sample_order %in% sub_meta$included_samples))
  # mds1 requires the structure style and an explicit group.
  expect_error(plot_factor_activity(view, metadata, group = "cohort", style = "heatmap",
    order_within_group = "mds1"), "structure")
  expect_error(plot_factor_activity(view, style = "structure", display_normalization = "row_sum",
    order_within_group = "mds1"), "group")
})
