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
  mds <- stats::cmdscale(d, k = 1)[, 1]
  mds <- mds * sign(mds[which(mds != 0)[1]])
  reference <- paste0("s", 1:4)[order(mds, 1:4)]
  expect_equal(meta$displayed_sample_order[1:4], reference)
  expect_equal(meta$ordering_sign, "first_nonzero_positive")
  expect_equal(meta$group_order_source, "first_appearance")
  # Canonical sign: the same samples supplied in reverse row order keep the order.
  reversed_view <- .make_representation(A[7:1, , drop = FALSE], view$effects, view$manifest)
  reversed <- plot_factor_activity(reversed_view, metadata, group = "cohort", style = "structure",
    display_normalization = "row_sum", order_within_group = "mds1")
  reversed_meta <- attr(reversed, "analysis_metadata")
  expect_equal(reversed_meta$group_order, c("c", "b", "a"))
  mds_rev <- stats::cmdscale(stats::dist(A[4:1, , drop = FALSE]), k = 1)[, 1]
  mds_rev <- mds_rev * sign(mds_rev[which(mds_rev != 0)[1]])
  expect_equal(reversed_meta$displayed_sample_order[4:7], paste0("s", 4:1)[order(mds_rev, 1:4)])
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
  # Declared group order applies with or without MDS; it must name every group once.
  declared <- plot_factor_activity(view, metadata, group = "cohort", style = "structure",
    display_normalization = "row_sum", group_order = c("c", "a", "b"))
  declared_meta <- attr(declared, "analysis_metadata")
  expect_equal(declared_meta$group_order, c("c", "a", "b"))
  expect_equal(declared_meta$group_order_source, "declared")
  expect_equal(declared_meta$displayed_sample_order, paste0("s", c(6, 7, 1:4, 5)))
  expect_null(declared_meta$ordering_embedding)
  both <- plot_factor_activity(view, metadata, group = "cohort", style = "structure",
    display_normalization = "row_sum", order_within_group = "mds1", group_order = c("c", "a", "b"))
  expect_equal(attr(both, "analysis_metadata")$displayed_sample_order, c("s6", "s7", reference, "s5"))
  expect_error(plot_factor_activity(view, metadata, group = "cohort", style = "structure",
    display_normalization = "row_sum", group_order = c("a", "b")), "every displayed group")
  expect_error(plot_factor_activity(view, metadata, group = "cohort", style = "structure",
    display_normalization = "row_sum", group_order = c("a", "a", "b", "c")), "unique")
  expect_error(plot_factor_activity(view, metadata, group = "cohort", style = "structure",
    display_normalization = "row_sum", group_order = c("a", "b", "z")), "Unknown group order")
  expect_error(plot_factor_activity(view, metadata, group = "cohort", group_order = c("a", "b", "c")), "structure")
})
test_that("activity plots record metadata rows absent from the fit", {
  metadata <- data.frame(sample=c("s3","s1","s2","s99"),cohort=c("b","a","a","a"))
  p <- plot_factor_activity(fixture_view(),metadata,group="cohort",seed=3)
  expect_equal(attr(p,"analysis_metadata")$exclusions$sample,"s99")
})

test_that("ranked weights keep every feature and label only the strongest", {
  view <- fixture_view()
  p <- plot_factor_weights(view, c("F1", "F2"), n_label = 2L)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 8L)
  expect_equal(sum(!is.na(p$data$label)), 4L)
  expect_true(all(diff(p$data$estimate[p$data$factor == "F1"]) <= 0))
  expect_true(any(p$data$estimate < 0))
  neg <- plot_factor_weights(view, "F1", n_label = 1L, direction = "negative")
  expect_equal(neg$data$label[!is.na(neg$data$label)], "g2")
  scaled <- plot_factor_weights(view, "F1", scale = "max_abs")
  expect_equal(max(abs(scaled$data$estimate)), 1)
  expect_error(plot_factor_weights(view, 1L), "character")
})

test_that("factor scatter needs two factors and joins metadata by sample", {
  view <- fixture_view()
  metadata <- data.frame(sample = c("s3", "s1", "s2"), cohort = c("b", "a", "a"), age = c(3, 1, 2))
  p <- plot_factor_scatter(view, c("F1", "F2"), metadata, color = "age", shape = "cohort")
  expect_equal(p$data$cohort, c("a", "a", "b"))
  expect_equal(p$data$x, unname(view$activity[, "F1"]))
  expect_error(plot_factor_scatter(view, "F1"), "exactly two")
  expect_error(plot_factor_scatter(view, c("F1", "F2"), color = "cohort"), "require metadata")
  expect_error(plot_factor_scatter(view, c("F1", "F2"), metadata, shape = "age"), "discrete")
  bounded <- plot_factor_scatter(view, c("F2", "F1"), max_samples = 2L, seed = 1)
  expect_equal(nrow(bounded$data), 2L)
  expect_length(attr(bounded, "analysis_metadata")$excluded_samples, 1L)
})

test_that("PVE bars use the native table", {
  fit <- fixture_fit()
  p <- plot_factor_pve(fit)
  expect_equal(as.character(p$data$factor), factor_pve(fit)$factor)
  expect_equal(p$data$pve, factor_pve(fit)$pve)
  expect_error(plot_factor_pve(fixture_view()), "PVE|pve")
})

test_that("data plots order samples by activity and respect sample_side", {
  fit <- fixture_fit()
  view <- standardize_factors(fit)
  X <- fixture_matrix()
  p <- plot_factor_data(view, X, "F1", n = 3L, scale = "none")
  expect_equal(nlevels(p$data$feature), 3L)
  expect_equal(levels(p$data$sample), names(sort(view$activity[, "F1"], decreasing = TRUE)))
  expect_equal(p$data$value[p$data$feature == levels(p$data$feature)[1]],
    unname(X[levels(p$data$sample), levels(p$data$feature)[1]]))
  z <- plot_factor_data(view, X, "F1", n = 3L)
  expect_equal(as.vector(tapply(z$data$value, z$data$feature, mean)), rep(0, 3), tolerance = 1e-12)
  flipped <- plot_factor_data(view, t(X), "F1", n = 3L, style = "scatter", sample_side = "columns")
  expect_true(all(c("activity", "value") %in% names(flipped$data)))
  expect_error(plot_factor_data(view, X[, 1:2], "F1", n = 3L), "lacks selected features")
  expect_error(plot_factor_data(view, X, "F1", group = "g"), "requires metadata")
})

test_that("embedding display accepts a plain coordinate matrix", {
  view <- fixture_view()
  emb <- matrix(c(0, 1, 2, 0, -1, 1), 3, 2, dimnames = list(c("s2", "s3", "s1"), c("d1", "d2")))
  p <- plot_factor_embedding(view, emb, c("F1", "F2"))
  expect_equal(nrow(p$data), 6L)
  expect_equal(p$data$x[p$data$sample == "s1"], c(2, 2))
  pct <- plot_factor_embedding(view, emb, "F1", display = "percentile", reference_samples = c("s1", "s2"))
  expect_equal(pct$data$estimate[pct$data$sample == "s1"], 100 * ecdf(view$activity[c("s1", "s2"), "F1"])(view$activity["s1", "F1"]))
  expect_error(plot_factor_embedding(view, emb[, 1, drop = FALSE], "F1"), "two embedding")
  expect_error(plot_factor_embedding(view, emb, "F1", activity_reduction = "x"), "Seurat")
  expect_error(plot_factor_embedding(view, "umap", "F1"), "named numeric")
})
