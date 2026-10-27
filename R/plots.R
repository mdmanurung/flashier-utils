# Independently written signed displays; tile/structure workflow credits in NOTICE.
#' Plot selected signed feature effects
#'
#' Heatmap or bar chart of the top features for one or more factors, selected
#' as in [factor_features()]. Independent signed-safe adaptation of the ZemmourLib tile workflow. The
#' displayed values use the declared matched representation, never positive-only.
#' @inheritParams factor_features
#' @param factors Character factor IDs to display, e.g. `c("F1", "F2")`.
#' @param style `"heatmap"` (default) or `"bar"`.
#' @param compare_factors Optional character factor IDs that compete for
#'   `"distinctive"`; `NULL` uses all other factors.
#' @param features Optional character feature IDs to show instead of a ranked
#'   selection.
#' @param direction `"absolute"` (default), `"positive"` or `"negative"`.
#' @param display_normalization `"none"` (default, effect units) or `"max_abs"`
#'   (scale each factor to max 1 for display only).
#' @param transpose Put features on the x axis (default `FALSE`).
#' @param annotation Optional named numeric vector, one value per factor, shown
#'   as a strip; e.g. correlations.
#' @param ... Optional `palette` (negative, zero, positive colours) and `limits`.
#' @return A ggplot object with display and basis metadata attributes.
#' @seealso [factor_features()], [plot_factor_activity()]
#' @importFrom ggplot2 .data
#' @export
#' @author David Zemmour.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: ADAPTED / none.
#'   Credits MyGeneTilePlot workflow; native flashier optional; use provenance("plot_factor_features")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examples
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' plot_factor_features(view, factors = "F1",n=3)
#' @param select Per-factor `largest`/`distinctive`/`both` (see
#'   [factor_features()]) or `global` for cross-factor selection by
#'   `max_k(abs(effect))` over the given factors (see
#'   [factor_features_global()]).
plot_factor_features <- function(fit, factors, representation = NULL, select = "largest", n = 20L,
    style = "heatmap", compare_factors = NULL, features = NULL, direction = "absolute",
    display_normalization = "none", transpose = FALSE, annotation = NULL, ...) {
  style <- match.arg(style, c("heatmap", "bar"))
  view <- .resolve_view(fit, representation)
  .require_biological_factors(view$manifest$factor_ids)
  .select_ids(factors, view$manifest$factor_ids, "factor")
  if (identical(select, "global")) {
    selected <- factor_features_global(view, factors, n = n, direction = direction)
  } else {
    selected <- lapply(factors, function(factor) factor_features(view, factor, select = select, n = n, direction = direction, compare_factors = compare_factors))
    selected <- do.call(rbind, selected)
  }
  if (is.null(features) && (is.null(selected) || !nrow(selected))) stop("No nonzero selected feature effects to plot")
  if (is.null(features)) features <- unique(selected$feature)
  table <- factor_effects(view, factors = factors, features = features, format = "long")
  display_normalization <- match.arg(display_normalization, c("none", "max_abs"))
  if (display_normalization == "max_abs") {
    B <- .display_scale(factor_effects(view, factors = factors), "max_abs")
    table$estimate <- as.vector(B[features, , drop = FALSE])
  }
  if (!is.null(annotation)) {
    annotation <- .factor_annotation(annotation, factors)
    table$annotation <- annotation[match(table$factor, factors)]
  }
  options <- list(...)
  if (any(!names(options) %in% c("palette", "limits"))) stop("Unknown display option")
  palette <- if (is.null(options$palette)) c("#2166ac", "white", "#b2182b") else options$palette
  if (length(palette) != 3L) stop("palette requires negative, zero and positive colors")
  p <- ggplot2::ggplot(table, ggplot2::aes(x = .data$factor, y = .data$feature, fill = .data$estimate))
  if (style == "heatmap") p <- p + ggplot2::geom_tile() else p <- ggplot2::ggplot(table, ggplot2::aes(x = .data$feature, y = .data$estimate, fill = .data$estimate)) + ggplot2::geom_col() + ggplot2::facet_wrap(~factor)
  p <- p + ggplot2::scale_fill_gradient2(low = palette[1], mid = palette[2], high = palette[3], midpoint = 0, limits = options$limits) + ggplot2::labs(fill = if(display_normalization == "none") "Feature effect" else "Effect / max absolute effect", subtitle = paste("Feature scale:", view$manifest$feature_scale))
  if (transpose) p <- p + ggplot2::coord_flip()
  if (!is.null(annotation)) p <- p + ggplot2::facet_grid(~annotation, scales = "free_x", space = "free_x")
  attr(p, "analysis_metadata") <- c(view$manifest, list(display_only = TRUE, display_normalization = display_normalization, selected_features = features, annotation = annotation, competition_factors = colnames(view$effects)[.select_ids(compare_factors, colnames(view$effects), "comparison factor")], direction = direction, transpose = transpose))
  .attach_provenance(p, "plot_factor_features", match.call(), .resolved_parameters("plot_factor_features", environment()))
}

#' Plot bounded sample activity with signed defaults
#'
#' Show factor activity across samples as a heatmap, points, distributions or a
#' structure plot, optionally grouped by a metadata column.
#' Samples before long-table expansion. Structure proportions require explicit
#' row_sum normalization and nonnegative scores; zero-sum rows are excluded and
#' recorded. Display normalization never changes inferential activities.
#' order_within_group = "mds1" reorders displayed samples within each group by
#' 1D classical MDS of the displayed row-sum-normalized factor proportions
#' (Euclidean distance); display-only, never an inferential result.
#' @param fit Native fit or a view from [standardize_factors()].
#' @param metadata Optional data frame with a `sample` column matching the
#'   activity row names; needed only with `group`.
#' @param group Optional name of a `metadata` column to group samples by.
#' @param style `"heatmap"` (default), `"point"`, `"distribution"` or
#'   `"structure"` (stacked proportions, nonnegative activities only).
#' @param factors Optional character factor IDs to display. Default all.
#' @param display_normalization `"none"` (default) or `"row_sum"` (structure
#'   style only).
#' @param max_samples Maximum samples drawn for display (default 2000).
#' @param seed Optional integer for the sample draw.
#' @param order_within_group `"input"` (default, current order) or `"mds1"`. mds1
#'   requires style = "structure" and a group column; groups keep their first
#'   displayed appearance, singleton and zero-distance groups keep input order,
#'   and remaining ties keep original displayed order. The MDS sign is
#'   canonical: the first nonzero coordinate in displayed order is positive.
#' @param group_order Optional character vector naming every displayed group
#'   once, in the order the structure plot should show them (default first
#'   appearance). Requires style = "structure" and a group column.
#' @param ... Optional `palette` and `limits`.
#' @return A ggplot object with included/excluded IDs and display metadata attributes.
#' @seealso [summarize_factor_activity()], [plot_factor_features()]
#' @export
#' @author David Zemmour.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: ADAPTED / none.
#'   Credits MyStructurePlot / ImmGenT display workflow; use provenance("plot_factor_activity")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examples
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' metadata <- data.frame(sample=rownames(X),donor=rownames(X),
#'   group=rep(c("A","B"),each=12))
#' plot_factor_activity(view, metadata = metadata,group="group",max_samples=12,seed=4)
plot_factor_activity <- function(fit, metadata = NULL, group = NULL, factors = NULL,
    style = "heatmap", display_normalization = "none", max_samples = 2000L,
    seed = NULL, order_within_group = c("input", "mds1"), group_order = NULL, ...) {
  style <- match.arg(style, c("heatmap", "point", "distribution", "structure"))
  display_normalization <- match.arg(display_normalization, c("none", "row_sum"))
  order_within_group <- match.arg(order_within_group)
  if (order_within_group == "mds1" && (style != "structure" || is.null(group))) {
    stop("mds1 ordering requires style='structure' and an explicit group column")
  }
  if (!is.null(group_order) && (style != "structure" || is.null(group))) {
    stop("group_order requires style='structure' and an explicit group column")
  }
  if (length(max_samples) != 1L || !is.finite(max_samples) || max_samples < 1 || max_samples != as.integer(max_samples)) stop("max_samples must be a positive integer")
  view <- .resolve_view(fit)
  A <- factor_activity(view, factors = factors)
  .require_biological_factors(colnames(A))
  ids <- rownames(A)
  if (is.null(metadata)) {
    if (!is.null(group)) stop("group requires metadata")
    metadata <- data.frame(sample = ids)
  }
  match <- .match_ids(ids, metadata, required_columns = c("sample", group))
  metadata <- match$metadata
  groups <- if (is.null(group)) rep("all", nrow(metadata)) else as.character(metadata[[group]])
  zero <- rep(FALSE, nrow(A))
  if (style == "structure") {
    if (display_normalization != "row_sum") stop("Structure display requires explicit display_normalization='row_sum'")
    if (any(A < 0)) stop("Structure proportions require non-negative activity")
    zero <- rowSums(A) == 0
    A <- A[!zero, , drop = FALSE]
    groups <- groups[!zero]
    if (!nrow(A)) stop("All samples have zero activity")
    A <- A / rowSums(A)
  } else if (display_normalization != "none") stop("row_sum is only supported as an explicit structure display")
  included <- seq_len(nrow(A))
  if (nrow(A) > max_samples) {
    strata <- split(included, groups)
    if (length(strata) > max_samples) stop("max_samples is smaller than number of groups")
    included <- .with_seed(seed, {
      first <- vapply(strata, function(rows) rows[sample.int(length(rows), 1L)], integer(1))
      remaining <- setdiff(seq_len(nrow(A)), first)
      sort(c(first, remaining[sample.int(length(remaining), max_samples - length(first))]))
    })
  }
  displayed <- rownames(A)[included]
  declared_order <- group_order
  group_order <- NULL
  if (style == "structure" && (order_within_group == "mds1" || !is.null(declared_order))) {
    displayed_groups <- groups[included]
    group_order <- unique(displayed_groups)
    if (!is.null(declared_order)) {
      .select_ids(declared_order, group_order, "group order")
      if (anyDuplicated(declared_order) || !setequal(declared_order, group_order)) {
        stop("group_order must name every displayed group exactly once")
      }
      group_order <- declared_order
    }
    ordered <- character()
    for (g in group_order) {
      rows <- which(displayed_groups == g)
      ids_g <- displayed[rows]
      if (order_within_group == "input" || length(rows) < 2L) {
        ordered <- c(ordered, ids_g)
        next
      }
      if (8 * length(rows)^2 > 256 * 1024^2) stop("Dense within-group distances exceed 256 MB; lower max_samples")
      distances <- as.matrix(stats::dist(A[included[rows], , drop = FALSE], method = "euclidean"))
      if (max(distances) == 0) {
        ordered <- c(ordered, ids_g)
        next
      }
      mds <- stats::cmdscale(stats::as.dist(distances), k = 1)[, 1]
      mds <- mds * .first_nonzero_sign(mds)
      ordered <- c(ordered, ids_g[order(mds, rows)])
    }
    displayed <- ordered
  }
  table <- .long_matrix(A[included, , drop = FALSE], "sample")
  table$group <- rep(groups[included], times = ncol(A))
  options <- list(...)
  if (any(!names(options) %in% c("palette", "limits"))) stop("Unknown display option")
  palette <- if (is.null(options$palette)) c("#2166ac", "white", "#b2182b") else options$palette
  if (style == "heatmap") {
    if (length(palette) != 3L) stop("Heatmap palette requires three colors")
    p <- ggplot2::ggplot(table, ggplot2::aes(x = .data$factor, y = .data$sample, fill = .data$estimate)) + ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = palette[1], mid = palette[2], high = palette[3], midpoint = 0, limits = options$limits) + ggplot2::labs(fill = "Activity")
    if (!is.null(group)) p <- p + ggplot2::facet_wrap(~group, scales = "free_y")
  } else if (style == "structure") {
    table$sample <- factor(table$sample, levels = displayed)
    p <- ggplot2::ggplot(table, ggplot2::aes(x = .data$sample, y = .data$estimate, fill = .data$factor)) + ggplot2::geom_col() + ggplot2::labs(y = "Display proportion")
    p <- p + if (!is.null(options$palette)) ggplot2::scale_fill_manual(values = options$palette) else .discrete_scale(table$factor, "fill")
    if (!is.null(options$limits)) p <- p + ggplot2::coord_cartesian(ylim = options$limits)
  } else {
    p <- ggplot2::ggplot(table, ggplot2::aes(x = .data$group, y = .data$estimate)) + ggplot2::facet_wrap(~factor) + ggplot2::labs(y = "Activity")
    p <- if (style == "point") p + ggplot2::geom_point() else p + ggplot2::geom_boxplot()
    if (!is.null(options$limits)) p <- p + ggplot2::coord_cartesian(ylim = options$limits)
  }
  p <- p + ggplot2::labs(subtitle = paste("Scaling:", view$manifest$scaling))
  ordering <- list(order_within_group = order_within_group)
  if (!is.null(group_order)) {
    ordering <- c(ordering, list(displayed_sample_order = displayed, group_order = group_order,
      group_order_source = if (is.null(declared_order)) "first_appearance" else "declared"))
    if (order_within_group == "mds1") ordering <- c(ordering, list(ordering_distance = "euclidean",
      ordering_embedding = "cmdscale_k1", ordering_sign = "first_nonzero_positive"))
  }
  attr(p, "analysis_metadata") <- c(view$manifest, list(display_only = TRUE, display_normalization = display_normalization,
    exclusions = match$exclusions, included_samples = rownames(A)[included], excluded_samples = setdiff(ids, rownames(A)[included]), zero_activity_samples = ids[zero],
    sampling = "one_per_group_then_uniform_remaining", max_samples = max_samples, seed = seed), ordering)
  .attach_provenance(p, "plot_factor_activity", match.call(), .resolved_parameters("plot_factor_activity", environment()))
}

# Sign that makes the first nonzero entry positive (1 for an all-zero vector).
.first_nonzero_sign <- function(x) {
  first <- which(x != 0)[1]
  if (is.na(first)) 1 else sign(x[first])
}
