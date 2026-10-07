# Independently written signed displays; tile/structure workflow credits in NOTICE.
#' Plot selected signed feature effects
#'
#' Independent signed-safe adaptation of the ZemmourLib tile workflow. The
#' displayed values use the declared matched representation, never positive-only.
#' @inheritParams factor_features
#' @param factors Unique displayed factors.
#' @param style heatmap or bar.
#' @param engine ggplot2.
#' @param compare_factors Optional competition pool.
#' @param features Explicit feature IDs instead of ranked selection.
#' @param direction absolute, positive or negative selection.
#' @param display_normalization none (effect units) or max_abs (display only).
#' @param transpose Flip display axes.
#' @param annotation One named factor annotation vector, including correlations.
#' @param ... Optional palette (negative, zero, positive) and limits.
#' @return A ggplot with display and basis metadata attributes.
#' @importFrom ggplot2 .data
#' @export
#' @author David Zemmour.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
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
#' plot_factor_features(view,"F1",n=3)
plot_factor_features <- function(fit, factors, select = "largest", n = 20L,
                                 style = "heatmap", engine = "ggplot2", representation = NULL, ..., compare_factors = NULL,
                                 features = NULL, direction = "absolute", display_normalization = "none",
                                 transpose = FALSE, annotation = NULL) {
  if (engine != "ggplot2") stop("Only the ggplot2 signed display engine is supported")
  style <- match.arg(style, c("heatmap", "bar"))
  view <- .resolve_view(fit, representation)
  .select_ids(factors, view$manifest$factor_ids, "factor")
  selected <- lapply(factors, function(factor) factor_features(view, factor, select = select, n = n, direction = direction, compare_factors = compare_factors))
  selected <- do.call(rbind, selected)
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
#' Samples before long-table expansion. Structure proportions require explicit
#' row_sum normalization and nonnegative scores; zero-sum rows are excluded and
#' recorded. Display normalization never changes inferential activities.
#' @param fit Native fit or matched view.
#' @param metadata Named sample metadata; NULL is allowed without grouping.
#' @param group Optional metadata grouping column.
#' @param style heatmap, point, distribution or structure.
#' @param factors Optional unique factors.
#' @param display_normalization none, or row_sum for structure only.
#' @param engine ggplot2.
#' @param max_samples Positive bound on displayed samples.
#' @param seed Optional scoped sampling seed.
#' @param ... Optional palette and limits.
#' @return ggplot with included/excluded IDs and display metadata.
#' @export
#' @author David Zemmour.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
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
#' plot_factor_activity(view,metadata,group="group",max_samples=12,seed=4)
plot_factor_activity <- function(fit, metadata = NULL, group = NULL, style = "heatmap", factors = NULL,
                                 display_normalization = "none", engine = "ggplot2", max_samples = 2000L, seed = NULL, ...) {
  if (engine != "ggplot2") stop("Only ggplot2 is supported")
  style <- match.arg(style, c("heatmap", "point", "distribution", "structure"))
  display_normalization <- match.arg(display_normalization, c("none", "row_sum"))
  if (length(max_samples) != 1L || !is.finite(max_samples) || max_samples < 1 || max_samples != as.integer(max_samples)) stop("max_samples must be a positive integer")
  view <- .resolve_view(fit)
  A <- factor_activity(view, factors = factors)
  if (!ncol(A)) stop("No factors to display")
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
  table <- .long_matrix(A[included, , drop = FALSE], "sample")
  table$group <- rep(groups[included], times = ncol(A))
  options <- list(...)
  if (any(!names(options) %in% c("palette", "limits"))) stop("Unknown display option")
  palette <- if (is.null(options$palette)) c("#2166ac", "white", "#b2182b") else options$palette
  if (style == "heatmap") {
    if (length(palette) != 3L) stop("Heatmap palette requires three colors")
    p <- ggplot2::ggplot(table, ggplot2::aes(x = .data$factor, y = .data$sample, fill = .data$estimate)) + ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = palette[1], mid = palette[2], high = palette[3], midpoint = 0, limits = options$limits)
    if (!is.null(group)) p <- p + ggplot2::facet_wrap(~group, scales = "free_y")
  } else if (style == "structure") {
    table$sample <- factor(table$sample, levels = rownames(A)[included])
    p <- ggplot2::ggplot(table, ggplot2::aes(x = .data$sample, y = .data$estimate, fill = .data$factor)) + ggplot2::geom_col() + ggplot2::labs(y = "Display proportion")
    if (!is.null(options$palette)) p <- p + ggplot2::scale_fill_manual(values = options$palette)
    if (!is.null(options$limits)) p <- p + ggplot2::coord_cartesian(ylim = options$limits)
  } else {
    p <- ggplot2::ggplot(table, ggplot2::aes(x = .data$group, y = .data$estimate)) + ggplot2::facet_wrap(~factor)
    p <- if (style == "point") p + ggplot2::geom_point() else p + ggplot2::geom_boxplot()
    if (!is.null(options$limits)) p <- p + ggplot2::coord_cartesian(ylim = options$limits)
  }
  attr(p, "analysis_metadata") <- c(view$manifest, list(display_only = TRUE, display_normalization = display_normalization,
    included_samples = rownames(A)[included], excluded_samples = setdiff(ids, rownames(A)[included]), zero_activity_samples = ids[zero],
    sampling = "one_per_group_then_uniform_remaining", max_samples = max_samples, seed = seed))
  .attach_provenance(p, "plot_factor_activity", match.call(), .resolved_parameters("plot_factor_activity", environment()))
}
