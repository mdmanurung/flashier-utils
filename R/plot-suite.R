# Additional display functions. Conceptual workflow credit: MOFA2 downstream
# plotting (Ricard Argelaguet, Damien Arnol, Danila Bredikhin, Britta Velten).
# Independently written on ggplot2; no MOFA2 source is reused.

.display_options <- function(options, allowed = c("palette", "limits")) {
  if (length(options) && (is.null(names(options)) || any(!names(options) %in% allowed))) stop("Unknown display option")
  options
}

.label_layer <- function(mapping, data) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(mapping, data = data, seed = 42, max.overlaps = Inf, size = 3)
  } else {
    ggplot2::geom_text(mapping, data = data, size = 3, vjust = -0.5)
  }
}

#' Plot every feature effect of a factor, ranked, with the top features labelled
#'
#' Rank all features of one or more factors by signed effect and label the
#' strongest ones, so you see the whole distribution and not only a selection.
#' This is the view of a factor's loadings that [plot_factor_features()] cuts
#' down to `n` features. Display only: effects keep their declared units.
#' @inheritParams factor_features
#' @param factors Character factor IDs to display, e.g. `c("F1", "F2")`; one
#'   panel per factor.
#' @param n_label Number of features to label per factor (default 10).
#' @param direction `"absolute"` (default), `"positive"` or `"negative"`: which
#'   features are labelled.
#' @param scale `"none"` (default, effect units) or `"max_abs"` (each factor
#'   divided by its largest absolute effect, display only).
#' @return A ggplot object; `$data` holds `feature`, `factor`, `estimate`, `rank`
#'   and `label`, with display and basis metadata attributes.
#' @export
#' @author Mikhael Manurung. Conceptual workflow: MOFA2 plot_weights (Ricard Argelaguet, Damien Arnol, Danila Bredikhin, Britta Velten).
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: NEW_QOL / none.
#'   Credits MOFA2 ranked-weight display workflow; use provenance("plot_factor_weights")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @seealso [plot_factor_features()], [factor_features()]
#' @examples
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' plot_factor_weights(view, factors = c("F1", "F2"), n_label = 4)
plot_factor_weights <- function(fit, factors, representation = NULL, n_label = 10L,
                                direction = "absolute", scale = "none") {
  direction <- match.arg(direction, c("absolute", "positive", "negative"))
  scale <- match.arg(scale, c("none", "max_abs"))
  .program_number(n_label, "n_label", TRUE)
  view <- .resolve_view(fit, representation)
  B <- factor_effects(view, factors = factors)
  .require_biological_factors(colnames(B))
  B <- .display_scale(B, scale)
  table <- do.call(rbind, lapply(colnames(B), function(k) {
    b <- B[, k]
    ord <- order(-b, rownames(B))
    score <- switch(direction, absolute = abs(b), positive = b, negative = -b)
    top <- rownames(B)[order(-score, rownames(B))][seq_len(min(n_label, length(b)))]
    top <- top[score[top] > 0]
    data.frame(feature = rownames(B)[ord], factor = k, estimate = b[ord], rank = seq_along(ord),
               label = ifelse(rownames(B)[ord] %in% top, rownames(B)[ord], NA_character_))
  }))
  rownames(table) <- NULL
  p <- ggplot2::ggplot(table, ggplot2::aes(x = .data$rank, y = .data$estimate)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_point(ggplot2::aes(colour = !is.na(.data$label)), size = 1.2, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey50", `TRUE` = "#b2182b")) +
    .label_layer(ggplot2::aes(label = .data$label), table[!is.na(table$label), , drop = FALSE]) +
    ggplot2::facet_wrap(~factor, scales = "free_x") +
    ggplot2::labs(x = "Feature rank", y = if (scale == "none") "Feature effect" else "Effect / max absolute effect",
                  subtitle = paste("Feature scale:", view$manifest$feature_scale))
  attr(p, "analysis_metadata") <- c(view$manifest, list(display_only = TRUE, display_scale = scale,
    direction = direction, n_label = n_label, labelled_features = table$feature[!is.na(table$label)]))
  .attach_provenance(p, "plot_factor_weights", match.call(), .resolved_parameters("plot_factor_weights", environment()))
}

#' Plot two factor activities against each other
#'
#' Scatter plot of sample activity on one factor against another, coloured or
#' shaped by sample metadata, to see how programs co-occur across samples.
#' @inheritParams plot_factor_activity
#' @param representation For a native fit: a matched view or [standardize_factors()] settings.
#' @param factors Exactly two character factor IDs, x then y.
#' @param color,shape Optional names of `metadata` columns for point colour and
#'   point shape. Shape needs a discrete column.
#' @return A ggplot object; `$data` has one row per displayed sample with
#'   `sample`, `x`, `y` and any metadata columns used, plus metadata attributes.
#' @export
#' @author Mikhael Manurung. Conceptual workflow: MOFA2 plot_factors (Ricard Argelaguet, Damien Arnol, Danila Bredikhin, Britta Velten).
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: NEW_QOL / none.
#'   Credits MOFA2 factor-scatter display workflow; use provenance("plot_factor_scatter")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @seealso [plot_factor_activity()], [plot_factor_embedding()]
#' @examples
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' metadata <- data.frame(sample=rownames(X),group=rep(c("A","B"),each=12))
#' plot_factor_scatter(view, c("F1", "F2"), metadata = metadata, color = "group")
plot_factor_scatter <- function(fit, factors, metadata = NULL, color = NULL, shape = NULL,
                                representation = NULL, max_samples = 2000L, seed = NULL) {
  if (!is.character(factors) || length(factors) != 2L || anyDuplicated(factors)) stop("factors must be exactly two distinct factor IDs")
  .program_number(max_samples, "max_samples", TRUE, 1)
  view <- .resolve_view(fit, representation)
  A <- factor_activity(view, factors = factors)
  .require_biological_factors(colnames(A))
  ids <- rownames(A)
  if (is.null(metadata)) {
    if (!is.null(color) || !is.null(shape)) stop("color and shape require metadata")
    metadata <- data.frame(sample = ids)
  }
  match <- .match_ids(ids, metadata, required_columns = c("sample", color, shape))
  metadata <- match$metadata
  rows <- .with_seed(seed, if (nrow(A) > max_samples) sort(sample.int(nrow(A), max_samples)) else seq_len(nrow(A)))
  table <- data.frame(sample = ids[rows], x = A[rows, 1], y = A[rows, 2])
  for (column in unique(c(color, shape))) table[[column]] <- metadata[[column]][rows]
  if (!is.null(shape) && is.numeric(table[[shape]])) stop("shape requires a discrete metadata column")
  mapping <- ggplot2::aes(x = .data$x, y = .data$y)
  if (!is.null(color)) mapping$colour <- ggplot2::aes(colour = .data[[color]])$colour
  if (!is.null(shape)) mapping$shape <- ggplot2::aes(shape = .data[[shape]])$shape
  p <- ggplot2::ggplot(table, mapping) + ggplot2::geom_point(alpha = 0.8) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") + ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
    ggplot2::labs(x = paste(factors[1], "activity"), y = paste(factors[2], "activity"), colour = color, shape = shape)
  attr(p, "analysis_metadata") <- c(view$manifest, list(display_only = TRUE, factors = factors,
    included_samples = ids[rows], excluded_samples = ids[-rows], max_samples = max_samples, seed = seed,
    exclusions = match$exclusions))
  .attach_provenance(p, "plot_factor_scatter", match.call(), .resolved_parameters("plot_factor_scatter", environment()))
}

#' Plot native proportion of variance explained per factor
#'
#' Bar chart of each factor's native PVE from [factor_pve()]. Nonorthogonal
#' factors do not give an additive partition, so bars need not sum to the
#' total variance explained.
#' @inheritParams factor_pve
#' @return A ggplot object whose `$data` is the [factor_pve()] table.
#' @export
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: WRAPPER / none.
#'   Credits native pve field; use provenance("plot_factor_pve")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @seealso [factor_pve()]
#' @examples
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' plot_factor_pve(fit)
plot_factor_pve <- function(fit, factors = NULL) {
  table <- factor_pve(fit, factors = factors)
  table$factor <- factor(table$factor, levels = unique(table$factor))
  p <- ggplot2::ggplot(table, ggplot2::aes(x = .data$factor, y = .data$pve)) +
    ggplot2::geom_col(fill = "#4477AA") +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
    ggplot2::labs(x = NULL, y = "Native proportion of variance explained")
  attr(p, "analysis_metadata") <- c(attr(table, "analysis_metadata"), list(display_only = TRUE))
  .attach_provenance(p, "plot_factor_pve", match.call(), .resolved_parameters("plot_factor_pve", environment()))
}

#' Plot the observed data behind one factor
#'
#' Show the input values of a factor's strongest features across samples
#' ordered by that factor's activity, as a heatmap or as per-feature scatter
#' plots. Use it to check that a factor reflects a pattern visible in the data
#' and not only in the fitted coordinates. The data are displayed in the units
#' you supply; z-scoring is a display choice recorded in the metadata.
#' @inheritParams factor_features
#' @param X Named numeric matrix or double sparse Matrix of the observed data,
#'   normally the matrix the model was fitted on. Samples are on the side named
#'   by `sample_side`.
#' @param factor One factor ID, e.g. `"F1"`.
#' @param n Number of strongest features to show (default 20).
#' @param style `"heatmap"` (default) or `"scatter"` (activity against each
#'   feature's values, one panel per feature).
#' @param sample_side `"rows"` or `"columns"` of `X`. `NULL` (default) uses the
#'   side recorded in the view.
#' @param metadata Optional data frame with a `sample` column; needed with `group`.
#' @param group Optional `metadata` column used to facet the heatmap or colour the scatter.
#' @param direction `"absolute"` (default), `"positive"` or `"negative"`: which
#'   features count as strongest.
#' @param scale `"feature_z"` (default; centre and scale each feature for the
#'   heatmap, display only) or `"none"` (supplied units). Ignored for scatter.
#' @param max_samples Maximum samples drawn for display (default 2000).
#' @param seed Optional integer for the sample draw.
#' @return A ggplot object; `$data` has `sample`, `feature`, `value` (and
#'   `activity` for scatter), with display and basis metadata attributes.
#' @export
#' @author Mikhael Manurung. Conceptual workflow: MOFA2 plot_data_heatmap and plot_data_scatter (Ricard Argelaguet, Damien Arnol, Danila Bredikhin, Britta Velten).
#'   Full source and dependency credits: \url{https://github.com/mdmanurung/flashier-utils/blob/main/CREDITS.md}.
#' @details Origin: NEW_QOL / none.
#'   Credits MOFA2 data-level display workflow; use provenance("plot_factor_data")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @seealso [factor_features()], [plot_factor_weights()]
#' @examples
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' plot_factor_data(view, X, factor = "F1", n = 6)
#' plot_factor_data(view, X, factor = "F1", n = 4, style = "scatter")
plot_factor_data <- function(fit, X, factor, n = 20L, style = "heatmap", sample_side = NULL,
                             representation = NULL, metadata = NULL, group = NULL,
                             direction = "absolute", scale = "feature_z", max_samples = 2000L, seed = NULL) {
  style <- match.arg(style, c("heatmap", "scatter"))
  direction <- match.arg(direction, c("absolute", "positive", "negative"))
  scale <- match.arg(scale, c("feature_z", "none"))
  .program_number(n, "n", TRUE, 1); .program_number(max_samples, "max_samples", TRUE, 1)
  view <- .resolve_view(fit, representation)
  if (!(is.matrix(X) || methods::is(X, "dMatrix")) || is.null(rownames(X)) || is.null(colnames(X))) stop("X must be a named numeric matrix or double sparse Matrix")
  side <- if (is.null(sample_side)) view$manifest$sample_side else match.arg(sample_side, c("rows", "columns"))
  if (side == "columns") X <- Matrix::t(X)
  selected <- factor_features(view, factor = factor, select = "largest", direction = direction, n = n)
  features <- as.character(selected$feature)
  if (!length(features)) stop("No eligible features for this factor and direction")
  a <- factor_activity(view, factors = factor)[, 1]
  ids <- names(a)
  if (any(!ids %in% rownames(X))) stop("X lacks samples present in the fit; check sample_side and IDs")
  if (any(!features %in% colnames(X))) stop("X lacks selected features; check sample_side and IDs")
  if (is.null(metadata)) {
    if (!is.null(group)) stop("group requires metadata")
    metadata <- data.frame(sample = ids)
  }
  match <- .match_ids(ids, metadata, required_columns = c("sample", group))
  groups <- if (is.null(group)) NULL else as.character(match$metadata[[group]])
  rows <- .with_seed(seed, if (length(ids) > max_samples) sort(sample.int(length(ids), max_samples)) else seq_along(ids))
  ordered <- rows[order(-a[rows])]
  values <- as.matrix(X[ids[ordered], features, drop = FALSE])
  if (!all(is.finite(values))) stop("Selected X values must be finite")
  if (style == "heatmap" && scale == "feature_z") {
    sds <- apply(values, 2, stats::sd)
    values <- sweep(sweep(values, 2, colMeans(values), "-"), 2, ifelse(sds == 0, 1, sds), "/")
  }
  table <- .long_matrix(values, "sample", "value")
  names(table)[2] <- "feature"
  table$sample <- factor(table$sample, levels = ids[ordered])
  table$feature <- factor(table$feature, levels = rev(features))
  if (!is.null(groups)) table$group <- rep(groups[ordered], times = ncol(values))
  if (style == "heatmap") {
    p <- ggplot2::ggplot(table, ggplot2::aes(x = .data$sample, y = .data$feature, fill = .data$value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = if (scale == "feature_z") 0 else mean(values)) +
      ggplot2::labs(x = paste("Samples ordered by", factor, "activity"), y = NULL,
                    fill = if (scale == "feature_z") "Feature z-score" else "Value",
                    subtitle = paste("Feature scale:", view$manifest$feature_scale)) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    if (!is.null(groups)) p <- p + ggplot2::facet_grid(~group, scales = "free_x", space = "free_x")
  } else {
    table$activity <- rep(a[ordered], times = ncol(values))
    mapping <- ggplot2::aes(x = .data$activity, y = .data$value)
    if (!is.null(groups)) mapping$colour <- ggplot2::aes(colour = .data$group)$colour
    p <- ggplot2::ggplot(table, mapping) + ggplot2::geom_point(alpha = 0.7, size = 1.2) +
      ggplot2::facet_wrap(~feature, scales = "free_y") +
      ggplot2::labs(x = paste(factor, "activity"), y = "Observed value", colour = group,
                    subtitle = paste("Feature scale:", view$manifest$feature_scale))
  }
  attr(p, "analysis_metadata") <- c(view$manifest, list(display_only = TRUE, style = style, display_scale = if (style == "heatmap") scale else "none",
    factor = factor, selected_features = features, direction = direction, included_samples = ids[ordered],
    excluded_samples = ids[-rows], max_samples = max_samples, seed = seed, exclusions = match$exclusions))
  .attach_provenance(p, "plot_factor_data", match.call(), .resolved_parameters("plot_factor_data", environment()))
}
