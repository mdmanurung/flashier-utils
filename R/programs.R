# Independently written program interpretation. Conceptual credit: ImmGen-T
# workflows and Peter Carbonetto's pancreas annotation tutorial; no source copied.
.display_scale <- function(X, scale) {
  scale <- match.arg(scale, c("none", "max_abs"))
  if (scale == "none" || !ncol(X)) return(X)
  maximum <- apply(abs(X), 2, max)
  sweep(X, 2, ifelse(maximum == 0, 1, maximum), "/")
}
.program_result <- function(out, view, name, call, settings = list()) {
  attr(out, "analysis_metadata") <- utils::modifyList(view$manifest, settings)
  .attach_provenance(out, name, call, settings)
}
.program_number <- function(x, name, integer = FALSE, minimum = 0) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < minimum || (integer && x != as.integer(x))) stop(name, " must be a finite ", if(integer) "integer " else "number ", ">= ", minimum)
}
.factor_annotation <- function(annotation, factors) {
  .validate_ids(names(annotation), "Annotation factor")
  if (!all(factors %in% names(annotation)) || anyNA(annotation)) stop("Missing factor annotations")
  unname(annotation[factors])
}
.program_labels <- function() {
  if(requireNamespace("ggrepel",quietly=TRUE)) ggrepel::geom_text_repel(ggplot2::aes(label=.data$label),seed=42,max.overlaps=Inf)
  else ggplot2::geom_text(ggplot2::aes(label=.data$label))
}
.need_program_engine <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Install optional package '", package, "' for this workflow")
}

#' Build a signed bipartite program network
#' @param fit Native fit or representation.
#' @param factors Selected factor IDs.
#' @param n Maximum features per factor and requested sign.
#' @param direction positive, negative or both.
#' @param display_normalization max_abs or none.
#' @param threshold Minimum absolute displayed weight (inclusive).
#' @return Node and edge tables, with source basis and selection metadata.
#' @author Mikhael Manurung. Conceptual workflow: ImmGen-T GP analysis
#'   (Ziang Zhang / AgueroZZ repository; per-script authors unresolved).
#' @export
factor_network <- function(fit, factors = NULL, n = 5L, direction = "positive", display_normalization = "max_abs", threshold = 0.1) {
  view <- .resolve_view(fit)
  B <- factor_effects(view, factors = factors)
  W <- .display_scale(B, display_normalization)
  direction <- match.arg(direction, c("positive", "negative", "both"))
  .program_number(n, "n", TRUE); .program_number(threshold, "threshold")
  edges <- data.frame(from=character(), to=character(), factor=character(), feature=character(), estimate=numeric(), weight=numeric(), sign=character(), rank=integer(), basis_id=character())
  for (k in seq_len(ncol(B))) for (sign in if(direction == "both") c("positive", "negative") else direction) {
    eligible <- which(if(sign == "positive") W[, k] > 0 & W[, k] >= threshold else W[, k] < 0 & -W[, k] >= threshold)
    selected <- utils::head(eligible[order(-abs(W[eligible, k]), rownames(B)[eligible])], n)
    if (length(selected)) edges <- rbind(edges, data.frame(from=paste0("factor:",colnames(B)[k]), to=paste0("feature:",rownames(B)[selected]), factor=colnames(B)[k], feature=rownames(B)[selected], estimate=B[selected,k], weight=W[selected,k], sign=sign, rank=seq_along(selected), basis_id=view$manifest$basis_id))
  }
  genes <- unique(edges$feature)
  nodes <- data.frame(id=c(if(ncol(B)) paste0("factor:",colnames(B)) else character(),if(length(genes)) paste0("feature:",genes) else character()), label=c(colnames(B),genes), type=c(rep("factor",ncol(B)),rep("feature",length(genes))), basis_id=rep(view$manifest$basis_id,ncol(B)+length(genes)))
  .program_result(list(nodes=nodes,edges=edges),view,"factor_network",match.call(),list(n=n,direction=direction,display_normalization=display_normalization,threshold=threshold,selection_factors=colnames(B)))
}

#' Plot a program network with a scoped native layout
#' @param network Result of factor_network.
#' @param layout fr or stress (ggraph native stress implementation).
#' @param seed Scoped layout and label seed.
#' @param highlight Factor IDs to emphasize.
#' @param groups Named factor group vector.
#' @param gene_labels all, highlighted, none, or explicit gene IDs.
#' @return Plot with deterministic layout coordinates and provenance.
#' @author Mikhael Manurung; conceptual ImmGen-T workflow (Ziang Zhang / AgueroZZ repository).
#'   Layouts by igraph and ggraph authors; see credits.
#' @export
plot_factor_network <- function(network, layout = "fr", seed = 42L, highlight = NULL, groups = NULL, gene_labels = "highlighted") {
  .need_program_engine("igraph"); .need_program_engine("ggraph")
  layout <- match.arg(layout,c("fr","stress"))
  nodes <- network$nodes; edges <- network$edges
  if(!nrow(nodes)) stop("No network nodes to plot")
  factors <- nodes$label[nodes$type == "factor"]
  .select_ids(highlight, factors,"highlight factor")
  nodes$group <- nodes$type
  if (!is.null(groups)) nodes$group[nodes$type=="factor"] <- .factor_annotation(groups,factors)
  nodes$highlight <- nodes$type == "factor" & nodes$label %in% highlight
  labels <- if (identical(gene_labels,"all")) nodes$label[nodes$type=="feature"] else if (identical(gene_labels,"highlighted")) edges$feature[edges$factor %in% highlight] else if (identical(gene_labels,"none")) character() else gene_labels
  .select_ids(unique(labels),nodes$label[nodes$type=="feature"],"gene label")
  nodes$display_label <- ifelse(nodes$type=="factor" | nodes$label %in% labels,nodes$label,"")
  graph <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
  p <- .with_seed(seed, {
    coordinates <- if(layout=="fr") igraph::layout_with_fr(graph, weights=NA) else as.matrix(ggraph::create_layout(graph,layout="stress")[,c("x","y")])
    ggraph::ggraph(graph, layout="manual", x=coordinates[,1], y=coordinates[,2]) +
      (if(nrow(edges)) list(ggraph::geom_edge_link(ggplot2::aes(edge_colour=.data$sign, edge_linetype=.data$sign), alpha=0.5),
      ggraph::scale_edge_colour_manual(values=c(positive="#b2182b",negative="#2166ac")),
      ggraph::scale_edge_linetype_manual(values=c(positive="solid",negative="dashed"))) else NULL) +
      ggraph::geom_node_point(ggplot2::aes(fill=.data$group,size=.data$highlight),shape=21) +
      ggplot2::scale_size_manual(values=c(`FALSE`=2,`TRUE`=4)) +
      ggraph::geom_node_text(ggplot2::aes(label=.data$display_label)) + ggplot2::theme_void() + ggplot2::theme(plot.background=ggplot2::element_rect(fill="white",colour=NA)) + ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=0.1)) + ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=0.1))
  })
  attr(p,"analysis_metadata") <- utils::modifyList(attr(network,"analysis_metadata"),list(engine="ggraph",engine_version=as.character(utils::packageVersion("ggraph")),layout_engine=if(layout=="fr") "igraph" else "ggraph",layout_engine_version=as.character(utils::packageVersion(if(layout=="fr") "igraph" else "ggraph"))))
  .attach_provenance(p,"plot_factor_network",match.call(),list(layout=layout,seed=seed,highlight=highlight,groups=groups,gene_labels=gene_labels))
}

#' Compare feature effects with untruncated competitor differences
#' @inheritParams factor_features
#' @param compare_factors Competition pool; NULL uses all factors.
#' @return Scatter plot; data include original signed effects, raw differences
#'   and labels from the union of top effect and top difference candidates.
#' @details Unlike least-extreme margins, differences can oppose the target
#'   effect's sign. Negative direction reverses both plotting coordinates.
#' @author Mikhael Manurung. Conceptual credit: Peter Carbonetto's pancreas tutorial.
#' @export
plot_feature_rankings <- function(fit, factor, n = 15L, direction = "positive", representation = NULL, compare_factors = NULL) {
  view <- .resolve_view(fit,representation)
  .select_ids(factor,colnames(view$effects),"factor")
  if(length(factor)!=1) stop("Select exactly one factor")
  direction <- match.arg(direction,c("positive","negative")); .program_number(n,"n",TRUE)
  B <- view$effects; pool <- colnames(B)[.select_ids(compare_factors,colnames(B),"comparison factor")]
  competitors <- setdiff(pool,factor)
  competitor <- if(length(competitors)) apply(B[,competitors,drop=FALSE],1,if(direction=="positive") max else min) else rep(0,nrow(B))
  sign <- if(direction=="positive") 1 else -1
  table <- data.frame(feature=rownames(B),estimate=B[,factor],difference=B[,factor]-competitor)
  table$effect_score <- sign*table$estimate; table$difference_score <- sign*table$difference
  top <- function(x) utils::head(table$feature[order(-x,table$feature)],n)
  labels <- union(top(table$effect_score),top(table$difference_score))
  table$label <- ifelse(table$feature %in% labels,table$feature,"")
  p <- ggplot2::ggplot(table,ggplot2::aes(x=.data$effect_score,y=.data$difference_score)) + ggplot2::geom_point() + .program_labels() + ggplot2::labs(x=paste(direction,"feature effect"),y=paste(direction,"raw competitor difference"))
  .program_result(p,view,"plot_feature_rankings",match.call(),list(competition_factors=pool,direction=direction,n=n,difference_definition="untruncated_directional_extreme"))
}
