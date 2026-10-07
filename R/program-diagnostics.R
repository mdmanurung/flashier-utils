# Descriptive program diagnostics, independently written. Workflow credit:
# ImmGen-T Figure 2; native classification metrics are delegated to pROC.
#' Describe program activity and signature sizes
#' @param fit Native fit or representation.
#' @param factors Optional factor IDs.
#' @param activity_threshold Activity threshold, strictly exceeded.
#' @param feature_threshold Absolute effect threshold, strictly exceeded.
#' @param activity_direction positive, negative or absolute.
#' @param display_normalization max_abs (per column) or none.
#' @return Per-factor and per-sample tables, thresholds and source metadata.
#' @author Mikhael Manurung; conceptual ImmGen-T workflow (Ziang Zhang / AgueroZZ repository).
#' @export
factor_diagnostics <- function(fit, factors = NULL, activity_threshold = 0.1, feature_threshold = 0.25, activity_direction = "positive", display_normalization = "max_abs") {
  view <- .resolve_view(fit)
  A <- .display_scale(factor_activity(view,factors=factors),display_normalization)
  B <- .display_scale(factor_effects(view,factors=factors),display_normalization)
  activity_direction <- match.arg(activity_direction,c("positive","negative","absolute"))
  .program_number(activity_threshold,"activity_threshold"); .program_number(feature_threshold,"feature_threshold")
  score <- switch(activity_direction,positive=A,negative=-A,absolute=abs(A))
  active <- score > activity_threshold
  per_factor <- data.frame(factor=colnames(A),active_proportion=colMeans(active),active_cells=colSums(active),active_features=colSums(abs(B)>feature_threshold),n_samples=rep(nrow(A),ncol(A)),n_features=rep(nrow(B),ncol(A)))
  per_sample <- data.frame(sample=rownames(A),programs=rowSums(active))
  .program_result(list(factors=per_factor,samples=per_sample),view,"factor_diagnostics",match.call(),list(activity_threshold=activity_threshold,feature_threshold=feature_threshold,activity_direction=activity_direction,display_normalization=display_normalization,threshold_operator=">",population=rownames(A),selection_factors=colnames(A),inference="descriptive"))
}

#' Plot the complete program diagnostic workflow
#' @param fit Native fit or representation.
#' @param metadata Optional sample metadata.
#' @param group Optional grouping column.
#' @param marker Optional named external-marker vector.
#' @param factor Factor for external-marker scatter.
#' @param max_samples Maximum scatter rendering population; statistics use all samples.
#' @param seed Scoped rendering seed.
#' @param ... Arguments passed to factor_diagnostics.
#' @return Named plots: prevalence, features, programs_per_cell, groups and marker.
#' @author Mikhael Manurung; conceptual ImmGen-T Figure 2 workflow.
#' @export
plot_factor_diagnostics <- function(fit, metadata = NULL, group = NULL, marker = NULL, factor = NULL, max_samples = 2000L, seed = 42L, ...) {
  view <- .resolve_view(fit); diagnostics <- factor_diagnostics(view,...)
  .program_number(max_samples,"max_samples",TRUE,1)
  f <- diagnostics$factors; s <- diagnostics$samples
  plots <- list(prevalence=ggplot2::ggplot(f,ggplot2::aes(x=.data$active_proportion))+ggplot2::geom_histogram(bins=30)+ggplot2::labs(x="Active sample proportion"),
    features=ggplot2::ggplot(f,ggplot2::aes(x=.data$active_features))+ggplot2::geom_histogram(bins=30)+ggplot2::labs(x="Active feature count"),
    programs_per_cell=ggplot2::ggplot(s,ggplot2::aes(x=.data$programs))+ggplot2::geom_histogram(binwidth=1)+ggplot2::labs(x="Programs per sample"))
  if(!is.null(group)) {
    md <- .match_ids(s$sample,metadata,required_columns=c("sample",group))$metadata
    s$group <- md[[group]]
    plots$groups <- ggplot2::ggplot(s,ggplot2::aes(x=.data$group,y=.data$programs))+ggplot2::geom_boxplot()
  }
  if(!is.null(marker)) {
    if(is.null(factor) || length(factor)!=1) stop("Choose one factor for the marker scatter")
    A <- factor_activity(view,factors=factor)
    .validate_ids(names(marker),"Marker sample")
    if(!setequal(names(marker),rownames(A)) || !is.numeric(marker) || any(!is.finite(marker))) stop("Marker must be finite and match all sample IDs")
    marker <- marker[rownames(A)]
    population <- data.frame(sample=rownames(A),activity=A[,1],marker=unname(marker))
    index <- .with_seed(seed,sort(sample.int(nrow(A),min(nrow(A),max_samples))))
    plots$marker <- ggplot2::ggplot(population[index,,drop=FALSE],ggplot2::aes(x=.data$activity,y=.data$marker))+ggplot2::geom_point()
    attr(plots$marker,"full_population_correlation") <- if(stats::sd(A[,1])>0 && stats::sd(marker)>0) stats::cor(A[,1],marker) else NA_real_
  }
  plots <- lapply(plots,function(p) .program_result(p,view,"plot_factor_diagnostics",match.call(),list(diagnostics=attr(diagnostics,"analysis_metadata"),max_samples=max_samples,seed=seed,rendering_only=TRUE)))
  .program_result(plots,view,"plot_factor_diagnostics",match.call(),list(diagnostics=attr(diagnostics,"analysis_metadata")))
}

#' Plot signed signature weights against mean expression
#' @param fit Native fit or representation.
#' @param factor One factor ID.
#' @param mean_expression Finite named feature means.
#' @param expression_scale Required description of expression units.
#' @param threshold Minimum absolute max-absolute-scaled weight for labels.
#' @param n_label Maximum labels by absolute weight then feature ID.
#' @return Plot with expression on the vertical axis, never significance.
#' @author Mikhael Manurung; conceptual ImmGen-T Figure 2 workflow.
#' @export
plot_factor_signature <- function(fit, factor, mean_expression, expression_scale, threshold = 0.1, n_label = 15L) {
  view <- .resolve_view(fit); B <- factor_effects(view,factors=factor)
  if(ncol(B)!=1) stop("Choose one factor")
  .program_number(threshold,"threshold"); .program_number(n_label,"n_label",TRUE)
  .validate_ids(names(mean_expression),"Expression feature")
  if(!is.numeric(mean_expression) || any(!is.finite(mean_expression)) || !all(rownames(B) %in% names(mean_expression))) stop("Supply finite expression for every feature")
  .measurement_scale(expression_scale)
  W <- .display_scale(B,"max_abs")
  table <- data.frame(feature=rownames(B),weight=W[,1],estimate=B[,1],expression=mean_expression[rownames(B)])
  eligible <- which(abs(table$weight)>=threshold & table$weight!=0)
  labels <- utils::head(eligible[order(-abs(table$weight[eligible]),table$feature[eligible])],n_label)
  table$label <- ""; table$label[labels] <- table$feature[labels]
  p <- ggplot2::ggplot(table,ggplot2::aes(x=.data$weight,y=.data$expression))+ggplot2::geom_point()+.program_labels()+ggplot2::labs(x="Signed weight / maximum absolute weight",y=paste("Mean expression:",expression_scale))
  .program_result(p,view,"plot_factor_signature",match.call(),list(expression_scale=expression_scale,threshold=threshold,n_label=n_label,display_normalization="max_abs"))
}

#' Measure descriptive one-versus-rest factor specificity
#' @param fit Native fit or representation.
#' @param metadata Metadata matching every sample.
#' @param group Grouping column.
#' @param groups Optional target labels, including absent classes.
#' @param factors Optional factor IDs.
#' @param direction positive (higher predicts target) or negative.
#' @param threshold_method closest.topleft or youden.
#' @param unit_col Optional unit IDs for counts; no inferential p-values.
#' @return AUC, finite optimal threshold, sensitivity, specificity, tie count,
#'   observation/unit counts and explicit status for every factor/group pair.
#' @author Mikhael Manurung; pROC authors Xavier Robin and collaborators.
#' @export
factor_specificity <- function(fit, metadata, group, groups = NULL, factors = NULL, direction = "positive", threshold_method = "closest.topleft", unit_col = NULL) {
  .check_engine("pROC"); view <- .resolve_view(fit)
  A <- factor_activity(view,factors=factors)
  md <- .match_ids(rownames(A),metadata,required_columns=c("sample",group,unit_col))$metadata
  direction <- match.arg(direction,c("positive","negative")); threshold_method <- match.arg(threshold_method,c("closest.topleft","youden"))
  labels <- as.character(md[[group]])
  if(is.null(groups)) groups <- unique(labels)
  .validate_ids(groups,"Group")
  out <- list()
  for(k in seq_len(ncol(A))) for(g in groups) {
    target <- labels==g; x <- A[,k]
    row <- data.frame(factor=colnames(A)[k],group=g,auc=NA_real_,threshold=NA_real_,sensitivity=NA_real_,specificity=NA_real_,threshold_ties=0L,n_samples=length(x),n_positive=sum(target),n_negative=sum(!target),n_units=if(is.null(unit_col)) NA_integer_ else length(unique(md[[unit_col]])),n_positive_units=if(is.null(unit_col)) NA_integer_ else length(unique(md[[unit_col]][target])),n_negative_units=if(is.null(unit_col)) NA_integer_ else length(unique(md[[unit_col]][!target])),status="unavailable_absent_class")
    if(any(target) && any(!target)) {
      if(length(unique(x))==1) { row$auc <- 0.5; row$status <- "constant_predictor" } else {
        roc <- pROC::roc(target,x,levels=c(FALSE,TRUE),direction=if(direction=="positive") "<" else ">",quiet=TRUE)
        best <- as.data.frame(pROC::coords(roc,"best",best.method=threshold_method,ret=c("threshold","sensitivity","specificity"),transpose=FALSE))
        row$auc <- as.numeric(pROC::auc(roc)); row$threshold_ties <- nrow(best)
        finite <- best[is.finite(best$threshold),,drop=FALSE]
        if(nrow(finite)) {
          chosen <- finite[which.min(finite$threshold),,drop=FALSE]
          row[,c("threshold","sensitivity","specificity")] <- chosen[,c("threshold","sensitivity","specificity")]
          row$status <- "available"
        } else row$status <- "unavailable_no_finite_optimum"
      }
    }
    out[[length(out)+1L]] <- row
  }
  .program_result(do.call(rbind,out),view,"factor_specificity",match.call(),list(engine="pROC",direction=direction,threshold_method=threshold_method,tie_rule="smallest_finite_threshold",unit_col=unit_col,inference="descriptive",groups=groups))
}

#' Summarize defining features and descriptive program diagnostics
#' @param fit Native fit or representation.
#' @param factors Optional factor IDs.
#' @param n Defining features per sign.
#' @param annotations Optional named user-authored factor labels.
#' @param compare_factors Optional comparison pool.
#' @param ... Threshold settings passed to factor_diagnostics.
#' @return Summary table with list columns of positive and negative features.
#' @author Mikhael Manurung; feature ranking credits Peter Carbonetto and Matthew Stephens.
#' @export
summarize_factors <- function(fit, factors = NULL, n = 5L, annotations = NULL, compare_factors = NULL, ...) {
  view <- .resolve_view(fit); out <- factor_diagnostics(view,factors=factors,...)$factors
  for(sign in c("positive","negative")) out[[paste0(sign,"_features")]] <- lapply(out$factor,function(factor) factor_features(view,factor,n=n,direction=sign,compare_factors=compare_factors))
  if(!is.null(annotations)) out$annotation <- .factor_annotation(annotations,out$factor)
  .program_result(out,view,"summarize_factors",match.call(),list(n=n,annotations=annotations,competition_factors=colnames(view$effects)[.select_ids(compare_factors,colnames(view$effects),"comparison factor")],diagnostics=attr(factor_diagnostics(view,factors=factors,...),"analysis_metadata")))
}
