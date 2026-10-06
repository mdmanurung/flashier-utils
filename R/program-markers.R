# Independent descriptive adaptations of ImmGen-T protein interpretation.
.aligned_measurements <- function(measurements, ids) {
  measurements <- .validate_matrix(measurements,allow_missing=FALSE)
  if(!setequal(rownames(measurements),ids)) stop("Measurement sample IDs must match the declared population exactly")
  as.matrix(measurements[ids,,drop=FALSE])
}
#' Correlate activities with aligned external markers
#' @param fit Native fit or representation.
#' @param measurements Named sample-by-marker matrix.
#' @param factors Optional factor IDs.
#' @param method pearson or spearman.
#' @return Descriptive correlations, counts and constant-predictor status.
#' @author Mikhael Manurung; conceptual ImmGen-T Figure 7 workflow.
#' @export
correlate_factor_markers <- function(fit, measurements, factors = NULL, method = "pearson") {
  view <- .resolve_view(fit); A <- factor_activity(view,factors=factors)
  M <- .aligned_measurements(measurements,rownames(A)); method <- match.arg(method,c("pearson","spearman"))
  out <- list()
  for(k in seq_len(ncol(A))) for(j in seq_len(ncol(M))) {
    variable <- nrow(A)>1 && stats::sd(A[,k])>0 && stats::sd(M[,j])>0
    out[[length(out)+1L]] <- data.frame(factor=colnames(A)[k],marker=colnames(M)[j],correlation=if(variable) stats::cor(A[,k],M[,j],method=method) else NA_real_,n_samples=nrow(A),status=if(variable) "available" else "constant_predictor")
  }
  .program_result(do.call(rbind,out),view,"correlate_factor_markers",match.call(),list(method=method,population=rownames(A),inference="descriptive_no_p_values"))
}

#' Compare two existing descriptive factor contrasts
#' @param x,y Results from compare_factor_groups.
#' @param n_label Number of factors labeled by largest absolute effect.
#' @return Effect-versus-effect plot requiring the same factor coordinates and units.
#' @details Marker-positive minus marker-negative contrasts can be compared
#'   within different lineages. Each result must retain the same source basis.
#' @author Mikhael Manurung; conceptual ImmGen-T Figure 7 workflow.
#' @export
plot_factor_contrasts <- function(x, y, n_label = 15L) {
  .program_number(n_label,"n_label",TRUE)
  a <- x$factor_contrast; b <- y$factor_contrast
  if(is.null(a) || is.null(b)) stop("Supply two compare_factor_groups results")
  for(key in c("basis_id","factor_ids","feature_scale","scaling","orientation")) if(is.null(x$analysis_metadata[[key]]) || !identical(x$analysis_metadata[[key]],y$analysis_metadata[[key]])) stop("Contrast factor coordinates or units differ: ",key)
  if(!setequal(a$factor,b$factor) || anyDuplicated(a$factor) || anyDuplicated(b$factor)) stop("Contrast factor IDs differ or are duplicated")
  table <- data.frame(factor=a$factor,x=a$estimate,y=b$estimate[match(a$factor,b$factor)])
  if(any(!is.finite(table$x)) || any(!is.finite(table$y))) stop("Contrast effects must be finite")
  labels <- utils::head(order(-pmax(abs(table$x),abs(table$y)),table$factor),n_label)
  table$label <- ""; table$label[labels] <- table$factor[labels]
  p <- ggplot2::ggplot(table,ggplot2::aes(x=.data$x,y=.data$y))+ggplot2::geom_point()+ggplot2::geom_abline(slope=1,intercept=0)+.program_labels()+ggplot2::labs(x=unique(a$contrast),y=unique(b$contrast))
  attr(p,"analysis_metadata") <- utils::modifyList(x$analysis_metadata,list(engine="ggplot2",engine_version=as.character(utils::packageVersion("ggplot2"))))
  .attach_provenance(p,"plot_factor_contrasts",match.call(),list(n_label=n_label,contrasts=c(unique(a$contrast),unique(b$contrast))))
}

#' Estimate marker threshold candidates with native Gaussian mixtures
#' @param measurements Named sample-by-marker matrix.
#' @param eligible_above Only observations strictly above this value enter mixtures.
#'   The default of 0.5 suits background-subtracted counts, where values at or
#'   below it are treated as absent; it is not a general default. Set it from the
#'   assay scale (a finite value below every observation keeps all values) and
#'   review the `n_excluded` column. It is recorded in the analysis
#'   metadata.
#' @param min_observations Minimum eligible observations.
#' @param seed Scoped native fitting seed.
#' @param thresholds Optional supplied named thresholds, which remain authoritative.
#' @return Candidate and authoritative thresholds, mixture parameters, sizes,
#'   warnings, separation in pooled SD units, BIC support and review flags.
#' @details Candidate is the maximum observation classified in the lower-mean
#'   component below the lowest upper-component value. Failed, empty, poorly
#'   separated, nonmonotone or BIC-unsupported fits require manual review. Native mclustBIC is explicitly bound for namespace-only use.
#' @author Mikhael Manurung; mixture engine by Luca Scrucca, Chris Fraley,
#'   Adrian E. Raftery and collaborators. Conceptual ImmGen-T Figure 7 workflow.
#' @export
estimate_marker_thresholds <- function(measurements, eligible_above = 0.5, min_observations = 50L, seed = 42L, thresholds = NULL) {
  .need_program_engine("mclust"); M <- .validate_matrix(measurements,allow_missing=FALSE)
  .program_number(min_observations,"min_observations",TRUE,2)
  if(length(eligible_above)!=1 || !is.finite(eligible_above)) stop("eligible_above must be finite")
  if(!is.null(thresholds)) {
    .marker_thresholds(thresholds,names(thresholds))
    .select_ids(names(thresholds),colnames(M),"supplied threshold marker")
  }
  rows <- fits <- list()
  # Mclust evaluates this symbol in its caller environment when not attached.
  mclustBIC <- mclust::mclustBIC
  for(j in seq_len(ncol(M))) {
    values <- as.numeric(M[,j]); x <- values[values>eligible_above]; warnings <- character()
    candidate <- separation <- bic1 <- bic2 <- NA_real_; sizes <- integer(); parameters <- NULL
    flags <- character(); model <- NULL
    if(length(x)<min_observations) flags <- "insufficient_observations" else if(length(unique(x))<2) flags <- "constant_marker" else {
      models <- tryCatch(.with_seed(seed,withCallingHandlers({
        list(one=mclust::Mclust(x,G=1,verbose=FALSE),two=mclust::Mclust(x,G=2,verbose=FALSE))
      },warning=function(w) { warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning") })),error=function(e) { flags <<- c(flags,paste0("fit_failed: ",conditionMessage(e))); NULL })
      if(!is.null(models) && !is.null(models$two)) {
        model <- models$two; parameters <- model$parameters
        sizes <- tabulate(model$classification,nbins=2)
        means <- as.numeric(parameters$mean); lower <- which.min(means)
        if(length(means)!=2 || any(sizes==0)) flags <- c(flags,"empty_component") else {
          low <- x[model$classification==lower]; high <- x[model$classification!=lower]
          # Unequal variances let the wide component claim the upper tail; cut below the upper component.
          if(any(low>min(high))) {
            flags <- c(flags,"nonmonotone_classification")
            low <- low[low<min(high)]
          }
          candidate <- if(length(low)) max(low) else NA_real_
          variance <- rep(as.numeric(parameters$variance$sigmasq),length.out=2)
          separation <- abs(diff(means))/sqrt(sum(parameters$pro*variance))
          if(!is.finite(separation) || separation<1) flags <- c(flags,"low_separation")
        }
        bic2 <- unname(model$bic)
        if(!is.null(models$one)) bic1 <- unname(models$one$bic)
        if(!is.finite(bic1) || !is.finite(bic2) || bic2<=bic1) flags <- c(flags,"no_two_component_BIC_support")
      } else if(!length(flags)) flags <- "fit_failed"
    }
    if(length(warnings)) flags <- c(flags,"native_warning")
    marker <- colnames(M)[j]
    supplied <- if(is.null(thresholds)) NA_real_ else unname(thresholds[marker])
    rows[[j]] <- data.frame(marker=marker,candidate=candidate,threshold=if(is.finite(supplied)) supplied else candidate,source=if(is.finite(supplied)) "supplied" else "candidate",n_total=length(values),n_eligible=length(x),n_excluded=length(values)-length(x),separation=separation,bic_one=bic1,bic_two=bic2,bic_support=bic2-bic1,manual_review=length(flags)>0,status=if(length(flags)) paste(flags,collapse="; ") else "candidate_available")
    fits[[marker]] <- list(parameters=parameters,component_sizes=sizes,warnings=warnings)
  }
  out <- list(thresholds=do.call(rbind,rows),mixtures=fits)
  attr(out,"analysis_metadata") <- list(engine="mclust",eligible_above=eligible_above,min_observations=min_observations,seed=seed,excluded_definition="values <= eligible_above",supplied_thresholds=thresholds,inference="candidate_requires_review")
  .attach_provenance(out,"estimate_marker_thresholds",match.call(),list(eligible_above=eligible_above,min_observations=min_observations,seed=seed,supplied_thresholds=thresholds))
}
.marker_thresholds <- function(thresholds, markers) {
  .validate_ids(names(thresholds),"Threshold marker")
  if(!is.numeric(thresholds) || any(!is.finite(thresholds)) || !all(markers %in% names(thresholds))) stop("Missing or nonfinite marker thresholds")
  thresholds[markers]
}

#' Display candidate and supplied marker thresholds
#' @param result Result of estimate_marker_thresholds.
#' @param measurements Named sample-by-marker matrix used for display.
#' @return Faceted distributions with candidate and authoritative threshold lines.
#' @author Mikhael Manurung; conceptual ImmGen-T protein workflow; mclust authors.
#' @export
plot_marker_thresholds <- function(result, measurements) {
  M <- .validate_matrix(measurements,allow_missing=FALSE); table <- result$thresholds
  if(!setequal(table$marker,colnames(M))) stop("Marker IDs differ from threshold result")
  long <- data.frame(marker=rep(colnames(M),each=nrow(M)),value=as.vector(as.matrix(M)))
  lines <- rbind(data.frame(marker=table$marker,value=table$candidate,type="candidate"),data.frame(marker=table$marker,value=table$threshold,type="authoritative"))
  lines <- lines[is.finite(lines$value),,drop=FALSE]
  p <- ggplot2::ggplot(long,ggplot2::aes(x=.data$value))+ggplot2::geom_histogram(bins=40)+ggplot2::facet_wrap(~marker,scales="free")+ggplot2::geom_vline(data=lines,ggplot2::aes(xintercept=.data$value,colour=.data$type,linetype=.data$type))
  attr(p,"analysis_metadata") <- utils::modifyList(attr(result,"analysis_metadata"),list(engine="ggplot2",engine_version=as.character(utils::packageVersion("ggplot2"))))
  .attach_provenance(p,"plot_marker_thresholds",match.call(),list(markers=colnames(M)))
}

#' Compare explicit marker gates with an equal-size activity population
#' @param fit Native fit or representation.
#' @param measurements Named sample-by-marker matrix.
#' @param factor One factor ID.
#' @param positive,negative Marker IDs requiring values above or at/below thresholds.
#' @param thresholds Explicit named authoritative marker thresholds.
#' @param direction `"high"` (default) compares the gate with the same number of
#'   highest-activity samples; `"low"` uses the lowest-activity samples, for
#'   factors expected to mark marker-negative populations.
#' @return Overlap counts/scores, mean activities and selected IDs. Empty gates
#'   have unavailable scores. Activity ties are broken by sample ID.
#' @author Mikhael Manurung; conceptual ImmGen-T Figure 7 gate workflow.
#' @export
factor_gate_alignment <- function(fit, measurements, factor, positive = character(), negative = character(), thresholds, direction = c("high","low")) {
  direction <- match.arg(direction)
  view <- .resolve_view(fit); A <- factor_activity(view,factors=factor)
  if(ncol(A)!=1) stop("Choose exactly one factor")
  M <- .aligned_measurements(measurements,rownames(A)); markers <- c(positive,negative)
  if(!length(markers)) stop("Supply at least one positive or negative marker")
  .validate_ids(markers,"Gate marker"); .select_ids(markers,colnames(M),"marker")
  threshold <- .marker_thresholds(thresholds,markers)
  keep <- rep(TRUE,nrow(A))
  for(marker in positive) keep <- keep & M[,marker]>threshold[marker]
  for(marker in negative) keep <- keep & M[,marker]<=threshold[marker]
  gated <- rownames(A)[keep]; n <- length(gated)
  selected <- utils::head(rownames(A)[order(if(direction=="high") -A[,1] else A[,1],rownames(A))],n)
  overlap <- length(intersect(gated,selected))
  summary <- data.frame(factor=factor,direction=direction,n_population=nrow(A),n_gate=n,n_activity=length(selected),overlap=overlap,overlap_fraction=if(n) overlap/n else NA_real_,jaccard=if(n) overlap/(2*n-overlap) else NA_real_,mean_gate=if(n) mean(A[gated,1]) else NA_real_,mean_activity=if(n) mean(A[selected,1]) else NA_real_,status=if(n) "available" else "unavailable_empty_gate")
  .program_result(list(summary=summary,gate_ids=gated,activity_ids=selected,population_ids=rownames(A)),view,"factor_gate_alignment",match.call(),list(positive=positive,negative=negative,thresholds=threshold,direction=direction,operators=c(positive=">",negative="<="),tie_rule=if(direction=="high") "descending_activity_then_sample_ID" else "ascending_activity_then_sample_ID",inference="descriptive_no_p_values"))
}

#' Plot paired marker-gate and activity-ranked populations
#' @param result Result of factor_gate_alignment.
#' @param embedding Named sample-by-two-coordinate matrix.
#' @param display point or density (binned selected-cell counts).
#' @param bins Bins per axis for density.
#' @return Paired facets with common geometry and the entire declared population.
#' @author Mikhael Manurung; conceptual ImmGen-T Figure 7 workflow.
#' @export
plot_factor_gates <- function(result, embedding, display = "point", bins = 50L) {
  E <- .aligned_measurements(embedding,result$population_ids)
  if(ncol(E)!=2) stop("Embedding needs exactly two coordinates")
  display <- match.arg(display,c("point","density")); .program_number(bins,"bins",TRUE,1)
  table <- rbind(data.frame(sample=rownames(E),x=E[,1],y=E[,2],panel="marker gate",selected=rownames(E) %in% result$gate_ids),data.frame(sample=rownames(E),x=E[,1],y=E[,2],panel=if(identical(attr(result,"analysis_metadata")$direction,"low")) "lowest activity" else "highest activity",selected=rownames(E) %in% result$activity_ids))
  p <- ggplot2::ggplot(table,ggplot2::aes(x=.data$x,y=.data$y))+ggplot2::geom_point(colour="grey85",size=0.5)+ggplot2::facet_wrap(~panel)+ggplot2::coord_equal()
  p <- if(display=="point") p+ggplot2::geom_point(data=table[table$selected,,drop=FALSE],colour="#b2182b",size=0.7) else p+ggplot2::geom_bin_2d(data=table[table$selected,,drop=FALSE],bins=bins)
  attr(p,"analysis_metadata") <- utils::modifyList(attr(result,"analysis_metadata"),list(engine="ggplot2",engine_version=as.character(utils::packageVersion("ggplot2"))))
  .attach_provenance(p,"plot_factor_gates",match.call(),list(display=display,bins=bins,common_geometry=TRUE))
}
