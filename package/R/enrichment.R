# Method-specific public enrichment engines; no generic statistical claims.
.enrichment_input <- function(fit, pathways, factors, representation, signal, feature_map) {
  view <- .resolve_view(fit, representation)
  signal <- match.arg(signal, c("effects", "distinctiveness"))
  B <- if (signal == "effects") view$effects else factor_distinctiveness(view, format = "matrix")
  B <- B[, .select_ids(factors, view$manifest$factor_ids, "factor"), drop = FALSE]
  network <- NULL
  if (is.data.frame(pathways)) {
    if (!all(c("source", "target", "mor") %in% names(pathways))) stop("Network requires source, target and mor columns")
    network <- as.data.frame(pathways[, c("source", "target", "mor")])
    .validate_ids(unique(network$source), "Network source")
    .validate_ids(unique(network$target), "Network target")
    if (!is.numeric(network$mor) || any(!is.finite(network$mor))) stop("Network mor weights must be finite numeric values")
    if (anyDuplicated(network[, c("source", "target")])) stop("Network source-target edges must be unique")
    pathways <- stats::setNames(lapply(unique(network$source), function(source) network$target[network$source == source]), unique(network$source))
  }
  if (!is.list(pathways) || is.null(names(pathways))) stop("pathways must be a named list of feature IDs or a decoupleR network")
  .validate_ids(names(pathways), "Pathway")
  for (set in pathways) if (!is.character(set) || anyNA(set) || any(!nzchar(set))) stop("Pathway members must be nonempty feature IDs")
  pathways <- lapply(pathways, unique)
  excluded <- character()
  if (!is.null(feature_map)) {
    if (!is.character(feature_map)) stop("feature_map must be a named one-to-one character mapping")
    .validate_ids(names(feature_map), "Source feature map")
    .validate_ids(unname(feature_map), "Mapped feature")
    selected <- intersect(rownames(B), names(feature_map))
    excluded <- setdiff(rownames(B), selected)
    B <- B[selected, , drop = FALSE]
    rownames(B) <- unname(feature_map[selected])
  }
  sizes <- data.frame(pathway = names(pathways), input_size = lengths(pathways),
                       overlap_size = vapply(pathways, function(set) length(intersect(set, rownames(B))), integer(1)))
  list(view = view, B = B, pathways = pathways, network = network, sizes = sizes, mapping_excluded = excluded, signal = signal)
}

#' Enrich or score feature programs with a chosen native method
#'
#' fgsea produces ranked enrichment statistics and native p-values. decoupleR
#' ULM/MLM yields network activity scores. These are different estimands. SuSiE
#' is unsupported until its missing backend and audit gate are resolved.
#' @param fit Native fit or matched representation.
#' @param pathways Named feature-set list with organism/annotation attributes,
#'   or, for decoupleR, a data frame with source, target and finite signed mor
#'   weights. Source-target edges must be unique; target IDs must match features
#'   after feature_map. Extra network columns are ignored.
#' @param engine Required fgsea or decoupleR (susie is explicitly unsupported).
#' @param factors Optional factor IDs.
#' @param representation Matched view or settings.
#' @param signal effects or distinctiveness.
#' @param feature_map Optional explicitly one-to-one named feature mapping.
#' @param control Native settings; seed is scoped, decoupleR method is ulm or mlm.
#' @return Native method table plus pathway universe, native result and provenance.
#' @export
#' @details Origin: WRAPPER / call_only.
#'   Credits perform_gsea; later fgsea/decoupleR adapters; use provenance("enrich_factors")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#'   Network scores annotate estimated feature weights, not sample activities.
#'   Their sign depends on factor orientation and native p-values do not account
#'   for uncertainty in the EBMF fit. Network targets absent from the feature
#'   universe are excluded; the full mapped feature universe is retained in mat.
#'   By default center=FALSE: centering across factors changes the annotation
#'   signal and makes scores depend on the selected factors.
#' @examplesIf requireNamespace("fgsea",quietly=TRUE) && requireNamespace("BiocParallel",quietly=TRUE)
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' pathways <- list(panelA=paste0("g",1:6),panelB=paste0("g",7:12))
#' enrich_factors(view,pathways,"fgsea",factors="F1",
#'   control=list(seed=4,minSize=3L,nPermSimple=100L,eps=1e-3))$table
enrich_factors <- function(fit, pathways, engine, factors = NULL, representation = NULL,
                            signal = "effects", feature_map = NULL, control = list()) {
  if (missing(engine)) stop("Choose enrichment engine explicitly")
  engine <- match.arg(engine, c("fgsea", "decoupleR", "susie"))
  if (engine == "susie") stop("SuSiE mode is unsupported: install and audit susieR/singlecelljamboreeR before enabling this capability")
  if (is.data.frame(pathways) && engine != "decoupleR") stop("Weighted networks require engine='decoupleR'")
  if (!is.list(control) || (length(control) && (is.null(names(control)) || anyDuplicated(names(control))))) stop("control must be a uniquely named list")
  data <- .enrichment_input(fit, pathways, factors, representation, signal, feature_map)
  eligible <- data$sizes$overlap_size > 0L
  sets <- data$pathways[eligible]
  seed <- control$seed; control$seed <- NULL
  metadata <- c(data$view$manifest, list(engine = engine, signal = data$signal, universe = rownames(data$B),
    mapping_excluded = data$mapping_excluded, duplicate_policy = "unique_pathway_members_one_to_one_feature_map",
    pathway_sizes = data$sizes, organism = attr(pathways, "organism"), annotation_version = attr(pathways, "annotation_version"),
    control = control, seed = seed, ranking_ties = "native_fgsea_warning_retained"))
  if (!is.null(data$network)) {
    metadata$network_hash <- digest::digest(data$network, algo = "sha256")
    metadata$network_targets_excluded <- setdiff(data$network$target, rownames(data$B))
    metadata$duplicate_policy <- "unique_source_target_edges_one_to_one_feature_map"
  }
  if (!length(sets) || !nrow(data$B) || !ncol(data$B)) {
    table <- if (engine == "fgsea") data.frame(pathway = character(), factor = character(), NES = numeric(), pval = numeric(), padj = numeric(), size = integer()) else data.frame(source = character(), condition = character(), statistic = character(), score = numeric(), p_value = numeric())
    return(.attach_provenance(list(table = table, native = list(), analysis_metadata = c(metadata, list(status = "no_overlapping_pathways"))), "enrich_factors", match.call(), .resolved_parameters("enrich_factors", environment())))
  }
  if (engine == "fgsea") {
    .check_engine("fgsea", "fgseaMultilevel")
    .check_engine("BiocParallel", "SerialParam")
    if (is.null(control$BPPARAM) && is.null(control$nproc)) control$BPPARAM <- BiocParallel::SerialParam()
    if (any(!names(control) %in% setdiff(names(formals(fgsea::fgseaMultilevel)), c("pathways", "stats")))) stop("Unsupported fgseaMultilevel control")
    native <- .with_seed(seed, lapply(seq_len(ncol(data$B)), function(k) {
      stats <- stats::setNames(data$B[, k], rownames(data$B))
      do.call(fgsea::fgseaMultilevel, c(list(pathways = sets, stats = stats), if (!"nproc" %in% names(control)) list(nproc = 1L) else list(), control))
    }))
    table <- do.call(rbind, lapply(seq_along(native), function(k) {
      out <- as.data.frame(native[[k]])
      out$factor <- rep(colnames(data$B)[k], nrow(out))
      out
    }))
    metadata$method <- "fgseaMultilevel"
  } else {
    method <- if (is.null(control$method)) "ulm" else match.arg(control$method, c("ulm", "mlm"))
    control$method <- NULL
    symbol <- paste0("run_", method)
    .check_engine("decoupleR", symbol)
    fn <- getExportedValue("decoupleR", symbol)
    if (any(!names(control) %in% setdiff(names(formals(fn)), c("mat", "network", ".source", ".target", ".mor", ".likelihood")))) stop("Unsupported decoupleR control")
    network <- if (is.null(data$network)) do.call(rbind, lapply(names(sets), function(name) data.frame(source = name, target = intersect(sets[[name]], rownames(data$B)), mor = 1))) else
      data$network[data$network$target %in% rownames(data$B), , drop = FALSE]
    native <- do.call(fn, c(list(mat = data$B, network = network), control))
    table <- as.data.frame(native)
    if (!is.null(data$network) || !is.null(factors)) {
      table <- table[order(match(table$condition, colnames(data$B))), , drop = FALSE]
      rownames(table) <- NULL
    }
    metadata$method <- symbol
    metadata$interpretation <- "network_activity_score"
  }
  table$basis_id <- rep(data$view$manifest$basis_id, nrow(table))
  metadata$engine_version <- as.character(utils::packageVersion(engine))
  attr(table, "analysis_metadata") <- metadata
  .attach_provenance(list(table = table, native = native, analysis_metadata = metadata), "enrich_factors", match.call(), .resolved_parameters("enrich_factors", environment()))
}
