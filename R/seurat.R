# Public SeuratObject/uwot calls; independent glue credits ImmGenT in NOTICE.
#' Add exactly aligned cell activities to a Seurat reduction
#' @param object Seurat object with cell IDs equal to activity sample IDs.
#' @param fit Native fit or matched representation, at cell level.
#' @param reduction Unique alphanumeric reduction name starting with a letter.
#' @param assay Assay name; NULL uses DefaultAssay.
#' @param overwrite Explicitly replace an existing reduction.
#' @return Seurat object with basis metadata in reduction Misc. Loadings are added
#'   only when the assay and effects have exactly the same feature IDs.
#' @export
#' @author David Zemmour / Zemmour Lab (ImmGenT workflows); Paul Hoffman; Rahul Satija; David Collins; Yuhan Hao; Austin Hartman; Gesmira Molla; Andrew Butler; Tim Stuart; Madeline Kowalski; Saket Choudhary; Skylar Li; Longda Jiang; Anagha Shenoy; Jeff Farrell; Shiwei Zheng; Christoph Hafemeister; Patrick Roelli (SeuratObject).
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: ADAPTED / call_only.
#'   Credits CreateDimReducObject workflow; AddLatentData conceptual context; use provenance("add_factors_to_seurat")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examplesIf requireNamespace("SeuratObject",quietly=TRUE)
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' # Adapter mechanics on simulated cells; no biological inference.
#' counts <- Matrix::Matrix(matrix(1,16,24,dimnames=list(colnames(X),rownames(X))),sparse=TRUE)
#' object <- SeuratObject::CreateSeuratObject(counts,assay="panel")
#' object <- add_factors_to_seurat(object,view)
#' SeuratObject::Embeddings(object[["ebmf"]])[1:3,,drop=FALSE]
add_factors_to_seurat <- function(object, fit, reduction = "ebmf", assay = NULL, overwrite = FALSE) {
  .check_engine("SeuratObject", c("CreateDimReducObject", "DefaultAssay", "Misc", "Reductions"))
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  .reduction_name(reduction)
  if (reduction %in% SeuratObject::Reductions(object) && !isTRUE(overwrite)) stop("Reduction already exists; request overwrite explicitly")
  view <- .resolve_view(fit)
  cells <- colnames(object)
  if (!setequal(cells, rownames(view$activity))) stop("Cell IDs must exactly match activity sample IDs; donor-to-cell expansion is unsupported")
  if (!ncol(view$activity)) stop("A reduction requires at least one factor")
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
  if (length(assay) != 1L || !assay %in% SeuratObject::Assays(object)) stop("Unknown assay")
  A <- view$activity[cells, , drop = FALSE]
  key <- paste0(reduction, "_")
  dimensions <- paste0(key, seq_len(ncol(A)))
  colnames(A) <- dimensions
  features <- rownames(object[[assay]])
  compatible <- setequal(features, rownames(view$effects))
  B <- if (compatible) view$effects[features, , drop = FALSE] else matrix(numeric(), 0L, 0L)
  if (compatible) colnames(B) <- dimensions
  metadata <- c(view$manifest, list(display_only = FALSE, cell_alignment = cells, assay = assay,
    loading_status = if (compatible) "matching_assay_feature_ids" else "omitted_incompatible_feature_ids",
    dimension_factor_map = stats::setNames(colnames(view$activity), dimensions), reduction = reduction))
  object[[reduction]] <- SeuratObject::CreateDimReducObject(embeddings = A, loadings = B,
    assay = assay, key = key, misc = list(flashier.utils = metadata))
  .attach_provenance(object, "add_factors_to_seurat", match.call(), .resolved_parameters("add_factors_to_seurat", environment()))
}

.reduction_name <- function(name) {
  if (length(name) != 1L || is.na(name) || !grepl("^[A-Za-z][A-Za-z0-9]*$", name)) stop("Reduction name must be alphanumeric and start with a letter")
}

#' Embed declared factor coordinates using native uwot
#' @param object Seurat object containing a factor reduction.
#' @param reduction Source factor reduction.
#' @param method Only umap is supported.
#' @param dims Explicit unique positive source dimension indices; NULL uses all.
#' @param seed Optional scoped RNG seed, also passed to native uwot.
#' @param overwrite Replace an existing output reduction explicitly.
#' @param ... Native uwot settings, plus reduction_name (default source name + umap).
#' @return Seurat object with exploratory embedding and resolved settings in Misc.
#' @export
#' @author David Zemmour / Zemmour Lab (ImmGenT workflows); James Melville and uwot contributors (full list in credits article).
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: WRAPPER / call_only.
#'   Credits Seurat/uwot public embedding APIs; use provenance("embed_factors")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examplesIf requireNamespace("SeuratObject",quietly=TRUE) && requireNamespace("uwot",quietly=TRUE)
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' # Adapter mechanics on simulated cells; no biological inference.
#' counts <- Matrix::Matrix(matrix(1,16,24,dimnames=list(colnames(X),rownames(X))),sparse=TRUE)
#' object <- SeuratObject::CreateSeuratObject(counts,assay="panel")
#' object <- add_factors_to_seurat(object,view)
#' object <- embed_factors(object,seed=4,n_neighbors=4,n_epochs=20,init="random")
#' SeuratObject::Embeddings(object[["ebmfumap"]])[1:3,,drop=FALSE]
embed_factors <- function(object, reduction = "ebmf", method = "umap", dims = NULL, seed = NULL, overwrite = FALSE, ...) {
  if (method != "umap") stop("Only native uwot umap is supported")
  .check_engine("SeuratObject", c("Embeddings", "Misc", "CreateDimReducObject"))
  .check_engine("uwot", "umap")
  if (!inherits(object, "Seurat") || !reduction %in% SeuratObject::Reductions(object)) stop("Missing factor reduction")
  metadata <- SeuratObject::Misc(object[[reduction]], "flashier.utils")
  if (!is.list(metadata) || isTRUE(metadata$display_only)) stop("Embedding requires original factor activity coordinates")
  A <- SeuratObject::Embeddings(object[[reduction]])
  if (is.null(dims)) dims <- seq_len(ncol(A))
  if (!is.numeric(dims) || !length(dims) || anyNA(dims) || anyDuplicated(dims) || any(dims != as.integer(dims)) || any(!dims %in% seq_len(ncol(A)))) stop("dims must identify unique source dimensions")
  settings <- list(...)
  output <- if (is.null(settings$reduction_name)) paste0(reduction, "umap") else settings$reduction_name
  settings$reduction_name <- NULL
  .reduction_name(output)
  if (output %in% SeuratObject::Reductions(object) && !isTRUE(overwrite)) stop("Embedding reduction already exists")
  allowed <- setdiff(names(formals(uwot::umap)), c("X", "seed", "ret_model", "ret_nn", "ret_extra", "y", "epoch_callback"))
  if (length(settings) && (is.null(names(settings)) || anyDuplicated(names(settings)) || any(!names(settings) %in% allowed))) stop("Unsupported native embedding setting")
  defaults <- list(metric = "euclidean", n_neighbors = min(15L, nrow(A)-1L), n_components = 2L,
                   n_threads = 1L, n_sgd_threads = 1L, verbose = FALSE)
  settings <- utils::modifyList(defaults, settings)
  if (nrow(A) < 4L) stop("UMAP requires at least four cells")
  coordinates <- .with_seed(seed, do.call(uwot::umap, c(list(X = A[, dims, drop = FALSE], seed = seed), settings)))
  rownames(coordinates) <- rownames(A)
  key <- paste0(output, "_"); colnames(coordinates) <- paste0(key, seq_len(ncol(coordinates)))
  embed_metadata <- list(display_only = TRUE, source_reduction = reduction, source_basis_id = metadata$basis_id,
    source_factor_ids = unname(metadata$dimension_factor_map[dims]), dims = dims, seed = seed,
    engine = "uwot", engine_version = as.character(utils::packageVersion("uwot")), settings = settings,
    interpretation = "exploratory_geometry_not_batch_correction_or_replication")
  object[[output]] <- SeuratObject::CreateDimReducObject(embeddings = coordinates,
    assay = SeuratObject::DefaultAssay(object[[reduction]]), key = key, misc = list(flashier.utils = embed_metadata))
  .attach_provenance(object, "embed_factors", match.call(), .resolved_parameters("embed_factors", environment()))
}

#' Display factor activity on a cell embedding
#' @param object Seurat object with original factor and embedding reductions.
#' @param reduction Embedding reduction name.
#' @param factors Original factor IDs to display.
#' @param display activity or percentile (explicitly display-only).
#' @param ... activity_reduction, reference_samples for percentile, max_samples,
#'   seed, palette (negative/zero/positive), limits. Defaults bound 2000 cells.
#' @return ggplot with sampled IDs, reference IDs and basis metadata.
#' @export
#' @author David Zemmour / Zemmour Lab (ImmGenT workflows).
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: NEW_QOL / call_only.
#'   Credits factor-space activity display workflow; use provenance("plot_factor_embedding")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examplesIf requireNamespace("SeuratObject",quietly=TRUE) && requireNamespace("uwot",quietly=TRUE)
#' set.seed(3)
#' X <- tcrossprod(matrix(rnorm(48),24,2), matrix(rnorm(32),16,2)) +
#'   matrix(rnorm(384,sd=0.2),24,16)
#' dimnames(X) <- list(paste0("s",1:24),paste0("g",1:16))
#' fit <- fit_ebmf(X,"rows",max_factors=2,seed=4,verbose=0,
#'   ebnm_fn=list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace),
#'   nullcheck=FALSE,backfit=TRUE,feature_scale="simulated_centered_intensity")
#' view <- standardize_factors(fit)
#' # Adapter mechanics on simulated cells; no biological inference.
#' counts <- Matrix::Matrix(matrix(1,16,24,dimnames=list(colnames(X),rownames(X))),sparse=TRUE)
#' object <- SeuratObject::CreateSeuratObject(counts,assay="panel")
#' object <- add_factors_to_seurat(object,view)
#' object <- embed_factors(object,seed=4,n_neighbors=4,n_epochs=20,init="random")
#' plot_factor_embedding(object,"ebmfumap","F1",display="percentile",max_samples=12,seed=4)
plot_factor_embedding <- function(object, reduction, factors, display = "activity", ...) {
  .check_engine("SeuratObject", c("Embeddings", "Misc"))
  display <- match.arg(display, c("activity", "percentile"))
  if (!inherits(object, "Seurat") || !reduction %in% SeuratObject::Reductions(object)) stop("Missing embedding reduction")
  options <- list(...)
  if (length(options) && (is.null(names(options)) || any(!names(options) %in% c("activity_reduction", "reference_samples", "max_samples", "seed", "palette", "limits")))) stop("Unknown display setting")
  geometry <- SeuratObject::Embeddings(object[[reduction]])
  if (ncol(geometry) < 2L) stop("Display requires two embedding coordinates")
  source <- options$activity_reduction
  if (is.null(source)) source <- SeuratObject::Misc(object[[reduction]], "flashier.utils")$source_reduction
  if (is.null(source) || !source %in% SeuratObject::Reductions(object)) stop("Declare activity_reduction explicitly")
  metadata <- SeuratObject::Misc(object[[source]], "flashier.utils")
  if (!is.list(metadata) || isTRUE(metadata$display_only)) stop("Source must contain original activity")
  A <- SeuratObject::Embeddings(object[[source]])
  colnames(A) <- unname(metadata$dimension_factor_map)
  A <- A[, .select_ids(factors, colnames(A), "factor"), drop = FALSE]
  if (!setequal(rownames(A), rownames(geometry))) stop("Activity and embedding cells differ")
  reference <- options$reference_samples
  if (is.null(reference)) reference <- rownames(A)
  if (display == "percentile") {
    .select_ids(reference, rownames(A), "reference sample")
    A <- vapply(seq_len(ncol(A)), function(k) 100*stats::ecdf(A[reference,k])(A[,k]), numeric(nrow(A))) |>
      matrix(nrow = nrow(A), dimnames = dimnames(A))
  }
  bound <- if (is.null(options$max_samples)) 2000L else options$max_samples
  if (length(bound) != 1L || !is.finite(bound) || bound < 1L || bound != as.integer(bound)) stop("max_samples must be a positive integer")
  rows <- .with_seed(options$seed, if (nrow(A) > bound) sort(sample.int(nrow(A),bound)) else seq_len(nrow(A)))
  ids <- rownames(A)[rows]
  table <- .long_matrix(A[rows,,drop=FALSE],"sample")
  table$x <- rep(geometry[ids,1],ncol(A)); table$y <- rep(geometry[ids,2],ncol(A))
  palette <- if (is.null(options$palette)) c("#2166ac","white","#b2182b") else options$palette
  if (length(palette) != 3L) stop("palette requires three colors")
  p <- ggplot2::ggplot(table,ggplot2::aes(.data$x,.data$y,color=.data$estimate)) + ggplot2::geom_point() + ggplot2::facet_wrap(~factor) +
    ggplot2::labs(color = if (display == "activity") "Activity" else "Display percentile")
  p <- if (display == "activity") p + ggplot2::scale_color_gradient2(low=palette[1],mid=palette[2],high=palette[3],midpoint=0,limits=options$limits) else
    p + ggplot2::scale_color_gradient(low=palette[1],high=palette[3],limits=c(0,100))
  attr(p,"analysis_metadata") <- utils::modifyList(metadata,list(display_only=TRUE,display=display,reference_sample_ids=if(display=="percentile")reference else NULL,
    included_samples=ids,excluded_samples=setdiff(rownames(A),ids),max_samples=bound,seed=options$seed))
  .attach_provenance(p, "plot_factor_embedding", match.call(), .resolved_parameters("plot_factor_embedding", environment()))
}

#' Compute explicitly display-only empirical activity percentiles
#' @param fit Native fit or matched representation.
#' @param reference_samples Named reference sample IDs; NULL uses all samples.
#' @param factors Selected original factor IDs.
#' @param na_action Only fail is supported for finite factor representations.
#' @return Matrix with display-only class and reference/basis metadata. Values are
#'   100 * ecdf(reference)(activity); ties retain native ECDF behavior.
#' @export
#' @author David Zemmour / Zemmour Lab (ImmGenT workflows).
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: ADAPTED / call_only.
#'   Credits stats::ecdf display workflow; use provenance("factor_activity_percentile")
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
#' head(factor_activity_percentile(view,reference_samples=rownames(X)[1:12],factors="F1"))
factor_activity_percentile <- function(fit, reference_samples = NULL, factors = NULL, na_action = "fail") {
  if (na_action != "fail") stop("Only finite activities with na_action='fail' are supported")
  view <- .resolve_view(fit)
  A <- factor_activity(view,factors=factors)
  if (!ncol(A)) stop("No selected factors")
  if (is.null(reference_samples)) reference_samples <- rownames(A)
  .select_ids(reference_samples,rownames(A),"reference sample")
  result <- matrix(vapply(seq_len(ncol(A)),function(k)100*stats::ecdf(A[reference_samples,k])(A[,k]),numeric(nrow(A))),
                   nrow=nrow(A),dimnames=dimnames(A))
  attr(result,"analysis_metadata") <- c(view$manifest,list(display_only=TRUE,reference_sample_ids=reference_samples,
    ties="native_ecdf_right_continuous",interpretation="display_percentiles_not_composition_or_inference"))
  class(result) <- c("factor_activity_percentile",class(result))
  .attach_provenance(result, "factor_activity_percentile", match.call(), .resolved_parameters("factor_activity_percentile", environment()))
}
