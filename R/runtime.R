# Shared runtime records stay outside the deterministic coordinate basis hash.
.resolved_parameters <- function(function_name, frame) {
  arguments <- names(formals(get(function_name, envir = asNamespace("flashier.utils"))))
  # ponytail: large scientific inputs are identified by the basis/configuration;
  # do not duplicate complete data, native fits or annotations into every result.
  arguments <- setdiff(arguments, c("...", "X", "fit", "object", "metadata", "results", "result", "reference", "target", "reference_fit", "representation", "pathways", "factor_effects", "reconstructed", "observed"))
  out <- lapply(arguments, function(name) {
    tryCatch({
      value <- get(name, envir = frame, inherits = FALSE)
      if (is.function(value)) list(callback_hash = digest::digest(value, algo = "sha256")) else value
    }, error = function(e) "not_supplied")
  })
  stats::setNames(out, arguments)
}

.attach_provenance <- function(result, function_name, call, parameters) {
  metadata <- if (is.list(result) && !is.null(result$analysis_metadata)) result$analysis_metadata else attr(result,"analysis_metadata")
  if (is.null(metadata) && inherits(result,"factor_representation")) metadata <- result$manifest
  if (is.null(metadata) && inherits(result,"flash")) metadata <- attr(result,"flashbridge_metadata")
  if (is.null(metadata)) metadata <- list()
  if (anyDuplicated(names(metadata))) metadata <- metadata[!duplicated(names(metadata),fromLast=TRUE)]
  if (is.null(metadata$engine)) {
    native_functions <- c("fit_ebmf","fit_nonnegative_ebmf","standardize_factors","factor_uncertainty","factor_intervals","factor_pve","fix_factors","refit_factors","remove_factors","reorder_factors")
    metadata$engine <- if (function_name %in% native_functions) "flashier" else if (grepl("^plot_",function_name)) "ggplot2" else "flashbridge_glue"
    metadata$engine_version <- if(metadata$engine %in% c("flashier","ggplot2")) as.character(utils::packageVersion(metadata$engine)) else as.character(utils::packageVersion("flashier.utils"))
  }
  if(is.null(metadata$engine_version)) {
    engine_package <- switch(metadata$engine,corshrink="CorShrink",pearson="stats",spearman="stats",metadata$engine)
    metadata$engine_version <- if(nzchar(system.file(package=engine_package))) as.character(utils::packageVersion(engine_package)) else "unknown"
  }
  runtime <- list(schema_version="1.0.0",package_version=as.character(utils::packageVersion("flashier.utils")),
    "function"=function_name,engine=if(is.null(metadata$engine))"flashbridge_glue" else metadata$engine,
    engine_version=if(is.null(metadata$engine_version))"not_applicable" else metadata$engine_version,
    call=call,parameters=parameters,date=format(Sys.time(),tz="UTC",usetz=TRUE),rng_kind=RNGkind(),
    basis_id=if(is.null(metadata$basis_id))"not_available" else metadata$basis_id,
    sample_side=if(is.null(metadata$sample_side))"not_available" else metadata$sample_side,
    factor_order=if(is.null(metadata$factor_ids))"not_available" else metadata$factor_ids,
    feature_scale=if(is.null(metadata$feature_scale))"unknown" else metadata$feature_scale,
    independent_unit=if(!is.null(metadata$unit_col))metadata$unit_col else if(is.null(metadata$independent_unit))"not_declared" else metadata$independent_unit)
  metadata <- utils::modifyList(metadata,runtime,keep.null=TRUE)
  if (is.list(result) && !is.null(result$analysis_metadata)) result$analysis_metadata <- metadata
  else attr(result,"analysis_metadata") <- metadata
  result
}
