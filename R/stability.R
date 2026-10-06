# Newly written biological-unit bootstrap and atomic checkpoint orchestration.
.bootstrap_units <- function(metadata, unit_col, strata, seed, id_col = "sample") {
  if (!unit_col %in% names(metadata)) stop("Missing biological unit column")
  units <- as.character(metadata[[unit_col]])
  if (anyNA(units) || any(!nzchar(trimws(units)))) stop("Unit IDs must be nonempty")
  if (!all(strata %in% names(metadata))) stop("Missing strata columns")
  unique_units <- unique(units)
  labels <- if (!length(strata)) rep("all", nrow(metadata)) else as.character(do.call(interaction, c(metadata[strata], list(drop = TRUE, lex.order = TRUE))))
  if (anyNA(labels)) stop("Strata cannot be missing")
  unit_strata <- vapply(unique_units, function(unit) {
    values <- unique(labels[units == unit])
    if (length(values) != 1L) stop("A biological unit crosses strata; define a valid donor-level stratum")
    values
  }, character(1))
  drawn <- .with_seed(seed, unlist(lapply(split(unique_units, unit_strata), function(group) group[sample.int(length(group), length(group), replace = TRUE)]), use.names = FALSE))
  md <- list(); mapping <- list()
  for (i in seq_along(drawn)) {
    rows <- which(units == drawn[i])
    selected <- metadata[rows, , drop = FALSE]
    original <- as.character(selected[[id_col]])
    occurrence <- paste0("bootstrap_unit_", i)
    selected[[unit_col]] <- occurrence
    selected[[id_col]] <- paste0("bootstrap_sample_", i, "_", seq_along(rows))
    md[[i]] <- selected
    mapping[[i]] <- data.frame(original_sample = original, sample = selected[[id_col]], original_unit = drawn[i], unit_occurrence = occurrence, source_row = rows, stratum = labels[rows])
  }
  list(metadata = do.call(rbind, md), mapping = do.call(rbind, mapping))
}

.resample_aux <- function(aux, mapping, X, sample_side) {
  if (!is.list(aux)) stop("observation_aux must be a named list")
  sample_ids <- if (sample_side == "rows") rownames(X) else colnames(X)
  feature_ids <- if (sample_side == "rows") colnames(X) else rownames(X)
  lapply(aux, function(value) {
    axis <- NULL
    if (is.list(value) && all(c("value", "axis") %in% names(value))) { axis <- value$axis; value <- value$value }
    if (length(value) == 1L) return(value)
    if (is.matrix(value)) {
      if (!identical(dimnames(value), dimnames(X))) stop("Matrix observation auxiliary must have exactly X's named axes")
      if (sample_side == "rows") {
        result <- value[match(mapping$original_sample, sample_ids), , drop = FALSE]; rownames(result) <- mapping$sample
      } else {
        result <- value[, match(mapping$original_sample, sample_ids), drop = FALSE]; colnames(result) <- mapping$sample
      }
      return(result)
    }
    if (!is.numeric(value)) stop("Observation auxiliary must be numeric and named by its axis")
    .validate_ids(names(value), "Observation auxiliary")
    sample_match <- setequal(names(value), sample_ids); feature_match <- setequal(names(value), feature_ids)
    if (is.null(axis)) {
      if (sample_match == feature_match) stop("Ambiguous observation auxiliary axis; supply list(value=...,axis='sample'/'feature')")
      axis <- if (sample_match) "sample" else "feature"
    }
    if (axis == "feature" && feature_match) return(value[feature_ids])
    if (axis != "sample" || !sample_match) stop("Observation auxiliary IDs do not match declared axis")
    stats::setNames(value[match(mapping$original_sample, names(value))], mapping$sample)
  })
}

.checkpoint_replicate <- function(path, config_hash, compute) {
  if (is.null(path)) return(compute())
  if (file.exists(path)) {
    saved <- readRDS(path)
    if (!identical(saved$config_hash, config_hash)) stop("Checkpoint configuration/input mismatch; use a new checkpoint directory")
    return(saved$result)
  }
  result <- compute()
  temporary <- tempfile(pattern = ".replicate-", tmpdir = dirname(path))
  on.exit(unlink(temporary))
  saveRDS(list(config_hash = config_hash, result = result), temporary)
  if (!file.rename(temporary, path)) stop("Atomic checkpoint rename failed")
  result
}

.default_fit_callback <- function(X, metadata, sample_side, seed, settings, observation_aux = list()) {
  args <- settings$configuration
  if (!is.null(observation_aux$S)) args$S <- observation_aux$S
  args$verbose <- 0L
  do.call(fit_ebmf, c(list(X = X, sample_side = sample_side, seed = seed,
    ebnm_fn = settings$priors, feature_scale = settings$feature_scale,
    independent_unit = settings$unit_col), args))
}

#' Refitting stability under explicit biological-unit resampling
#'
#' Donor occurrences preserve all their visits/cells and get unique synthetic
#' IDs. Similarity quantiles summarize successful perturbations, not posterior
#' uncertainty. No threshold gives candidate similarities without recovery rates.
#' @param X Named input matrix in recorded native sample orientation.
#' @param reference_fit Native reference fit with wrapper metadata.
#' @param metadata Explicit sample and biological unit metadata.
#' @param unit_col Required donor/subject column.
#' @param B Positive replicate count.
#' @param resample Only unit resampling is supported.
#' @param strata Optional donor-level stratification columns.
#' @param fit_fn Callback(X, metadata, sample_side, seed, settings, observation_aux).
#' @param preprocess_fn Optional per-replicate callback returning X, metadata, manifest.
#' @param alignment Only native congruence assignment.
#' @param min_congruence Optional caller acceptance threshold.
#' @param seed Required deterministic master seed.
#' @param workers Positive worker count; Unix fork backend for workers greater than one.
#' @param checkpoint_dir Optional compatible atomic checkpoint directory.
#' @param keep_fits Retain large native fits only on explicit request.
#' @param control settings, observation_aux and explicit callback_configuration.
#' @return stability_result with summaries, replicates, matches, resamples, failures.
#' @export
#' @author Jason Willwerscheid; Peter Carbonetto; Wei Wang; Matthew Stephens.
#'   Full source and dependency credits: \url{https://mdmanurung.github.io/flashier-utils/articles/credits.html}.
#' @details Origin: NEW_ORCHESTRATION / call_only.
#'   Credits flash refits with existing congruence/assignment APIs; use provenance("factor_stability")
#'   for audited sources and runtime versions. Input units and basis are retained;
#'   display/conditional results do not establish biological replication.
#' @examplesIf requireNamespace("psych",quietly=TRUE) && requireNamespace("clue",quietly=TRUE)
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
#' # Two resamples illustrate execution, not reliable recovery estimates.
#' factor_stability(X,fit,metadata,"donor",B=2,seed=4,min_congruence=0.8)$summary
factor_stability <- function(X, reference_fit, metadata, unit_col, B = 100L, resample = "unit", strata = NULL,
                             fit_fn = NULL, preprocess_fn = NULL, alignment = "congruence", min_congruence = NULL,
                             seed, workers = 1L, checkpoint_dir = NULL, keep_fits = FALSE, control = list()) {
  if (missing(seed) || missing(unit_col)) stop("Choose seed and biological unit explicitly")
  if (resample != "unit" || alignment != "congruence") stop("Only biological-unit resampling and native congruence are supported")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B) || length(workers) != 1L || !is.finite(workers) || workers < 1 || workers != as.integer(workers)) stop("B and workers must be positive integers")
  if (workers > 1 && .Platform$OS.type != "unix") stop("Multiple workers require the supported Unix fork backend; use workers=1 here")
  if (!is.list(control) || any(!names(control) %in% c("settings", "observation_aux", "callback_configuration"))) stop("Unsupported stability control")
  X <- .validate_matrix(X)
  reference <- .resolve_view(reference_fit)
  side <- reference$manifest$sample_side
  samples <- if (side == "rows") rownames(X) else colnames(X)
  features <- if (side == "rows") colnames(X) else rownames(X)
  if (!setequal(samples, reference$manifest$sample_ids) || !setequal(features, reference$manifest$feature_ids)) stop("Input sample/feature IDs differ from the reference fit")
  md <- .match_ids(samples, metadata, required_columns = c("sample", unit_col, strata))$metadata
  custom <- !is.null(fit_fn)
  settings <- control$settings
  if (is.null(fit_fn)) {
    native_metadata <- attr(reference_fit, "flashbridge_metadata")
    prior_names <- native_metadata$prior_identifiers
    if (is.null(prior_names) || any(prior_names == "custom_callback_required")) stop("Reference custom/unknown priors require an explicit fit_fn")
    priors <- lapply(prior_names, function(name) getExportedValue("ebnm", strsplit(name, ";", fixed = TRUE)[[1L]][1L]))
    if (is.null(settings)) settings <- list(configuration = native_metadata$configuration, priors = priors, feature_scale = native_metadata$feature_scale, unit_col = unit_col)
    if (!is.null(settings$configuration$S) && length(settings$configuration$S) > 1L && is.null(control$observation_aux$S)) stop("Non-scalar S requires named observation_aux$S for donor resampling")
    fit_fn <- .default_fit_callback
  }
  if (!is.function(fit_fn) || (!is.null(preprocess_fn) && !is.function(preprocess_fn))) stop("Callbacks must be functions")
  if (!is.null(checkpoint_dir) && custom && is.null(control$callback_configuration)) stop("Custom callback checkpoints require callback_configuration including captured inputs")
  seeds <- .with_seed(seed, sample.int(.Machine$integer.max - 1L, B))
  callback_hash <- digest::digest(list(formals(fit_fn), body(fit_fn), if (is.null(preprocess_fn)) NULL else list(formals(preprocess_fn), body(preprocess_fn)), control$callback_configuration), algo = "sha256")
  config_hash <- digest::digest(list(X, md, reference$manifest$basis_id, settings, control$observation_aux, callback_hash,
    seed, B, strata, min_congruence, keep_fits, .engine_metadata(), body(.bootstrap_units), body(.align_programs)), algo = "sha256")
  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  run <- function(i) {
    path <- if (is.null(checkpoint_dir)) NULL else file.path(checkpoint_dir, paste0("replicate-", i, ".rds"))
    .checkpoint_replicate(path, config_hash, function() {
      resampled <- .bootstrap_units(md, unit_col, strata, seeds[i])
      mapping <- resampled$mapping
      Y <- if (side == "rows") X[match(mapping$original_sample, samples), , drop = FALSE] else X[, match(mapping$original_sample, samples), drop = FALSE]
      if (side == "rows") rownames(Y) <- mapping$sample else colnames(Y) <- mapping$sample
      auxiliary <- .resample_aux(if (is.null(control$observation_aux)) list() else control$observation_aux, mapping, X, side)
      started <- proc.time()[[3L]]
      result <- tryCatch({
        prepared <- if (is.null(preprocess_fn)) list(X = Y, metadata = resampled$metadata, manifest = list(scope = "conditional_on_supplied_preprocessing"), observation_aux = auxiliary) else
          preprocess_fn(Y, resampled$metadata, training_ids = mapping$sample, seed = seeds[i], settings = settings)
        if (!is.list(prepared) || !all(c("X", "metadata", "manifest") %in% names(prepared))) stop("Preprocessing must expose X, metadata and fitted manifest")
        fitted <- fit_fn(prepared$X, prepared$metadata, side, seeds[i], settings, observation_aux = if (is.null(prepared$observation_aux)) auxiliary else prepared$observation_aux)
        if (!inherits(fitted, "flash")) stop("fit_fn must return a native flash object")
        view <- standardize_factors(fitted, sample_side = side)
        if (!identical(view$manifest$sample_ids, as.character(prepared$metadata$sample))) stop("Refit callback did not preserve prepared sample IDs")
        alignment <- .align_programs(reference, view, min_congruence)
        list(status = "supported", matches = alignment$matches, pve = factor_pve(fitted), rank = fitted$n_factors,
             preparation = prepared$manifest, fit = if (keep_fits) fitted else NULL, alignment_metadata = alignment$analysis_metadata)
      }, error = function(e) list(status = "failed", reason = conditionMessage(e)))
      c(result, list(replicate = i, seed = seeds[i], resamples = mapping, elapsed_seconds = proc.time()[[3L]] - started))
    })
  }
  results <- if (workers == 1L) lapply(seq_len(B), run) else parallel::mclapply(seq_len(B), run, mc.cores = workers, mc.set.seed = FALSE)
  if (any(vapply(results, inherits, logical(1), "try-error"))) stop("Parallel worker/checkpoint failure; inspect checkpoint logs")
  successful <- vapply(results, function(r) identical(r$status, "supported"), logical(1))
  matches <- do.call(rbind, lapply(results[successful], function(r) cbind(replicate = r$replicate, r$matches)))
  resamples <- do.call(rbind, lapply(results, function(r) cbind(replicate = r$replicate, seed = r$seed, r$resamples)))
  failures <- data.frame(replicate = which(!successful), reason = vapply(results[!successful], function(r) r$reason, character(1)))
  summary <- do.call(rbind, lapply(reference$manifest$factor_ids, function(factor) {
    rows <- if (is.null(matches)) NULL else matches[matches$reference_factor == factor & !is.na(matches$reference_factor), , drop = FALSE]
    similarity <- if (is.null(rows)) numeric() else rows$similarity[is.finite(rows$similarity)]
    recovered <- if (is.null(rows)) 0L else sum(rows$matched %in% TRUE)
    eligible <- if (is.null(rows)) 0L else sum(rows$eligible)
    quantiles <- if (length(similarity)) stats::quantile(similarity, c(0.025, 0.5, 0.975), names = FALSE) else rep(NA_real_, 3)
    pv <- vapply(results[successful], function(r) {
      target <- r$matches$target_factor[match(factor, r$matches$reference_factor)]
      if (is.na(target)) NA_real_ else r$pve$pve[match(target, r$pve$factor)]
    }, numeric(1))
    data.frame(factor = factor, n_requested = B, n_successful = sum(successful), n_failed = sum(!successful), n_eligible = eligible,
      recovery_rate = if (is.null(min_congruence) || !eligible) NA_real_ else recovered / eligible,
      recovery_rate_all_requested = if (is.null(min_congruence)) NA_real_ else recovered / B,
      median_congruence = quantiles[2], q025_congruence = quantiles[1], q975_congruence = quantiles[3],
      native_pve_mean = if (any(is.finite(pv))) mean(pv, na.rm = TRUE) else NA_real_, native_pve_sd = if (sum(is.finite(pv)) > 1L) stats::sd(pv, na.rm = TRUE) else NA_real_,
      orientation_flip_rate = if (length(similarity)) mean(rows$sign[is.finite(rows$similarity)] == -1) else NA_real_)
  }))
  replicates <- do.call(rbind, lapply(results, function(r) data.frame(replicate = r$replicate, seed = r$seed, status = r$status, rank = if (is.null(r$rank)) NA_real_ else r$rank, elapsed_seconds = r$elapsed_seconds)))
  .attach_provenance(structure(list(schema_version = "1.0.0", summary = summary, replicates = replicates, matches = matches, resamples = resamples, failures = failures,
    fits = if (keep_fits) lapply(results, function(r) r$fit) else NULL,
    analysis_metadata = list(reference_basis = reference$manifest, config_hash = config_hash, callback_hash = callback_hash, seeds = seeds,
      unit_col = unit_col, strata = strata, n_requested = B, seed = seed, workers = workers, min_congruence = min_congruence,
      preparation = if (is.null(preprocess_fn)) "conditional_on_supplied_preprocessing" else "refit_per_replicate",
      eligibility = "successful_refit_with_nonzero_reference_and_at_least_one_nonzero_target", similarity_quantiles = "conditional_on_success_and_finite_assigned_similarity", checkpoint_dir = checkpoint_dir)), class = "stability_result"), "factor_stability", match.call(), .resolved_parameters("factor_stability", environment()))
}
