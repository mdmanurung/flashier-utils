# Native fixed-program projection on a NEW object; no same-dimension update shortcut.
.prediction_metrics <- function(observed, predicted) {
  keep <- is.finite(observed)
  o <- observed[keep]; p <- predicted[keep]
  if (!length(o) || any(!is.finite(p))) stop("No finite scored entries or nonfinite predictions")
  sst <- sum((o - mean(o))^2)
  c(rmse = sqrt(mean((o-p)^2)), mae = mean(abs(o-p)), prediction_r2 = if (sst > 0) 1-sum((o-p)^2)/sst else NA_real_,
    correlation = if (length(o)>1L && stats::sd(o)>0 && stats::sd(p)>0) stats::cor(o,p) else NA_real_, n_scored_entries = length(o))
}

.native_projection <- function(Y, training, sample_side, priors, maxiter = 100L, S = NULL, var_type = 0L) {
  K <- ncol(training$effects)
  if (!K) return(matrix(numeric(), nrow(Y), 0L, dimnames = list(rownames(Y), character())))
  A0 <- matrix(0.1, nrow(Y), K, dimnames = list(rownames(Y), colnames(training$effects)))
  data <- if (sample_side == "rows") Y else t(Y)
  S_dim <- NULL
  if (!is.null(S) && is.null(dim(S)) && length(S) > 1L) {
    row_match <- setequal(names(S), rownames(data)); column_match <- setequal(names(S), colnames(data))
    if (row_match == column_match) stop("Projection S vector has ambiguous axis IDs")
    S_dim <- if (row_match) 1L else 2L
  }
  object <- flashier::flash_init(data, S = S, var_type = var_type, S_dim = S_dim)
  pair <- if (sample_side == "rows") list(A0, training$effects) else list(training$effects, A0)
  object <- flashier::flash_factors_init(object, pair, ebnm_fn = priors)
  object <- flashier::flash_factors_fix(object, seq_len(K), which_dim = if (sample_side == "rows") "factors" else "loadings")
  object <- flashier::flash_backfit(object, maxiter = maxiter, extrapolate = FALSE, verbose = 0L)
  fixed <- if (sample_side == "rows") object$F_pm else object$L_pm
  if (!isTRUE(all.equal(unname(fixed), unname(training$effects), tolerance = 1e-12))) stop("Native projection changed fixed feature programs")
  A <- if (sample_side == "rows") object$L_pm else object$F_pm
  colnames(A) <- colnames(training$effects)
  A
}

#' Evaluate training programs on held-out biological units
#'
#' Projection scores use observations seen during activity estimation. Masked
#' scores use disjoint entries/features and never enter activity/noise fitting.
#' Learned hyperparameters may recalibrate on permitted test observations.
#' Paired preprocessing callbacks fit on training data then apply frozen state;
#' masked inference receives NA before transform. Global preprocessing is labelled
#' conditional and cannot establish independent full-pipeline validation.
#' @param X Named numerical input in explicitly declared native sample orientation.
#' @param metadata Sample and donor metadata.
#' @param unit_col Required biological unit column.
#' @param folds Data frame with sample or unit IDs and fold labels; donors cannot cross folds.
#' @param fit_fn Required native fitting callback following the shared plan contract.
#' @param preprocess_fn NULL or list(fit,transform); fit returns X,metadata,manifest,state,
#'   transform(X,metadata,state,seed,settings) returns transformed X and metadata.
#' @param evaluation projection or masked.
#' @param activity_features,scoring_features Optional disjoint named feature sets for masked mode.
#' @param mask Optional named logical matrix; TRUE entries are withheld for scoring.
#' @param freeze_noise TRUE is unsupported until a separate native freeze probe passes.
#' @param seed Required master seed.
#' @param workers Only serial fold execution is supported.
#' @param checkpoint_dir Fold checkpoints require explicit callback_configuration.
#' @param control sample_side (required), settings, observation_aux, preprocessing_manifest (required without callbacks),
#'   projection_ebnm_fn, maxiter, dense_budget_bytes, cohort_col, evaluated_folds,
#'   callback_configuration, keep_predictions, keep_fits.
#' @return holdout_result with metrics, folds, failures, activities and information boundaries.
#' @export
#' @details Origin: NEW_ORCHESTRATION / call_only.
#'   Credits flash_init; flash_factors_init; flash_factors_fix; flash_backfit; use provenance("factor_holdout")
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
#' folds <- data.frame(sample=metadata$sample,fold=rep(c("train","test"),c(16,8)))
#' callback <- function(X,metadata,sample_side,seed,settings,observation_aux=list()) {
#'   fit_ebmf(X,sample_side,max_factors=2,seed=seed,verbose=0,nullcheck=FALSE,
#'     feature_scale="simulated_centered_intensity")
#' }
#' control <- list(sample_side="rows",preprocessing_manifest=list(scope="identity"),evaluated_folds="test")
#' factor_holdout(X,metadata,"donor",folds,callback,evaluation="masked",seed=4,
#'   activity_features=paste0("g",1:8),scoring_features=paste0("g",9:16),control=control)$metrics
factor_holdout <- function(X, metadata, unit_col, folds, fit_fn, preprocess_fn = NULL,
                           evaluation = "projection", activity_features = NULL, scoring_features = NULL,
                           mask = NULL, freeze_noise = FALSE, seed, workers = 1L, checkpoint_dir = NULL, control = list()) {
  if (missing(seed) || missing(unit_col) || missing(fit_fn) || !is.function(fit_fn)) stop("Choose biological unit, fit_fn and seed explicitly")
  if (freeze_noise) stop("Freezing all trained noise/hyperparameters is unsupported; native recalibration uses permitted observations only")
  if (workers != 1L) stop("Only serial fold execution is supported")
  allowed <- c("sample_side", "settings", "observation_aux", "preprocessing_manifest", "projection_ebnm_fn", "maxiter", "dense_budget_bytes", "cohort_col", "evaluated_folds", "callback_configuration", "keep_predictions", "keep_fits")
  if (!is.list(control) || any(!names(control) %in% allowed)) stop("Unsupported holdout control")
  if (is.null(control$sample_side)) stop("Supply control$sample_side explicitly")
  side <- match.arg(control$sample_side, c("rows", "columns"))
  evaluation <- match.arg(evaluation, c("projection", "masked"))
  X <- .validate_matrix(X)
  samples <- if (side == "rows") rownames(X) else colnames(X)
  features <- if (side == "rows") colnames(X) else rownames(X)
  md <- .match_ids(samples, metadata, required_columns = c("sample", unit_col, control$cohort_col))$metadata
  units <- as.character(md[[unit_col]])
  if (any(!nzchar(trimws(units)))) stop("Biological unit IDs cannot be empty")
  if (!is.data.frame(folds) || !"fold" %in% names(folds)) stop("folds requires named sample/unit and fold columns")
  if ("sample" %in% names(folds)) labels <- .match_ids(samples, folds, required_columns = c("sample", "fold"))$metadata$fold else {
    if (!unit_col %in% names(folds)) stop("folds has no sample or declared unit IDs")
    ids <- .validate_ids(as.character(folds[[unit_col]]), "Fold unit")
    if (!setequal(ids, unique(units))) stop("Fold unit IDs differ from metadata")
    labels <- folds$fold[match(units, ids)]
  }
  labels <- as.character(labels)
  if (anyNA(labels) || any(!nzchar(labels)) || length(unique(labels)) < 2L) stop("At least two nonempty fold labels are required")
  for (unit in unique(units)) if (length(unique(labels[units == unit])) != 1L) stop("A donor's visits/cells cross folds")
  evaluated <- if (is.null(control$evaluated_folds)) unique(labels) else control$evaluated_folds
  .select_ids(evaluated, unique(labels), "fold")
  if (is.null(preprocess_fn)) {
    preprocessing <- control$preprocessing_manifest
    if (!is.list(preprocessing) || is.null(preprocessing$scope) || !preprocessing$scope %in% c("identity", "within_sample", "external_reference", "conditional_global")) stop("Without callbacks declare preprocessing_manifest$scope: identity, within_sample, external_reference or conditional_global")
    if (evaluation == "masked" && preprocessing$scope == "within_sample") stop("Masked within-sample transforms need raw data and callbacks that exclude withheld measurements")
  } else {
    if (!is.list(preprocess_fn) || !all(c("fit", "transform") %in% names(preprocess_fn)) || !all(vapply(preprocess_fn[c("fit", "transform")], is.function, logical(1)))) stop("preprocess_fn requires paired fit/transform callbacks")
    preprocessing <- list(scope = "training_only_callbacks")
  }
  if (!is.null(mask)) {
    if (!is.matrix(mask) || !is.logical(mask) || anyNA(mask) || !identical(dimnames(mask), dimnames(X)) || !identical(dim(mask), dim(X))) stop("mask must be a logical matrix with exactly X's named axes")
  }
  if (evaluation == "projection" && (!is.null(mask) || !is.null(activity_features) || !is.null(scoring_features))) stop("Masks and feature splits require explicit masked evaluation")
  if (evaluation == "masked" && is.null(mask) && (is.null(activity_features) || is.null(scoring_features))) stop("Masked evaluation requires a mask or disjoint activity/scoring feature IDs")
  if (!is.null(activity_features)) .select_ids(activity_features, features, "activity feature")
  if (!is.null(scoring_features)) .select_ids(scoring_features, features, "scoring feature")
  if (length(intersect(activity_features, scoring_features))) stop("Activity and scoring features must be disjoint")
  budget <- if (is.null(control$dense_budget_bytes)) 256 * 1024^2 else control$dense_budget_bytes
  if (length(budget) != 1L || !is.finite(budget) || budget <= 0) stop("dense_budget_bytes must be positive")
  if (!is.null(checkpoint_dir) && is.null(control$callback_configuration)) stop("Checkpointing callbacks requires explicit callback_configuration including captured inputs")
  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  config_hash <- digest::digest(list(X, md, folds, control, preprocessing, seed, evaluation, mask, activity_features, scoring_features,
    body(fit_fn), formals(fit_fn), lapply(preprocess_fn, function(fn) list(body(fn), formals(fn))), .engine_metadata()), algo = "sha256")
  fold_seeds <- .with_seed(seed, sample.int(.Machine$integer.max - 1L, length(evaluated)))
  results <- lapply(seq_along(evaluated), function(i) {
    label <- evaluated[i]
    path <- if (is.null(checkpoint_dir)) NULL else file.path(checkpoint_dir, paste0("fold-", i, ".rds"))
    .checkpoint_replicate(path, config_hash, function() tryCatch({
      train <- which(labels != label); test <- which(labels == label)
      subset_samples <- function(index) if (side == "rows") X[index, , drop = FALSE] else X[, index, drop = FALSE]
      trainX <- subset_samples(train); testX <- subset_samples(test)
      auxiliary <- if (is.null(control$observation_aux)) list() else control$observation_aux
      train_aux <- .resample_aux(auxiliary, data.frame(original_sample = samples[train], sample = samples[train]), X, side)
      test_aux <- .resample_aux(auxiliary, data.frame(original_sample = samples[test], sample = samples[test]), X, side)
      if (8 * length(test) * length(features) * 4 > budget) stop("Held-out dense inference/scoring copies exceed the declared memory budget")
      truth <- if (side == "rows") as.matrix(testX) else t(as.matrix(testX))
      inference <- truth
      score <- is.finite(truth)
      if (evaluation == "masked") {
        masked <- if (is.null(mask)) matrix(FALSE, nrow(truth), ncol(truth)) else if (side == "rows") mask[test, , drop = FALSE] else t(mask[, test, drop = FALSE])
        if (!is.null(activity_features)) masked[, !colnames(truth) %in% activity_features] <- TRUE
        if (!is.null(scoring_features)) { score[, !colnames(truth) %in% scoring_features] <- FALSE; masked[, colnames(truth) %in% scoring_features] <- TRUE }
        score <- score & masked
        inference[masked] <- NA_real_
        if (!any(score) || any(rowSums(is.finite(inference)) == 0L)) stop("A masked fold has no scoring entries or no permitted activity observations for a sample")
      }
      inferenceX <- if (side == "rows") inference else t(inference)
      prepared <- if (is.null(preprocess_fn)) list(X = trainX, metadata = md[train, , drop = FALSE], manifest = preprocessing, observation_aux = train_aux) else
        preprocess_fn$fit(trainX, md[train, , drop = FALSE], training_ids = samples[train], seed = fold_seeds[i], settings = control$settings)
      if (!is.list(prepared) || !all(c("X", "metadata", "manifest") %in% names(prepared))) stop("Training preprocessing must expose X, metadata and fitted manifest")
      transformed <- if (is.null(preprocess_fn)) list(X = inferenceX, metadata = md[test, , drop = FALSE]) else
        preprocess_fn$transform(inferenceX, md[test, , drop = FALSE], state = prepared$state, seed = fold_seeds[i], settings = control$settings)
      train_fit <- fit_fn(prepared$X, prepared$metadata, side, fold_seeds[i], control$settings, observation_aux = if (is.null(prepared$observation_aux)) train_aux else prepared$observation_aux)
      view <- standardize_factors(train_fit, sample_side = side, scaling = "raw")
      if (!identical(view$manifest$sample_ids, samples[train]) || !identical(view$manifest$feature_ids, features)) stop("Training callback changed required sample/feature IDs")
      Y <- if (side == "rows") as.matrix(transformed$X) else t(as.matrix(transformed$X))
      if (!identical(dimnames(Y), dimnames(truth))) stop("Held-out transform changed sample/feature axes")
      if (evaluation == "masked" && any(is.finite(Y[score]))) stop("Held-out preprocessing filled scoring entries; masked information boundary violated")
      Y <- sweep(Y, 2L, view$offset, "-")
      priors <- control$projection_ebnm_fn
      if (is.null(priors) && ncol(view$effects)) {
        names <- attr(train_fit, "flashbridge_metadata")$prior_identifiers
        if (is.null(names) || any(names == "custom_callback_required")) stop("Projection requires explicit native prior functions for a custom/bare callback")
        priors <- lapply(names, function(name) getExportedValue("ebnm", strsplit(name, ";", fixed = TRUE)[[1L]][1L]))
      }
      A <- .native_projection(Y, view, side, priors, if (is.null(control$maxiter)) 100L else control$maxiter, S = test_aux$S, var_type = if (is.null(attr(train_fit, "flashbridge_metadata")$configuration$var_type)) 0L else attr(train_fit, "flashbridge_metadata")$configuration$var_type)
      prediction <- sweep(tcrossprod(A, view$effects), 2L, view$offset, "+")
      # Score transformation occurs after activity inference using the same frozen training state.
      if (!is.null(preprocess_fn)) {
        scoring <- preprocess_fn$transform(testX, md[test, , drop = FALSE], state = prepared$state, seed = fold_seeds[i], settings = control$settings)
        truth <- if (side == "rows") as.matrix(scoring$X) else t(as.matrix(scoring$X))
      }
      trainingY <- if (side == "rows") prepared$X else Matrix::t(prepared$X)
      baseline <- matrix(Matrix::colMeans(trainingY, na.rm = TRUE), length(test), length(features), byrow = TRUE)
      observed <- truth; observed[!score] <- NA_real_
      metric_rows <- list()
      score_rows <- function(indices, weighting, unit = NA_character_, cohort = NA_character_) {
        for (model in c("fixed_programs", "training_mean")) {
          predicted <- if (model == "fixed_programs") prediction else baseline
          values <- .prediction_metrics(observed[indices, , drop = FALSE], predicted[indices, , drop = FALSE])
          metric_rows[[length(metric_rows)+1L]] <<- data.frame(fold = label, evaluation_type = evaluation, model = model, weighting = weighting, unit = unit, cohort = cohort,
            rmse = values["rmse"], mae = values["mae"], prediction_r2 = values["prediction_r2"], correlation = values["correlation"], n_scored_entries = values["n_scored_entries"],
            n_samples = length(indices), n_units = length(unique(units[test][indices])), n_features = sum(colSums(score[indices, , drop = FALSE]) > 0))
        }
      }
      score_rows(seq_along(test), "observation")
      for (unit in unique(units[test])) score_rows(which(units[test] == unit), "unit", unit = unit)
      if (!is.null(control$cohort_col)) for (cohort in unique(as.character(md[[control$cohort_col]][test]))) score_rows(which(as.character(md[[control$cohort_col]][test]) == cohort), "cohort", cohort = cohort)
      metrics <- do.call(rbind, metric_rows)
      for (weighting in c("unit", if (!is.null(control$cohort_col)) "cohort")) for (model in c("fixed_programs", "training_mean")) {
        rows <- metrics$weighting == weighting & metrics$model == model
        macro <- metrics[which(rows)[1], , drop = FALSE]
        for (name in c("rmse", "mae", "prediction_r2", "correlation")) macro[[name]] <- if (all(is.na(metrics[[name]][rows]))) NA_real_ else mean(metrics[[name]][rows], na.rm = TRUE)
        macro$weighting <- paste0("macro_", weighting); macro$unit <- macro$cohort <- NA_character_
        macro$n_samples <- length(test); macro$n_units <- length(unique(units[test])); macro$n_scored_entries <- sum(metrics$n_scored_entries[rows])
        metrics <- rbind(metrics, macro)
      }
      attr(A, "analysis_metadata") <- c(view$manifest, list(training_basis_id = view$manifest$basis_id, heldout_sample_ids = rownames(A), program_values_fixed = TRUE, prior_noise_policy = "recalibrate_on_permitted_observed_entries"))
      list(status = "supported", fold = label, metrics = metrics, activity = A, effects = view$effects, preparation = prepared$manifest,
           training_ids = samples[train], scoring_ids = samples[test], train_basis = view$manifest,
           prediction = if (isTRUE(control$keep_predictions)) prediction else NULL, fit = if (isTRUE(control$keep_fits)) train_fit else NULL)
    }, error = function(e) list(status = "failed", fold = label, reason = conditionMessage(e))))
  })
  names(results) <- evaluated
  successful <- vapply(results, function(result) result$status == "supported", logical(1))
  failures <- data.frame(fold = evaluated[!successful], reason = vapply(results[!successful], function(result) result$reason, character(1)))
  .attach_provenance(structure(list(schema_version = "1.0.0", metrics = do.call(rbind, lapply(results[successful], function(result) result$metrics)),
    folds = data.frame(sample = samples, unit = units, fold = labels), failures = failures,
    feature_coverage = data.frame(feature = features, used_for_activity = if (is.null(activity_features)) NA else features %in% activity_features, scored = if (is.null(scoring_features)) NA else features %in% scoring_features),
    diagnostics = list(failures = failures, prior_noise_policy = "recalibrated_on_permitted_test_entries", auxiliary_policy = "named_sample_feature_subset"),
    activities = lapply(results[successful], function(result) result$activity), fold_results = results,
    analysis_metadata = list(evaluation_type = evaluation, preprocessing = preprocessing, sample_side = side, unit_col = unit_col, seed = seed,
      config_hash = config_hash, program_values_fixed = TRUE, hyperparameters = "recalibrated_on_permitted_test_entries", freeze_noise = FALSE,
      dense_budget_bytes = budget, projection_is_unseen_entry_prediction = FALSE)), class = "holdout_result"), "factor_holdout", match.call(), .resolved_parameters("factor_holdout", environment()))
}
