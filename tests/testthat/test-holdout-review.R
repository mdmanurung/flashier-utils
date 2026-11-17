# Statistical review of factor_holdout() and the evidence pairing in compare_prior_fits().
review_callback <- function(X, metadata, sample_side, seed, settings, observation_aux = list()) {
  priors <- list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)
  if (sample_side == "columns") {
    priors <- rev(priors)
  }
  fit_ebmf(
    X,
    sample_side,
    ebnm_fn = priors,
    max_factors = 2L,
    backfit = TRUE,
    nullcheck = FALSE,
    seed = seed,
    verbose = 0,
    feature_scale = "intensity",
    var_type = if (is.null(settings$var_type)) 0L else settings$var_type
  )
}

review_setup <- function(n_test = 4L, per_donor = 2L) {
  X <- fixture_matrix()
  md <- data.frame(
    sample = rownames(X),
    donor = paste0("d", rep(seq_len(12 / per_donor), each = per_donor)),
    cohort = rep(c("A", "B"), each = 6)
  )
  folds <- data.frame(sample = md$sample, fold = c(rep("train", 12 - n_test), rep("test", n_test)))
  list(X = X, md = md, folds = folds)
}

review_control <- function(...) {
  utils::modifyList(
    list(sample_side = "rows", preprocessing_manifest = list(scope = "identity"), evaluated_folds = "test"),
    list(...)
  )
}

test_that("metrics, the training-mean baseline and skill equal a hand computation", {
  s <- review_setup()
  out <- factor_holdout(
    s$X,
    s$md,
    "donor",
    s$folds,
    review_callback,
    seed = 12,
    control = review_control(keep_predictions = TRUE)
  )
  prediction <- out$fold_results$test$prediction
  truth <- s$X[9:12, ]
  baseline <- matrix(colMeans(s$X[1:8, ]), 4, 8, byrow = TRUE)
  row_of <- function(model) out$metrics[out$metrics$weighting == "observation" & out$metrics$model == model, ]
  fit_row <- row_of("fixed_programs")
  base_row <- row_of("training_mean")
  sse <- sum((truth - prediction)^2)
  sst <- sum((truth - mean(truth))^2)
  sse_base <- sum((truth - baseline)^2)
  expect_equal(fit_row$rmse, sqrt(mean((truth - prediction)^2)))
  expect_equal(fit_row$mae, mean(abs(truth - prediction)))
  expect_equal(fit_row$prediction_r2, 1 - sse / sst)
  expect_equal(fit_row$skill, 1 - sse / sse_base)
  expect_equal(fit_row$correlation, cor(as.vector(truth), as.vector(prediction)))
  expect_equal(base_row$rmse, sqrt(mean((truth - baseline)^2)))
  expect_equal(base_row$prediction_r2, 1 - sse_base / sst)
  expect_true(is.na(base_row$skill))
  # prediction_r2 is centred on the held-out grand mean, so it can be high for the baseline too
  expect_gt(base_row$prediction_r2, 0)
})

test_that("masked scoring entries never reach activity estimation", {
  s <- review_setup()
  run <- function(X) {
    factor_holdout(
      X,
      s$md,
      "donor",
      s$folds,
      review_callback,
      seed = 12,
      evaluation = "masked",
      activity_features = paste0("g", 1:4),
      scoring_features = paste0("g", 5:8),
      control = review_control()
    )
  }
  base <- run(s$X)
  altered <- s$X
  altered[9:12, paste0("g", 5:8)] <- altered[9:12, paste0("g", 5:8)] + 100
  changed <- run(altered)
  # the activity values are identical; only the fit-hash attributes differ with the altered input
  expect_equal(changed$activities$test, base$activities$test, ignore_attr = TRUE)
  expect_gt(
    changed$metrics$rmse[changed$metrics$model == "fixed_programs"][1],
    base$metrics$rmse[base$metrics$model == "fixed_programs"][1]
  )
})

test_that("holdout metrics do not depend on the row order of the inputs", {
  s <- review_setup()
  order <- c(12:1)
  shuffled <- factor_holdout(
    s$X[order, ],
    s$md[order, ],
    "donor",
    s$folds[order, ],
    review_callback,
    seed = 12,
    control = review_control()
  )
  base <- factor_holdout(s$X, s$md, "donor", s$folds, review_callback, seed = 12, control = review_control())
  key <- function(m) m[order(m$weighting, m$unit, m$model), c("rmse", "mae", "skill")]
  expect_equal(key(shuffled$metrics), key(base$metrics), ignore_attr = TRUE)
})

test_that("a one-sample fold is refused up front when samples are rows and works for columns", {
  s <- review_setup(n_test = 1L, per_donor = 1L)
  expect_error(
    factor_holdout(s$X, s$md, "donor", s$folds, review_callback, seed = 12, control = review_control()),
    "one sample cannot be projected"
  )
  out <- factor_holdout(
    t(s$X),
    s$md,
    "donor",
    s$folds,
    review_callback,
    seed = 12,
    control = review_control(sample_side = "columns")
  )
  expect_equal(nrow(out$failures), 0L)
})

test_that("a masked feature split runs when variance is estimated per feature or per sample", {
  set.seed(1)
  X <- fixture_matrix()
  X <- rbind(X, X + matrix(rnorm(length(X), sd = 0.3), nrow(X)))
  rownames(X) <- paste0("s", 1:24)
  md <- data.frame(sample = rownames(X), donor = paste0("d", 1:24))
  folds <- data.frame(sample = md$sample, fold = rep(c("train", "test"), c(18, 6)))
  for (var_type in list(1L, 2L)) {
    out <- allow_maxiter(factor_holdout(
      X,
      md,
      "donor",
      folds,
      review_callback,
      seed = 12,
      evaluation = "masked",
      activity_features = paste0("g", 1:4),
      scoring_features = paste0("g", 5:8),
      control = review_control(settings = list(var_type = var_type))
    ))
    expect_equal(nrow(out$failures), 0L)
    expect_true(all(is.finite(out$activities$test)))
  }
})

test_that("a unit with no scored entry under a custom mask has no metric row instead of failing the fold", {
  s <- review_setup()
  mask <- matrix(FALSE, nrow(s$X), ncol(s$X), dimnames = dimnames(s$X))
  mask[9:10, 1:3] <- TRUE
  out <- factor_holdout(
    s$X,
    s$md,
    "donor",
    s$folds,
    review_callback,
    seed = 12,
    evaluation = "masked",
    mask = mask,
    control = review_control()
  )
  expect_equal(nrow(out$failures), 0L)
  units <- out$metrics$unit[out$metrics$weighting == "unit"]
  expect_true("d5" %in% units)
  expect_false("d6" %in% units)
})

test_that("a failed fold is not checkpointed, so a rerun after the fix succeeds", {
  s <- review_setup()
  state <- new.env()
  state$fail <- TRUE
  callback <- function(X, metadata, sample_side, seed, settings, observation_aux = list()) {
    if (state$fail) {
      stop("transient failure")
    }
    review_callback(X, metadata, sample_side, seed, settings, observation_aux)
  }
  dir <- tempfile("ckpt")
  control <- review_control(callback_configuration = list(version = 1))
  first <- factor_holdout(s$X, s$md, "donor", s$folds, callback, seed = 12, checkpoint_dir = dir, control = control)
  expect_equal(nrow(first$failures), 1L)
  state$fail <- FALSE
  second <- factor_holdout(s$X, s$md, "donor", s$folds, callback, seed = 12, checkpoint_dir = dir, control = control)
  expect_equal(nrow(second$failures), 0L)
  expect_true("fixed_programs" %in% second$metrics$model)
})

test_that(".checkpoint_replicate stores successes and replays them but never stores failures", {
  path <- tempfile(fileext = ".rds")
  calls <- 0L
  failing <- function() {
    calls <<- calls + 1L
    list(status = "failed", reason = "x")
  }
  .checkpoint_replicate(path, "h", failing)
  .checkpoint_replicate(path, "h", failing)
  expect_equal(calls, 2L)
  expect_false(file.exists(path))
  ok <- function() {
    calls <<- calls + 1L
    list(status = "supported", value = 1)
  }
  .checkpoint_replicate(path, "h", ok)
  .checkpoint_replicate(path, "h", ok)
  expect_equal(calls, 3L)
})

test_that("holdout evidence is paired only for the same scored entries, folds and no failures", {
  B <- diag(2)
  dimnames(B) <- list(c("g1", "g2"), c("F1", "F2"))
  r <- list(
    activity = matrix(1:6, 3, 2, dimnames = list(paste0("s", 1:3), c("F1", "F2"))),
    effects = B,
    support = list(activity = "signed", effects = "signed"),
    feature_scale = "simulated",
    preprocessing = list(scope = "fixture", transformation = "none", compatibility = "fixture")
  )
  s <- .program_source(r)
  stub <- function(truth = "t1", failures = 0L, folds = "test", factors = 2L) {
    structure(
      list(
        metrics = data.frame(fold = folds, weighting = "observation", unit = NA_character_, cohort = NA_character_),
        failures = data.frame(fold = rep("x", failures), reason = rep("boom", failures)),
        analysis_metadata = list(
          input_digest = "input",
          evaluation_digest = "same",
          truth_digest = truth,
          evaluation_type = "masked",
          fold_factor_counts = c(test = factors)
        )
      ),
      class = "holdout_result"
    )
  }
  evidence <- function(a, b) {
    list(
      reference = list(holdout = a, source_basis_id = s$basis_id),
      target = list(holdout = b, source_basis_id = s$basis_id)
    )
  }
  sources <- list(reference = s, target = s)
  expect_equal(.program_evidence(evidence(stub(), stub()), sources)$state, "driver_bound_paired_evaluation_attached")
  expect_error(.program_evidence(evidence(stub("t1"), stub("t2")), sources), "paired")
  expect_error(.program_evidence(evidence(stub(failures = 1L), stub()), sources), "failed folds")
  expect_error(.program_evidence(evidence(stub(), stub(folds = "other")), sources), "different folds")
  expect_warning(.program_evidence(evidence(stub(factors = 5L), stub()), sources), "may not belong to this source")
})

test_that("a masked split with var_type 0 still works when samples are columns and one feature is observed", {
  s <- review_setup()
  out <- factor_holdout(
    t(s$X),
    s$md,
    "donor",
    s$folds,
    review_callback,
    seed = 12,
    evaluation = "masked",
    activity_features = "g1",
    scoring_features = paste0("g", 5:8),
    control = review_control(sample_side = "columns")
  )
  expect_equal(nrow(out$failures), 0L)
})

test_that("a masked feature split with a supplied S runs for var_type 1 and 2, samples as rows or columns", {
  set.seed(1)
  X <- fixture_matrix()
  X <- rbind(X, X + matrix(rnorm(length(X), sd = 0.3), nrow(X)))
  rownames(X) <- paste0("s", 1:24)
  md <- data.frame(sample = rownames(X), donor = paste0("d", 1:24))
  folds <- data.frame(sample = md$sample, fold = rep(c("train", "test"), c(18, 6)))
  S_callback <- function(X, metadata, sample_side, seed, settings, observation_aux = list()) {
    fit_ebmf(
      X,
      sample_side,
      ebnm_fn = if (sample_side == "rows") {
        list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)
      } else {
        list(ebnm::ebnm_point_laplace, ebnm::ebnm_normal)
      },
      max_factors = 2L,
      backfit = TRUE,
      nullcheck = FALSE,
      seed = seed,
      verbose = 0,
      feature_scale = "intensity",
      var_type = settings$var_type,
      S = observation_aux$S
    )
  }
  features <- colnames(X)
  S_vector <- stats::setNames(seq(0.2, 0.4, length.out = length(features)), features)
  for (side in c("rows", "columns")) {
    data <- if (side == "rows") X else t(X)
    for (S in list(0.3, S_vector)) {
      for (var_type in list(1L, 2L)) {
        out <- allow_maxiter(factor_holdout(
          data,
          md,
          "donor",
          folds,
          S_callback,
          seed = 12,
          evaluation = "masked",
          activity_features = paste0("g", 1:4),
          scoring_features = paste0("g", 5:8),
          control = review_control(
            sample_side = side,
            settings = list(var_type = var_type),
            observation_aux = list(S = S)
          )
        ))
        expect_equal(nrow(out$failures), 0L)
        expect_true(all(is.finite(out$activities$test)))
      }
    }
  }
})

test_that("a masked split with a supplied S and var_type 0 is unchanged by the S feature drop", {
  # 18 training samples against 8 features keeps the S axis unambiguous in the training fit
  set.seed(1)
  X <- fixture_matrix()
  X <- rbind(X, X + matrix(rnorm(length(X), sd = 0.3), nrow(X)))
  rownames(X) <- paste0("s", 1:24)
  md <- data.frame(sample = rownames(X), donor = paste0("d", 1:24))
  folds <- data.frame(sample = md$sample, fold = rep(c("train", "test"), c(18, 6)))
  S_vector <- stats::setNames(rep(0.3, ncol(X)), colnames(X))
  S_callback <- function(X, metadata, sample_side, seed, settings, observation_aux = list()) {
    fit_ebmf(
      X,
      sample_side,
      ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace),
      max_factors = 2L,
      backfit = TRUE,
      nullcheck = FALSE,
      seed = seed,
      verbose = 0,
      feature_scale = "intensity",
      var_type = 0L,
      S = observation_aux$S
    )
  }
  run <- function(control) {
    allow_maxiter(factor_holdout(
      X,
      md,
      "donor",
      folds,
      S_callback,
      seed = 12,
      evaluation = "masked",
      activity_features = paste0("g", 1:4),
      scoring_features = paste0("g", 5:8),
      control = control
    ))
  }
  out <- run(review_control(observation_aux = list(S = S_vector)))
  expect_equal(nrow(out$failures), 0L)
  expect_true(all(is.finite(out$activities$test)))
})

test_that("a masked feature split with a matrix S runs for per-feature and per-sample variance", {
  set.seed(1)
  X <- fixture_matrix()
  X <- rbind(X, X + matrix(rnorm(length(X), sd = 0.3), nrow(X)))
  rownames(X) <- paste0("s", 1:24)
  md <- data.frame(sample = rownames(X), donor = paste0("d", 1:24))
  folds <- data.frame(sample = md$sample, fold = rep(c("train", "test"), c(18, 6)))
  S_callback <- function(X, metadata, sample_side, seed, settings, observation_aux = list()) {
    fit_ebmf(
      X,
      sample_side,
      ebnm_fn = if (sample_side == "rows") {
        list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)
      } else {
        list(ebnm::ebnm_point_laplace, ebnm::ebnm_normal)
      },
      max_factors = 2L,
      backfit = TRUE,
      nullcheck = FALSE,
      seed = seed,
      verbose = 0,
      feature_scale = "intensity",
      var_type = settings$var_type,
      S = observation_aux$S
    )
  }
  for (side in c("rows", "columns")) {
    data <- if (side == "rows") X else t(X)
    # a matrix S in the orientation of the data, with a little variation so it is not a scalar in disguise
    S_matrix <- matrix(0.3 + 0.05 * seq_along(data) / length(data), nrow(data), ncol(data), dimnames = dimnames(data))
    for (var_type in list(1L, 2L)) {
      out <- allow_maxiter(factor_holdout(
        data,
        md,
        "donor",
        folds,
        S_callback,
        seed = 12,
        evaluation = "masked",
        activity_features = paste0("g", 1:4),
        scoring_features = paste0("g", 5:8),
        control = review_control(
          sample_side = side,
          settings = list(var_type = var_type),
          observation_aux = list(S = S_matrix)
        )
      ))
      expect_equal(nrow(out$failures), 0L, info = paste(side, var_type))
      expect_true(all(is.finite(out$activities$test)), info = paste(side, var_type))
    }
  }
})

test_that("with var_type 0 the masked activities equal a direct projection with no feature dropped", {
  s <- review_setup()
  out <- factor_holdout(
    s$X,
    s$md,
    "donor",
    s$folds,
    review_callback,
    seed = 12,
    evaluation = "masked",
    activity_features = paste0("g", 1:4),
    scoring_features = paste0("g", 5:8),
    control = review_control(keep_fits = TRUE)
  )
  fit <- out$fold_results$test$fit
  view <- standardize_factors(fit, sample_side = "rows", scaling = "raw")
  # the old path: every feature, scoring features as NA, offset removed, no drop
  Y <- s$X[9:12, ]
  Y[, paste0("g", 5:8)] <- NA_real_
  Y <- sweep(Y, 2L, view$offset, "-")
  priors <- lapply(attr(fit, "flashbridge_metadata")$prior_identifiers, .prior_from_identifier)
  direct <- .native_projection(Y, view, "rows", priors, 100L, var_type = 0L)
  expect_equal(out$activities$test, direct, ignore_attr = TRUE)
})
