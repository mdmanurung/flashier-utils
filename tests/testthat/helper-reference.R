# Fixtures for the consensus-reference layer (flashier-starcat R003).
# Each fixture is a deterministic factor_representation built from fixed matrices,
# so tests need no fitting. `support` sets the declared activity/effect support.
reference_view <- function(
  effects,
  activity = NULL,
  support = c(activity = "signed", effects = "signed"),
  fit_id = "reference_fixture",
  feature_scale = "simulated_intensity",
  offset = NULL
) {
  if (is.null(activity)) {
    activity <- matrix(
      seq_len(3 * ncol(effects)) / 3,
      3,
      ncol(effects),
      dimnames = list(paste0(fit_id, "_s", 1:3), colnames(effects))
    )
  }
  manifest <- list(
    fit_id = fit_id,
    schema_version = "1.0.0",
    sample_side = "rows",
    native_indices = seq_len(ncol(effects)),
    scaling = "raw",
    orientation = "as_fit",
    constraints = support,
    feature_scale = feature_scale
  )
  if (is.null(offset)) {
    .make_representation(activity, effects, manifest)
  } else {
    .make_representation(activity, effects, manifest, offset)
  }
}

# M0 legacy-output snapshot (R001/R002). The generator script in
# validation/starcat-20261118/benchmark/ saved these outputs at b03ba4c.
legacy_baseline_outputs <- function() {
  scaled_fit <- function(K) {
    fit_ebmf(
      fixture_matrix(),
      "rows",
      ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace),
      max_factors = K,
      nullcheck = FALSE,
      backfit = TRUE,
      verbose = 0L,
      seed = 12,
      feature_scale = "simulated_intensity"
    )
  }
  fit <- scaled_fit(2L)
  fit3 <- scaled_fit(3L)
  view <- standardize_factors(fit, scaling = "raw", feature_scale = "simulated_intensity")
  priors <- list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)
  X <- fixture_matrix()
  replicability <- factor_replicability(view, standardize_factors(fit3, feature_scale = "simulated_intensity"))
  activity <- project_factor_activity(
    view,
    X[1:8, ] * 1.1,
    "rows",
    view$manifest$feature_scale,
    preprocessing = list(scope = "identity"),
    ebnm_fn = priors,
    maxiter = 50
  )
  fixed <- project_factor_activity(
    fit,
    X[1:6, ],
    "rows",
    view$manifest$feature_scale,
    preprocessing = list(scope = "identity"),
    prior_policy = "fixed_reference",
    maxiter = 50
  )
  set.seed(21)
  new_features <- tcrossprod(view$activity, matrix(stats::rnorm(6), 3, 2)) +
    matrix(stats::rnorm(36, sd = 0.2), 12, 3)
  dimnames(new_features) <- list(rownames(X), paste0("p", 1:3))
  features <- project_factor_features(view, new_features, "rows", "simulated_protein", ebnm_fn = priors, maxiter = 50)
  list(
    replicability = list(
      matches = replicability$matches,
      similarity = replicability$similarity_matrix,
      congruence = replicability$congruence_matrix
    ),
    activity = list(
      activity = activity$representation$activity,
      effects = activity$representation$effects,
      reconstruction = activity$reconstruction
    ),
    fixed_reference = list(activity = fixed$representation$activity, reconstruction = fixed$reconstruction),
    features = list(
      effects = features$representation$effects,
      offset = features$representation$offset,
      reconstruction = features$reconstruction
    )
  )
}

# Planted multi-dataset programs: `shared` programs appear in every dataset,
# each dataset also has `specific` private programs. Effects are features x
# programs; `noise` adds Gaussian perturbation to the effects. Seeded by caller.
planted_programs <- function(
  n_datasets = 3L,
  n_features = 40L,
  shared = 3L,
  specific = 1L,
  noise = 0,
  support = c(activity = "signed", effects = "signed"),
  scale = NULL
) {
  features <- sprintf("g%02d", seq_len(n_features))
  block <- function(k) {
    out <- numeric(n_features)
    width <- 4L
    out[((k - 1L) * width + seq_len(width) - 1L) %% n_features + 1L] <- c(4, 3, 2, 1)
    out
  }
  truth <- vapply(seq_len(shared), block, numeric(n_features))
  dimnames(truth) <- list(features, paste0("P", seq_len(shared)))
  views <- list()
  for (d in seq_len(n_datasets)) {
    private <- vapply(
      seq_len(specific),
      function(j) block(shared + (d - 1L) * specific + j),
      numeric(n_features)
    )
    effects <- cbind(truth, matrix(private, n_features))
    if (support[["effects"]] == "nonnegative") {
      effects <- abs(effects + noise * abs(matrix(stats::rnorm(length(effects)), n_features)))
    } else {
      effects <- effects + noise * matrix(stats::rnorm(length(effects)), n_features)
    }
    multiplier <- if (is.null(scale)) rep(1, ncol(effects)) else rep_len(scale, ncol(effects))
    effects <- sweep(effects, 2L, multiplier, "*")
    dimnames(effects) <- list(features, paste0("F", seq_len(ncol(effects))))
    views[[paste0("D", d)]] <- reference_view(effects, support = support, fit_id = paste0("D", d))
  }
  list(views = views, truth = truth, shared = shared, specific = specific)
}
