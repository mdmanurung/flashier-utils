# FU-3A — RNA->ADT fixed-activity projection validation (native oracle,
# held-out-cell design). Mathematical configuration follows the ImmGenT-GP
# CITE-seq method: Y_ADT ~ A_RNA B_ADT' with A_RNA fixed, exact measured cell
# intersection, signed protein priors, feature variance estimated on the query.
# No exact numeric parity with ImmGenT is claimed (different backfit controls /
# iteration schedule); the oracle is the native fixed-activity projection itself.

# Shared tiny oracle: strong 2-factor RNA fit, ADT synthesized from the FITTED
# activities so the oracle is exact up to noise.
fu3a_fixture <- function(noise_sd = 0.05, n_cells = 60, seed = 7) {
  set.seed(seed)
  A_true <- matrix(rnorm(n_cells * 2, sd = c(1.5, 1)), n_cells, 2,
                   dimnames = list(paste0("cell", seq_len(n_cells)), c("F1", "F2")))
  B_rna <- matrix(0, 30, 2, dimnames = list(paste0("gene", 1:30), c("F1", "F2")))
  B_rna[1:10, 1] <- 1.5
  B_rna[11:20, 2] <- 2
  Y_rna <- tcrossprod(A_true, B_rna) + matrix(rnorm(n_cells * 30, sd = 0.1), n_cells, 30)
  fit <- fit_ebmf(Y_rna, "rows", max_factors = 2, seed = 3, verbose = 0,
                  ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace),
                  nullcheck = FALSE, backfit = TRUE, feature_scale = "simulated_rna")
  view <- standardize_factors(fit)
  A_rna <- factor_activity(view)
  K <- ncol(A_rna)
  # signed protein effects on the fitted basis (2 proteins per factor pattern)
  B_adt <- matrix(0, 4, K, dimnames = list(paste0("ADT", 1:4), colnames(A_rna)))
  B_adt[1, 1] <- 2.0; B_adt[2, 1] <- -1.0
  if (K >= 2) { B_adt[3, 2] <- 1.5; B_adt[4, 2] <- -0.8 }
  Y_adt <- tcrossprod(A_rna, B_adt) +
    matrix(rnorm(n_cells * 4, sd = noise_sd), n_cells, 4,
           dimnames = list(rownames(A_rna), rownames(B_adt)))
  list(view = view, A = A_rna, B_adt = B_adt, Y_adt = Y_adt)
}

test_that("fixed-activity projection recovers oracle ADT effects", {
  fx <- fu3a_fixture()
  proj <- project_factor_features(fx$view, fx$Y_adt, sample_side = "rows",
                                  feature_scale = "clr_adt_sim",
                                  ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace))
  B_est <- factor_effects(proj$representation)
  expect_identical(rownames(B_est), rownames(fx$B_adt))
  # per-protein recovery against the oracle
  for (p in rownames(fx$B_adt)) {
    expect_gt(cor(B_est[p, ], fx$B_adt[p, ]), 0.98)
  }
  expect_lt(sqrt(mean((B_est - fx$B_adt)^2)), 0.15)
  # fixed activities retained exactly (values; metadata differs by design)
  expect_equal(unname(factor_activity(proj$representation)), unname(fx$A), ignore_attr = TRUE)
  # in-sample reconstruction metrics are reported
  expect_true(all(c("rmse", "fit_r2") %in% names(proj$reconstruction)))
})

test_that("held-out cells are reconstructable from training-fit ADT effects", {
  fx <- fu3a_fixture()
  train <- paste0("cell", 1:45); heldout <- paste0("cell", 46:60)
  proj <- project_factor_features(fx$view, fx$Y_adt[train, ], sample_side = "rows",
                                  feature_scale = "clr_adt_sim",
                                  ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace))
  B_tr <- factor_effects(proj$representation)
  # held-out prediction uses ONLY their externally available RNA activities;
  # their ADT values were never observed by the fit
  Y_hat <- tcrossprod(fx$A[heldout, ], B_tr)
  for (p in rownames(B_tr)) {
    expect_gt(cor(Y_hat[, p], fx$Y_adt[heldout, p]), 0.85)
  }
})

test_that("rank-deficient fixed activity is refused explicitly", {
  fx <- fu3a_fixture()
  A <- fx$A
  A_rd <- cbind(A, Fdup = A[, 1])  # exact duplicate column
  B <- factor_effects(fx$view)
  B_rd <- cbind(B, Fdup = 0)
  v_rd <- .make_representation(A_rd, B_rd, list(fit_id = "rankdef", schema_version = "1.0.0",
        sample_side = "rows", native_indices = seq_len(ncol(A_rd)), scaling = "raw",
        orientation = "as_fit", constraints = c(activity = "signed", effects = "signed"),
        feature_scale = "simulated_rna"))
  expect_error(
    project_factor_features(v_rd, fx$Y_adt, sample_side = "rows", feature_scale = "clr_adt_sim",
                            ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)),
    "Rank-deficient"
  )
})

test_that("cells absent from the reference activities are refused", {
  fx <- fu3a_fixture()
  Y_bad <- fx$Y_adt
  rownames(Y_bad)[1] <- "unknown_cell"
  expect_error(
    project_factor_features(fx$view, Y_bad, sample_side = "rows", feature_scale = "clr_adt_sim",
                            ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)),
    "absent from reference activities"
  )
})

test_that("undeclared or unknown ADT units are refused", {
  fx <- fu3a_fixture()
  expect_error(
    project_factor_features(fx$view, fx$Y_adt, sample_side = "rows", feature_scale = "unknown",
                            ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)),
    "Declare non-unknown measurement units"
  )
})
