# flashier-starcat M1-M6: consensus reference construction and projection.
generic_manifest <- function(scale = "simulated_intensity") {
  list(mode = "generic", feature_scale = scale, preprocessing_compatibility = "user_declared")
}
projection_manifest <- function(prior = "ebnm_normal") {
  list(
    feature_scale = "simulated_intensity",
    scope = "external_reference",
    transformation = "identity",
    compatibility = "exact_declared_contract",
    activity_prior = prior
  )
}
build_generic <- function(views, ...) {
  build_factor_reference(views, comparison_manifest = generic_manifest(), min_similarity = 0.8, ...)
}
graph_members <- function(ids) {
  data.frame(member_id = ids, dataset_id = sub("[0-9]+$", "", ids), eligible = TRUE)
}
graph_similarity <- function(i, j, s) {
  data.frame(
    member_i = i,
    member_j = j,
    dataset_i = sub("[0-9]+$", "", i),
    dataset_j = sub("[0-9]+$", "", j),
    similarity = s,
    eligible = TRUE
  )
}
cluster_of <- function(i, j, s, rule = "merged_set", density = 2 / 3, reciprocal = FALSE, best = FALSE, rank = NULL) {
  members <- graph_members(sort(unique(c(i, j)), method = "radix"))
  graph <- .reference_build_graph(graph_similarity(i, j, s), members, 0.5, rank, reciprocal, best)
  out <- .reference_cluster_graph(graph, members, rule, density)
  split(out$membership$member_id, out$membership$cluster_id)
}

test_that("sources are immutable and every invariance holds", {
  set.seed(3)
  p <- planted_programs(noise = 0.05)
  before <- serialize(p$views, NULL)
  ref <- build_generic(p$views)
  expect_identical(serialize(p$views, NULL), before)
  expect_s3_class(ref, "factor_reference")
  expect_identical(build_generic(rev(p$views))$reference_id, ref$reference_id)
  expect_identical(build_generic(rev(p$views))$graph$accepted_edges, ref$graph$accepted_edges)
  permuted <- p$views
  permuted$D2 <- reference_view(p$views$D2$effects[, c(4, 2, 1, 3)], support = c(activity = "signed", effects = "signed"), fit_id = "D2")
  expect_identical(build_generic(permuted)$reference_id, ref$reference_id)
  set.seed(3)
  scaled <- planted_programs(noise = 0.05, scale = c(10, 0.1, 3, 7))
  rescaled <- build_generic(scaled$views)
  expect_identical(rescaled$clusters, ref$clusters)
  expect_equal(rescaled$programs, ref$programs, tolerance = 1e-10)
  path <- tempfile(fileext = ".rds")
  saveRDS(ref, path)
  again <- readRDS(path)
  expect_identical(again$reference_id, ref$reference_id)
  expect_identical(reference_programs(again), reference_programs(ref))
  tampered <- ref
  tampered$programs[1, 1] <- tampered$programs[1, 1] + 1
  expect_error(reference_programs(tampered), "content changed")
})

test_that("noiseless planted programs are recovered exactly and private ones stay singletons", {
  p <- planted_programs()
  ref <- build_generic(p$views)
  expect_identical(colnames(ref$programs), c("RP001", "RP002", "RP003"))
  expected <- sweep(p$truth, 2L, sqrt(colSums(p$truth^2)), "/")
  matched <- abs(crossprod(ref$programs, expected[rownames(ref$programs), ]))
  expect_equal(unname(apply(matched, 1, max)), rep(1, 3), tolerance = 1e-12)
  table <- reference_diagnostics(ref, "programs")
  expect_true(all(table$n_datasets[table$status == "program"] == 3L))
  expect_true(all(table$status[table$n_members == 1L] == "below_min_datasets"))
  expect_false(any(duplicated(paste(ref$clusters$cluster_id, ref$clusters$dataset_id))))
})

test_that("similarity follows the sign and metric contract", {
  B <- matrix(c(1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4), 8, 2, dimnames = list(paste0("g", 1:8), c("F1", "F2")))
  signed <- c(activity = "signed", effects = "signed")
  semi <- c(activity = "nonnegative", effects = "signed")
  collect <- function(views) .reference_collect_sources(views, names(views), .reference_manifest(generic_manifest()))
  views <- list(A = reference_view(B, support = signed, fit_id = "A"), B = reference_view(-B, support = signed, fit_id = "B"))
  table <- .reference_similarity_table(collect(views), "congruence")
  cross <- table[table$dataset_i != table$dataset_j, ]
  diagonal <- cross[cross$factor_i == cross$factor_j, ]
  expect_equal(diagonal$raw_similarity, c(-1, -1))
  expect_equal(diagonal$similarity, c(1, 1))
  expect_equal(diagonal$orientation, c(-1, -1))
  expect_equal(cross$raw_similarity[cross$factor_i != cross$factor_j], c(0, 0))
  expect_identical(table$reason[table$dataset_i == table$dataset_j], c("same_dataset", "same_dataset"))
  semi_views <- list(A = reference_view(B, support = semi, fit_id = "A"), B = reference_view(-B, support = semi, fit_id = "B"))
  semi_table <- .reference_similarity_table(collect(semi_views), "congruence")
  expect_equal(semi_table$similarity[semi_table$dataset_i != semi_table$dataset_j & semi_table$factor_i == semi_table$factor_j], c(-1, -1))
  shifted <- list(
    A = reference_view(B[, 1, drop = FALSE], support = signed, fit_id = "A"),
    B = reference_view(B[, 1, drop = FALSE] + 5, support = signed, fit_id = "B")
  )
  pearson <- .reference_similarity_table(collect(shifted), "pearson")
  congruence <- .reference_similarity_table(collect(shifted), "congruence")
  expect_equal(pearson$raw_similarity, 1)
  expect_lt(congruence$raw_similarity, 0.95)
  sparse <- list(
    A = reference_view(B[1:3, , drop = FALSE] + 1, support = signed, fit_id = "A"),
    B = reference_view(B[3:5, , drop = FALSE] + 1, support = signed, fit_id = "B")
  )
  low <- .reference_similarity_table(collect(sparse), "congruence", min_shared_features = 2L)
  expect_true(all(low$reason[low$dataset_i != low$dataset_j] == "insufficient_shared_features"))
})

test_that("all-pairs congruence equals factor_replicability before assignment and never calls LSAP", {
  skip_if_not_installed("psych")
  skip_if_not_installed("clue")
  set.seed(5)
  p <- planted_programs(n_datasets = 2L, noise = 0.3)
  rep <- factor_replicability(p$views$D1, p$views$D2)
  ref_table <- .reference_similarity_table(
    .reference_collect_sources(p$views, names(p$views), .reference_manifest(generic_manifest())),
    "congruence"
  )
  cross <- ref_table[ref_table$dataset_i != ref_table$dataset_j, ]
  expect_equal(cross$raw_similarity, rep$congruence_matrix[cbind(cross$factor_i, cross$factor_j)], tolerance = 1e-12)
  # Adversarial R023: D2 holds two near-copies of D1::F1. LSAP keeps one; the graph sees both.
  B1 <- p$views$D1$effects[, 1, drop = FALSE]
  B2 <- cbind(F1 = B1[, 1], F2 = B1[, 1] + 0.01 * seq_len(nrow(B1)))
  rownames(B2) <- rownames(B1)
  views <- list(D1 = reference_view(B1, fit_id = "D1"), D2 = reference_view(B2, fit_id = "D2"))
  testthat::local_mocked_bindings(solve_LSAP = function(...) stop("LSAP called"), .package = "clue")
  ref <- build_generic(views)
  edges <- reference_diagnostics(ref, "edges")
  expect_identical(nrow(edges), 2L)
  expect_true("duplicate_dataset" %in% reference_diagnostics(ref, "merges")$decision)
})

test_that("graph clustering handles the adversarial cases exactly", {
  expect_identical(unname(cluster_of("A1", "B1", 0.9)), list(c("A1", "B1")))
  triangle <- cluster_of(c("A1", "A1", "B1"), c("B1", "C1", "C1"), c(0.9, 0.8, 0.7))
  expect_identical(unname(triangle), list(c("A1", "B1", "C1")))
  # Chain with A1-C1 missing: merged-set density 2/3 merges at threshold 2/3, not at 0.7.
  chain <- list(c("A1", "B1"), c("B1", "C1"), c(0.9, 0.8))
  expect_length(cluster_of(chain[[1]], chain[[2]], chain[[3]]), 1L)
  expect_length(cluster_of(chain[[1]], chain[[2]], chain[[3]], density = 0.7), 2L)
  # starCAT cross-cluster rule on the same chain: 2 of 4 directed cells = 0.5, no merge.
  expect_length(cluster_of(chain[[1]], chain[[2]], chain[[3]], rule = "starcat_cross", density = 0.666), 2L)
  # Duplicate dataset: A1 and A2 both link to B1; only the stronger joins.
  dup <- cluster_of(c("A1", "A2"), c("B1", "B1"), c(0.9, 0.8))
  expect_identical(unname(dup), list(c("A1", "B1"), "A2"))
  # Equal-similarity tie: member order decides, independent of input order.
  tie_a <- cluster_of(c("A1", "A2"), c("B1", "B1"), c(0.8, 0.8))
  tie_b <- cluster_of(c("A2", "A1"), c("B1", "B1"), c(0.8, 0.8))
  expect_identical(tie_a, tie_b)
  expect_identical(tie_a[[1]], c("A1", "B1"))
  two <- cluster_of(c("A1", "A2", "C1"), c("B1", "B2", "D1"), c(0.9, 0.85, 0.6))
  expect_length(two, 3L)
  expect_true(all(lengths(two) == 2L))
})

test_that("starcat directed candidates keep the top rank_max, a strict threshold and one per dataset", {
  members <- graph_members(c("A1", "B1", "B2", "C1"))
  sim <- graph_similarity(c("A1", "A1", "A1", "B1", "B2"), c("B1", "B2", "C1", "C1", "C1"), c(0.9, 0.8, 0.5, 0.7, 0.6))
  graph <- .reference_build_graph(sim, members, 0.5, 2L, TRUE, TRUE)
  from_a <- graph$directed_candidates[graph$directed_candidates$from == "A1", ]
  expect_identical(from_a$status, c("candidate", "not_best_in_target_dataset", "outside_rank_max"))
  graph3 <- .reference_build_graph(sim, members, 0.5, 3L, TRUE, TRUE)
  expect_identical(
    graph3$directed_candidates$status[graph3$directed_candidates$from == "A1" & graph3$directed_candidates$to == "C1"],
    "not_above_min_similarity"
  )
  expect_true(all(graph$accepted_edges$reciprocal))
})

test_that("the builder fails closed on undeclared or incompatible inputs", {
  p <- planted_programs()
  expect_error(build_factor_reference(p$views, min_similarity = 0.8), "comparison_manifest")
  expect_error(build_factor_reference(p$views, comparison_manifest = generic_manifest()), "explicit min_similarity")
  expect_error(build_generic(p$views[1]), "at least two")
  bad <- generic_manifest()
  bad$preprocessing_compatibility <- "assumed"
  expect_error(build_factor_reference(p$views, comparison_manifest = bad, min_similarity = 0.8), "user_declared")
  expect_error(
    build_factor_reference(p$views, comparison_manifest = generic_manifest("counts"), min_similarity = 0.8),
    "feature_scale differs"
  )
  offset <- p$views
  offset$D1 <- reference_view(
    p$views$D1$effects,
    fit_id = "D1",
    offset = stats::setNames(rep(1, nrow(p$views$D1$effects)), rownames(p$views$D1$effects))
  )
  expect_error(build_generic(offset), "zero source offsets")
  starcat <- generic_manifest()
  starcat$mode <- "starcat2025"
  expect_error(build_factor_reference(p$views, comparison_manifest = starcat), "nonnegative feature-effect support")
  expect_error(
    build_generic(p$views, source_feature_sets = list()),
    "only to starcat2025"
  )
  standardized <- generic_manifest()
  standardized$feature_standardization <- "source_sd"
  expect_error(build_factor_reference(p$views, comparison_manifest = standardized, min_similarity = 0.8), "never applies")
})

test_that("unusable factors are recorded and incompatible clusters never become programs", {
  p <- planted_programs(support = c(activity = "nonnegative", effects = "nonnegative"))
  views <- p$views
  B <- views$D1$effects
  B[, 4] <- 0
  B[1, 3] <- -1
  views$D1 <- reference_view(B, support = c(activity = "nonnegative", effects = "nonnegative"), fit_id = "D1")
  ref <- build_generic(views)
  excluded <- reference_diagnostics(ref, "excluded")
  expect_setequal(excluded$reason, c("zero_effect", "negative_effect_under_nonnegative_support"))
  expect_false(any(c("D1::F3", "D1::F4") %in% ref$clusters$member_id))
  mixed <- planted_programs(n_datasets = 2L)
  mixed$views$D2 <- reference_view(
    mixed$views$D2$effects,
    support = c(activity = "nonnegative", effects = "signed"),
    fit_id = "D2"
  )
  ref_mixed <- build_generic(mixed$views)
  table <- reference_diagnostics(ref_mixed, "programs")
  expect_true(all(table$status[table$n_members == 2L] == "constraint_incompatible"))
  expect_identical(ncol(ref_mixed$programs), 0L)
})

test_that("starcat2025 mode reproduces the STARCAT pairwise correlation and consensus units", {
  set.seed(7)
  support <- c(activity = "nonnegative", effects = "nonnegative")
  p <- planted_programs(support = support, noise = 0.2)
  features <- rownames(p$views$D1$effects)
  sets <- list(D1 = features[1:20], D2 = features[11:30], D3 = features[5:35])
  sds <- lapply(1:3, function(d) stats::setNames(seq(0.5, 2, length.out = length(features)) * d, features))
  names(sds) <- names(sets)
  manifest <- list(mode = "starcat2025", feature_scale = "simulated_intensity", preprocessing_compatibility = "user_declared")
  ref <- build_factor_reference(p$views, comparison_manifest = manifest, source_feature_sets = sets, source_feature_sd = sds)
  expect_identical(ref$build_settings$rank_max, 3L)
  expect_identical(ref$build_settings$density_rule, "starcat_cross")
  pair <- ref$similarity[ref$similarity$member_i == "D1::F1" & ref$similarity$member_j == "D2::F1", ]
  keep <- intersect(features, union(sets$D1, sets$D2))
  spectrum <- function(d) {
    b <- p$views[[d]]$effects[, "F1"]
    (b / sum(b) * 1e6)[keep] / sds[[d]][keep]
  }
  expect_equal(pair$raw_similarity, stats::cor(spectrum("D1"), spectrum("D2")), tolerance = 1e-12)
  expect_identical(ref$feature_sets$reference_features, sort(intersect(features, unlist(sets)), method = "radix"))
  members <- ref$clusters$member_id[ref$clusters$cluster_id == colnames(ref$programs)[1]]
  profile <- function(member) {
    d <- sub("::.*", "", member)
    b <- p$views[[d]]$effects[, sub(".*::", "", member)]
    (b / sum(b) * 1e6)[rownames(ref$programs)] / sds[[d]][rownames(ref$programs)]
  }
  expect_equal(ref$programs[, 1], rowMeans(vapply(members, profile, numeric(nrow(ref$programs)))), tolerance = 1e-10)
  zero_sd <- sds
  zero_sd$D1[keep[1]] <- 0
  expect_error(
    build_factor_reference(p$views, comparison_manifest = manifest, source_feature_sets = sets, source_feature_sd = zero_sd),
    "finite and positive"
  )
})

test_that("project_reference recovers activities with the programs held fixed", {
  set.seed(11)
  p <- planted_programs(noise = 0.02)
  ref <- build_generic(p$views, projection_manifest = projection_manifest())
  expect_true(ref$projection_ready)
  truth <- matrix(stats::rnorm(30 * 3), 30, 3)
  Y <- tcrossprod(truth, ref$programs) + matrix(stats::rnorm(30 * nrow(ref$programs), sd = 0.01), 30)
  dimnames(Y) <- list(paste0("q", 1:30), rownames(ref$programs))
  out <- project_reference(ref, Y, "rows", "simulated_intensity", preprocessing = list(scope = "identity"), maxiter = 100)
  expect_s3_class(out, "reference_projection")
  expect_identical(out$representation$effects, ref$programs)
  expect_gt(min(diag(stats::cor(out$activity, truth))), 0.99)
  expect_gt(out$reconstruction[["fit_r2"]], 0.99)
  expect_identical(out$reference_id, ref$reference_id)
  columns <- project_reference(ref, t(Y), "columns", "simulated_intensity", preprocessing = list(scope = "identity"), maxiter = 100)
  expect_equal(columns$activity, out$activity, tolerance = 1e-6)
  partial <- project_reference(
    ref,
    Y[, -(1:5)],
    "rows",
    "simulated_intensity",
    preprocessing = list(scope = "identity"),
    feature_match = "intersection",
    maxiter = 100
  )
  expect_identical(partial$missing_reference_features, rownames(ref$programs)[1:5])
  expect_error(
    project_reference(ref, Y[, 1:10], "rows", "simulated_intensity", preprocessing = list(scope = "identity"), feature_match = "intersection"),
    "rank-deficient|Insufficient"
  )
  expect_error(project_reference(ref, Y, "rows", "counts", preprocessing = list(scope = "identity")), "projection_manifest")
  expect_error(project_reference(ref, Y, "rows", "simulated_intensity", preprocessing = list()), "explicit scope")
  expect_error(
    project_reference(ref, Y, "rows", "simulated_intensity", preprocessing = list(scope = "identity"), normalize_activity = "l1"),
    "nonnegative"
  )
  expect_error(
    project_reference(ref, Y, "rows", "simulated_intensity", preprocessing = list(scope = "identity"), ebnm_fn = ebnm::ebnm_point_exponential),
    "differs"
  )
  inspection <- build_generic(p$views)
  expect_false(inspection$projection_ready)
  expect_error(project_reference(inspection, Y, "rows", "simulated_intensity", preprocessing = list(scope = "identity")), "not projection-ready")
})

test_that("nonnegative references project with L1-normalized activity", {
  set.seed(12)
  support <- c(activity = "nonnegative", effects = "nonnegative")
  p <- planted_programs(support = support, noise = 0.02)
  ref <- build_generic(p$views, projection_manifest = projection_manifest("ebnm_point_exponential"))
  truth <- matrix(stats::rexp(20 * 3), 20, 3)
  truth[1, ] <- 0
  Y <- tcrossprod(truth, ref$programs)
  Y <- Y + matrix(abs(stats::rnorm(length(Y), sd = 1e-3)), nrow(Y))
  dimnames(Y) <- list(paste0("q", 1:20), rownames(ref$programs))
  out <- project_reference(ref, Y, "rows", "simulated_intensity", preprocessing = list(scope = "identity"), normalize_activity = "l1", maxiter = 100)
  expect_true(all(out$activity >= 0))
  ok <- !out$zero_sum_activity
  expect_equal(unname(rowSums(out$activity_normalized[ok, , drop = FALSE])), rep(1, sum(ok)), tolerance = 1e-12)
  expect_identical(dim(out$activity_normalized), dim(out$activity))
})

test_that("split and merged source factors cannot inflate support (R102)", {
  features <- sprintf("g%02d", 1:12)
  block <- function(at, values) replace(numeric(12), at, values)
  p1 <- block(1:4, c(4, 3, 2, 1))
  p2 <- block(5:8, c(4, 3, 2, 1))
  p3 <- block(9:12, c(4, 3, 2, 1))
  make <- function(...) {
    B <- cbind(...)
    dimnames(B) <- list(features, paste0("F", seq_len(ncol(B))))
    B
  }
  views <- list(
    D1 = reference_view(make(p1, p2, p3), fit_id = "D1"),
    D2 = reference_view(make(block(1:2, c(4, 3)), block(3:4, c(2, 1)), p2, p3), fit_id = "D2"),
    D3 = reference_view(make(p1 + p2, p3), fit_id = "D3")
  )
  ref <- build_generic(views)
  clusters <- ref$clusters
  expect_false(any(duplicated(paste(clusters$cluster_id, clusters$dataset_id))))
  split_halves <- clusters$cluster_id[clusters$member_id %in% c("D2::F1", "D2::F2")]
  expect_false(identical(split_halves[1], split_halves[2]))
  merged <- clusters$cluster_id[clusters$member_id == "D3::F1"]
  expect_false(merged %in% colnames(ref$programs))
  pair <- ref$similarity[ref$similarity$member_i == "D1::F1" & ref$similarity$member_j == "D3::F1", ]
  expect_equal(pair$similarity, 30 / sqrt(30 * 60))
  table <- reference_diagnostics(ref, "programs")
  expect_identical(table$n_datasets[table$status == "program"], c(3L, 2L, 2L))
})
