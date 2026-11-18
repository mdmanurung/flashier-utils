# Cross-dataset consensus program references (flashier-starcat M1-M5).
# Conceptual method: STARCAT catalog construction (Kotliar, Curtis et al. 2025,
# Nat Methods 22:1964-1980; https://github.com/immunogenomics/starCAT, MIT).
# Independent R implementation; no STARCAT source is copied. The starcat2025
# mode follows the audited build_consensus_reference.py contract at 5399183
# (decisions D021-D023 in flashier-starcat/07_SOURCE_AUDIT_AND_DECISIONS.md).

.reference_algorithm_version <- "1.0.0"
.reference_manifest_keys <- c(
  "mode",
  "feature_scale",
  "preprocessing_compatibility",
  "profile_normalization",
  "feature_standardization"
)

.reference_sort <- function(x) sort(unique(x), method = "radix")

# Validate the comparison contract; returns it in a fixed key order for hashing.
.reference_manifest <- function(manifest) {
  if (missing(manifest) || !is.list(manifest) || is.null(names(manifest)) || anyDuplicated(names(manifest))) {
    stop("Supply comparison_manifest as a named list with mode, feature_scale and preprocessing_compatibility")
  }
  unknown <- setdiff(names(manifest), .reference_manifest_keys)
  if (length(unknown)) {
    stop("Unknown comparison_manifest fields: ", paste(unknown, collapse = ", "))
  }
  mode <- manifest$mode
  if (!is.character(mode) || length(mode) != 1L || !mode %in% c("generic", "starcat2025")) {
    stop("comparison_manifest$mode must be \"generic\" or \"starcat2025\"")
  }
  .measurement_scale(manifest$feature_scale)
  if (!identical(manifest$preprocessing_compatibility, "user_declared")) {
    stop(
      "comparison_manifest$preprocessing_compatibility must be \"user_declared\"; ",
      "declare that the source datasets were preprocessed comparably"
    )
  }
  if (is.null(manifest$profile_normalization)) {
    manifest$profile_normalization <- "auto"
  }
  if (!identical(manifest$profile_normalization, "auto")) {
    stop("comparison_manifest$profile_normalization supports only \"auto\" in schema v1")
  }
  standard <- if (mode == "starcat2025") "source_sd" else "none"
  if (is.null(manifest$feature_standardization)) {
    manifest$feature_standardization <- standard
  }
  if (!identical(manifest$feature_standardization, standard)) {
    stop(
      "comparison_manifest$feature_standardization must be \"",
      standard,
      "\" in ",
      mode,
      " mode; generic mode never applies STARCAT source-SD scaling"
    )
  }
  manifest[.reference_manifest_keys]
}

# M1 R010/R011: resolve sources, build the compact source and member tables.
.reference_collect_sources <- function(fits, dataset_id, manifest) {
  if (!is.list(fits) || inherits(fits, "flash") || inherits(fits, "factor_representation")) {
    stop("fits must be a list of native fits or views, one per dataset")
  }
  if (length(fits) < 2L) {
    stop("A consensus reference needs at least two source datasets")
  }
  if (is.null(dataset_id)) {
    stop("Name fits or supply dataset_id")
  }
  .validate_ids(dataset_id, "Dataset")
  if (length(dataset_id) != length(fits)) {
    stop("dataset_id must name every element of fits")
  }
  if (any(grepl("::", dataset_id, fixed = TRUE))) {
    stop("Dataset IDs cannot contain \"::\"; it separates dataset and factor in member IDs")
  }
  views <- stats::setNames(lapply(fits, .resolve_view), dataset_id)
  views <- views[.reference_sort(dataset_id)]
  tolerance <- function(b) 1e-8 * max(1, max(abs(b)))
  sources <- list()
  members <- list()
  for (d in names(views)) {
    view <- views[[d]]
    if (!identical(view$manifest$feature_scale, manifest$feature_scale)) {
      stop("Dataset ", d, " feature_scale differs from comparison_manifest$feature_scale")
    }
    if (any(view$offset != 0)) {
      stop("Schema-v1 references need zero source offsets; dataset ", d, " has a nonzero offset")
    }
    support <- .factor_support(view$manifest)
    B <- view$effects
    A <- view$activity
    sources[[d]] <- data.frame(
      dataset_id = d,
      basis_id = view$manifest$basis_id,
      fit_id = view$manifest$fit_id,
      feature_scale = view$manifest$feature_scale,
      preprocessing_hash = digest::digest(view$manifest$preprocessing, algo = "sha256"),
      activity_support = unname(view$manifest$constraints["activity"]),
      effect_support = unname(view$manifest$constraints["effects"]),
      n_factors = ncol(B),
      n_features = nrow(B)
    )
    if (!ncol(B)) {
      next
    }
    zero_effect <- colSums(abs(B)) == 0
    negative <- vapply(
      seq_len(ncol(B)),
      function(k) support$effect_support[k] == "nonnegative" && any(B[, k] < -tolerance(B[, k])),
      logical(1)
    )
    reason <- ifelse(
      zero_effect,
      "zero_effect",
      ifelse(negative, "negative_effect_under_nonnegative_support", "eligible")
    )
    members[[d]] <- data.frame(
      member_id = paste0(d, "::", colnames(B)),
      dataset_id = d,
      factor_id = colnames(B),
      basis_id = view$manifest$basis_id,
      feature_scale = view$manifest$feature_scale,
      activity_support = support$activity_support,
      effect_support = support$effect_support,
      n_features = nrow(B),
      zero_activity = colSums(abs(A)) == 0,
      zero_effect = unname(zero_effect),
      eligible = reason == "eligible",
      reason = reason
    )
  }
  members <- do.call(rbind, c(members, list(make.row.names = FALSE)))
  if (is.null(members) || !sum(members$eligible)) {
    stop("No eligible source factors remain")
  }
  members <- members[order(members$member_id, method = "radix"), , drop = FALSE]
  rownames(members) <- NULL
  if (manifest$mode == "starcat2025" && any(members$effect_support != "nonnegative")) {
    stop("starcat2025 mode needs nonnegative feature-effect support on every source factor")
  }
  source_table <- do.call(rbind, c(unname(sources), list(make.row.names = FALSE)))
  list(views = views, sources = source_table, members = members)
}

# Validate STARCAT source statistics (spec section 9).
.reference_starcat_inputs <- function(views, sets, sds) {
  ids <- names(views)
  if (!is.list(sets) || !is.list(sds) || !setequal(names(sets), ids) || !setequal(names(sds), ids)) {
    stop("starcat2025 mode needs source_feature_sets and source_feature_sd as lists named by dataset_id")
  }
  for (d in ids) {
    set <- sets[[d]]
    if (!is.character(set) || !length(set) || anyNA(set) || anyDuplicated(set)) {
      stop("source_feature_sets for dataset ", d, " must be unique feature IDs")
    }
    if (!all(set %in% rownames(views[[d]]$effects))) {
      stop("source_feature_sets for dataset ", d, " name features absent from its effects")
    }
    sd <- sds[[d]]
    if (!is.numeric(sd) || is.null(names(sd)) || anyDuplicated(names(sd))) {
      stop("source_feature_sd for dataset ", d, " must be a numeric vector named by feature")
    }
  }
  list(sets = sets[ids], sds = sds[ids])
}

.reference_sd <- function(sds, d, features) {
  value <- sds[[d]][features]
  if (anyNA(value) || any(!is.finite(value)) || any(value <= 0)) {
    stop("source_feature_sd for dataset ", d, " must be finite and positive on every retained feature")
  }
  value
}

# STARCAT spectrum transform: sum to 1e6 over `renorm` features, divide by the
# dataset's own feature SD, keep `keep` features. Non-positive sums give NA columns.
.reference_starcat_profile <- function(B, renorm, keep, sd) {
  X <- B[renorm, , drop = FALSE]
  total <- colSums(X)
  X <- sweep(X, 2L, ifelse(total > 0, total, NA_real_), "/") * 1e6
  X[keep, , drop = FALSE] / sd
}

# M1 R012: generic comparison profiles; L1 for nonnegative effects, otherwise L2.
.reference_generic_profile <- function(B, effect_support) {
  if (!ncol(B)) {
    return(B)
  }
  norm <- ifelse(effect_support == "nonnegative", colSums(B), sqrt(colSums(B^2)))
  sweep(B, 2L, ifelse(norm > 0, norm, NA_real_), "/")
}

.reference_similarity_block <- function(P, Q, metric) {
  if (metric == "pearson") {
    constant_P <- apply(P, 2L, .is_constant)
    constant_Q <- apply(Q, 2L, .is_constant)
    P <- sweep(P, 2L, colMeans(P), "-")
    Q <- sweep(Q, 2L, colMeans(Q), "-")
  }
  nP <- sqrt(colSums(P^2))
  nQ <- sqrt(colSums(Q^2))
  out <- crossprod(P, Q) / outer(nP, nQ)
  out[!is.finite(nP) | nP <= 0, ] <- NA_real_
  out[, !is.finite(nQ) | nQ <= 0] <- NA_real_
  if (metric == "pearson") {
    out[constant_P, ] <- NA_real_
    out[, constant_Q] <- NA_real_
  }
  pmin(pmax(out, -1), 1)
}

# M2 R020-R023: all-pairs similarity. One row per unordered member pair; no
# assignment step, so every pair is available to the graph.
.reference_similarity_table <- function(
  collected,
  metric,
  min_shared_features = 2L,
  starcat = NULL
) {
  views <- collected$views
  members <- collected$members
  rank <- stats::setNames(seq_len(nrow(members)), members$member_id)
  datasets <- names(views)
  blocks <- list()
  for (ia in seq_along(datasets)) {
    for (ib in ia:length(datasets)) {
      a <- datasets[ia]
      b <- datasets[ib]
      ma <- members[members$dataset_id == a, , drop = FALSE]
      mb <- members[members$dataset_id == b, , drop = FALSE]
      if (!nrow(ma) || !nrow(mb)) {
        next
      }
      shared <- .reference_sort(intersect(rownames(views[[a]]$effects), rownames(views[[b]]$effects)))
      keep <- shared
      if (!is.null(starcat)) {
        keep <- .reference_sort(intersect(shared, union(starcat$sets[[a]], starcat$sets[[b]])))
      }
      raw <- matrix(NA_real_, nrow(ma), nrow(mb))
      if (length(keep) >= min_shared_features) {
        Ba <- views[[a]]$effects[, ma$factor_id, drop = FALSE]
        Bb <- views[[b]]$effects[, mb$factor_id, drop = FALSE]
        if (is.null(starcat)) {
          P <- .reference_generic_profile(Ba[keep, , drop = FALSE], ma$effect_support)
          Q <- .reference_generic_profile(Bb[keep, , drop = FALSE], mb$effect_support)
        } else {
          P <- .reference_starcat_profile(Ba, shared, keep, .reference_sd(starcat$sds, a, keep))
          Q <- .reference_starcat_profile(Bb, shared, keep, .reference_sd(starcat$sds, b, keep))
        }
        raw <- .reference_similarity_block(P, Q, metric)
      }
      grid <- expand.grid(i = seq_len(nrow(ma)), j = seq_len(nrow(mb)))
      if (a == b) {
        grid <- grid[grid$i < grid$j, , drop = FALSE]
      }
      if (!nrow(grid)) {
        next
      }
      blocks[[length(blocks) + 1L]] <- data.frame(
        id_i = ma$member_id[grid$i],
        id_j = mb$member_id[grid$j],
        raw_similarity = raw[cbind(grid$i, grid$j)],
        shared_features = length(keep),
        n_i = ma$n_features[grid$i],
        n_j = mb$n_features[grid$j]
      )
    }
  }
  pairs <- do.call(rbind, blocks)
  swap <- rank[pairs$id_i] > rank[pairs$id_j]
  member_i <- ifelse(swap, pairs$id_j, pairs$id_i)
  member_j <- ifelse(swap, pairs$id_i, pairs$id_j)
  n_i <- ifelse(swap, pairs$n_j, pairs$n_i)
  n_j <- ifelse(swap, pairs$n_i, pairs$n_j)
  row_i <- match(member_i, members$member_id)
  row_j <- match(member_j, members$member_id)
  signed <- function(row) members$activity_support[row] == "signed" & members$effect_support[row] == "signed"
  flip_allowed <- signed(row_i) & signed(row_j)
  raw <- pairs$raw_similarity
  reason <- ifelse(
    !members$eligible[row_i] | !members$eligible[row_j],
    "member_ineligible",
    ifelse(
      members$dataset_id[row_i] == members$dataset_id[row_j],
      "same_dataset",
      ifelse(
        pairs$shared_features < min_shared_features,
        "insufficient_shared_features",
        ifelse(is.na(raw), "undefined_similarity", "eligible")
      )
    )
  )
  out <- data.frame(
    member_i = member_i,
    dataset_i = members$dataset_id[row_i],
    factor_i = members$factor_id[row_i],
    member_j = member_j,
    dataset_j = members$dataset_id[row_j],
    factor_j = members$factor_id[row_j],
    raw_similarity = raw,
    similarity = ifelse(flip_allowed, abs(raw), raw),
    orientation = ifelse(flip_allowed & !is.na(raw) & raw < 0, -1, 1),
    sign_flip_allowed = flip_allowed,
    shared_features = pairs$shared_features,
    coverage_i = pairs$shared_features / n_i,
    coverage_j = pairs$shared_features / n_j,
    eligible = reason == "eligible",
    reason = reason
  )
  out <- out[order(rank[out$member_i], rank[out$member_j]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# M3 R030/R031: directed candidates and accepted undirected edges.
# Directed rule (STARCAT contract, spec section 11): rank all other-dataset
# members by similarity (ties by member order), keep the first rank_max, keep
# those strictly above min_similarity, then optionally the best per target dataset.
.reference_build_graph <- function(similarity, members, min_similarity, rank_max, reciprocal, best_per_dataset) {
  rank <- stats::setNames(seq_len(nrow(members)), members$member_id)
  edges <- similarity[similarity$eligible, , drop = FALSE]
  directed <- data.frame(
    from = c(edges$member_i, edges$member_j),
    to = c(edges$member_j, edges$member_i),
    from_dataset = c(edges$dataset_i, edges$dataset_j),
    to_dataset = c(edges$dataset_j, edges$dataset_i),
    similarity = c(edges$similarity, edges$similarity)
  )
  directed <- directed[order(rank[directed$from], -directed$similarity, rank[directed$to]), , drop = FALSE]
  directed$rank <- stats::ave(seq_len(nrow(directed)), directed$from, FUN = seq_along)
  status <- rep("candidate", nrow(directed))
  if (!is.null(rank_max)) {
    status[directed$rank > rank_max] <- "outside_rank_max"
  }
  status[status == "candidate" & !(directed$similarity > min_similarity)] <- "not_above_min_similarity"
  if (best_per_dataset) {
    live <- which(status == "candidate")
    seen <- duplicated(paste(directed$from[live], directed$to_dataset[live], sep = "\r"))
    status[live[seen]] <- "not_best_in_target_dataset"
  }
  directed$status <- status
  rownames(directed) <- NULL
  candidate <- directed[directed$status == "candidate", , drop = FALSE]
  key <- function(x, y) paste(x, y, sep = "\r")
  forward <- key(candidate$from, candidate$to)
  backward <- key(candidate$to, candidate$from)
  candidate$mutual <- backward %in% forward
  first <- rank[candidate$from] < rank[candidate$to]
  keep <- if (reciprocal) candidate$mutual & first else first | !(backward %in% forward)
  accepted <- candidate[keep, , drop = FALSE]
  i_first <- rank[accepted$from] < rank[accepted$to]
  accepted <- data.frame(
    member_i = ifelse(i_first, accepted$from, accepted$to),
    member_j = ifelse(i_first, accepted$to, accepted$from),
    similarity = accepted$similarity,
    reciprocal = accepted$mutual
  )
  accepted <- accepted[order(rank[accepted$member_i], rank[accepted$member_j]), , drop = FALSE]
  rownames(accepted) <- NULL
  list(directed_candidates = directed, accepted_edges = accepted)
}

# M3 R032-R034: deterministic descending-edge merging.
# "merged_set": undirected accepted edges inside the merged set / choose(n, 2) >= threshold.
# "starcat_cross": directed candidate cells between the two sets / (2 |X| |Y|) > threshold,
# as in starCAT cluster_geps() (D021). Ties: similarity, then member order.
.reference_cluster_graph <- function(graph, members, density_rule, min_cluster_density) {
  ids <- members$member_id[members$eligible]
  n <- length(ids)
  dataset <- stats::setNames(members$dataset_id[members$eligible], ids)
  position <- stats::setNames(seq_len(n), ids)
  undirected <- matrix(FALSE, n, n)
  directed <- matrix(FALSE, n, n)
  edges <- graph$accepted_edges
  undirected[cbind(position[edges$member_i], position[edges$member_j])] <- TRUE
  undirected <- undirected | t(undirected)
  live <- graph$directed_candidates[graph$directed_candidates$status == "candidate", , drop = FALSE]
  directed[cbind(position[live$from], position[live$to])] <- TRUE
  cluster <- seq_len(n)
  edges <- edges[order(-edges$similarity, position[edges$member_i], position[edges$member_j]), , drop = FALSE]
  decisions <- vector("list", nrow(edges))
  for (e in seq_len(nrow(edges))) {
    u <- position[[edges$member_i[e]]]
    v <- position[[edges$member_j[e]]]
    X <- which(cluster == cluster[u])
    Y <- which(cluster == cluster[v])
    density <- NA_real_
    if (cluster[u] == cluster[v]) {
      decision <- "already_joined"
    } else if (length(intersect(dataset[X], dataset[Y]))) {
      decision <- "duplicate_dataset"
    } else {
      if (density_rule == "merged_set") {
        Z <- c(X, Y)
        pairs <- length(Z) * (length(Z) - 1) / 2
        density <- sum(undirected[Z, Z]) / 2 / pairs
        pass <- sum(undirected[Z, Z]) / 2 + 1e-9 >= min_cluster_density * pairs
      } else {
        density <- (sum(directed[X, Y]) + sum(directed[Y, X])) / (2 * length(X) * length(Y))
        pass <- density > min_cluster_density
      }
      decision <- if (pass) "merged" else "density_below_threshold"
      if (pass) cluster[Y] <- cluster[u]
    }
    decisions[[e]] <- data.frame(
      step = e,
      member_i = edges$member_i[e],
      member_j = edges$member_j[e],
      similarity = edges$similarity[e],
      size_i = length(X),
      size_j = length(Y),
      density = density,
      decision = decision
    )
  }
  decisions <- if (length(decisions)) {
    do.call(rbind, decisions)
  } else {
    data.frame(
      step = integer(),
      member_i = character(),
      member_j = character(),
      similarity = numeric(),
      size_i = integer(),
      size_j = integer(),
      density = numeric(),
      decision = character()
    )
  }
  groups <- split(ids, cluster)
  groups <- lapply(groups, function(g) sort(g, method = "radix"))
  size <- lengths(groups)
  first <- vapply(groups, `[`, character(1), 1L)
  groups <- unname(groups[order(-size, position[first])])
  cluster_ids <- sprintf("RP%03d", seq_along(groups))
  membership <- data.frame(
    cluster_id = rep(cluster_ids, lengths(groups)),
    member_id = unlist(groups),
    dataset_id = unname(dataset[unlist(groups)])
  )
  list(membership = membership, merge_decisions = decisions)
}

# M4 R040-R043: consensus programs on the final reference feature universe.
.reference_estimate_consensus <- function(collected, clusters, similarity, reference_features, starcat, min_datasets) {
  views <- collected$views
  members <- collected$members
  membership <- clusters$membership
  ids <- unique(membership$cluster_id)
  programs <- list()
  rows <- list()
  all_shared <- Reduce(intersect, lapply(views, function(v) rownames(v$effects)))
  for (cid in ids) {
    m <- members[match(membership$member_id[membership$cluster_id == cid], members$member_id), , drop = FALSE]
    n_datasets <- length(unique(m$dataset_id))
    pair <- similarity[similarity$member_i %in% m$member_id & similarity$member_j %in% m$member_id, , drop = FALSE]
    activity <- unique(m$activity_support)
    effect <- unique(m$effect_support)
    status <- if (length(activity) > 1L || length(effect) > 1L) {
      "constraint_incompatible"
    } else if (n_datasets < min_datasets) {
      "below_min_datasets"
    } else {
      "program"
    }
    profiles <- NULL
    if (status == "program") {
      profiles <- vapply(
        seq_len(nrow(m)),
        function(r) {
          B <- views[[m$dataset_id[r]]]$effects[, m$factor_id[r], drop = FALSE]
          if (is.null(starcat)) {
            .reference_generic_profile(B[reference_features, , drop = FALSE], m$effect_support[r])[, 1L]
          } else {
            sd <- .reference_sd(starcat$sds, m$dataset_id[r], reference_features)
            .reference_starcat_profile(B, all_shared, reference_features, sd)[, 1L]
          }
        },
        numeric(length(reference_features))
      )
      dimnames(profiles) <- list(reference_features, m$member_id)
      if (anyNA(profiles)) {
        status <- "zero_profile_on_reference_features"
      }
    }
    consensus <- NULL
    if (status == "program") {
      flip <- activity == "signed" && effect == "signed"
      if (flip && ncol(profiles) > 1L) {
        anchor <- profiles[, 1L]
        signs <- ifelse(colSums(profiles * anchor) < 0, -1, 1)
        profiles <- sweep(profiles, 2L, signs, "*")
      }
      consensus <- rowMeans(profiles)
      if (is.null(starcat)) {
        norm <- if (effect == "nonnegative") sum(consensus) else sqrt(sum(consensus^2))
        if (!(norm > 0)) {
          status <- "zero_consensus"
        } else {
          consensus <- consensus / norm
        }
      }
      if (status == "program" && flip) {
        top <- max(abs(consensus))
        anchor <- sort(names(consensus)[abs(consensus) == top], method = "radix")[1L]
        if (consensus[[anchor]] < 0) consensus <- -consensus
      }
      if (status == "program") programs[[cid]] <- consensus
    }
    fit_congruence <- if (!is.null(consensus) && status == "program") {
      .reference_similarity_block(profiles, matrix(consensus), "congruence")[, 1L]
    } else {
      NA_real_
    }
    rows[[cid]] <- data.frame(
      program_id = cid,
      status = status,
      n_members = nrow(m),
      n_datasets = n_datasets,
      activity_support = paste(sort(activity), collapse = "|"),
      effect_support = paste(sort(effect), collapse = "|"),
      mean_pair_similarity = if (nrow(pair)) mean(pair$similarity) else NA_real_,
      min_pair_similarity = if (nrow(pair)) min(pair$similarity) else NA_real_,
      min_member_consensus_congruence = min(fit_congruence),
      members = paste(m$member_id, collapse = ";")
    )
  }
  table <- do.call(rbind, c(unname(rows), list(make.row.names = FALSE)))
  matrix_out <- if (length(programs)) {
    do.call(cbind, programs)
  } else {
    matrix(numeric(), length(reference_features), 0L)
  }
  dimnames(matrix_out) <- list(reference_features, if (length(programs)) names(programs) else character())
  list(programs = matrix_out, programs_table = table)
}

# Spec section 13: the explicit common feature basis.
.reference_feature_universe <- function(views, reference_features, starcat) {
  shared <- .reference_sort(Reduce(intersect, lapply(views, function(v) rownames(v$effects))))
  default <- if (is.null(starcat)) shared else .reference_sort(intersect(shared, unlist(starcat$sets)))
  if (is.null(reference_features)) {
    features <- default
  } else {
    .validate_ids(reference_features, "Reference feature")
    if (!all(reference_features %in% shared)) {
      stop("reference_features must be measured in every source dataset")
    }
    features <- .reference_sort(reference_features)
  }
  if (length(features) < 2L) {
    stop("The reference feature universe has fewer than two features")
  }
  features
}
