test_that("native activity transfer preserves programs in both orientations", {
  v <- fixture_view(); Y <- tcrossprod(v$activity,v$effects)
  priors <- list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace)
  for(side in c("rows","columns")) {
    p <- if(side=="rows") priors else rev(priors)
    query <- if(side=="rows") Y[3:1,4:1] else t(Y[3:1,4:1])
    r <- project_factor_activity(v,query,side,"arcsinh_intensity",list(scope="identity"),ebnm_fn=p,maxiter=5)
    expect_identical(r$representation$effects,v$effects)
    expect_identical(r$representation$offset,v$offset)
    expect_equal(rownames(r$representation$activity),c("s3","s2","s1"))
    expected <- .native_projection(Y[3:1,,drop=FALSE],v,side,p,maxiter=5)
    expect_equal(r$representation$activity,expected)
    expect_equal(r$analysis_metadata$source_basis_id,v$manifest$basis_id)
    reference <- .prediction_metrics(Y[3:1,],tcrossprod(expected,v$effects)); names(reference)[names(reference)=="prediction_r2"] <- "fit_r2"
    expect_equal(r$reconstruction,reference)
  }
  sparse <- Matrix::Matrix(Y,sparse=TRUE)
  r <- project_factor_activity(v,sparse,"rows","arcsinh_intensity",list(scope="identity"),ebnm_fn=priors,maxiter=5)
  expect_identical(r$representation$effects,v$effects)
})
test_that("new feature projection agrees with explicit native calls", {
  v <- fixture_view(); A <- v$activity
  M <- cbind(P1=c(s1=1,s2=2,s3=5),P2=c(3,0,-2)); priors <- list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace)
  for(side in c("rows","columns")) {
    p <- if(side=="rows") priors else rev(priors)
    Y <- M[3:1,,drop=FALSE]; a <- A[3:1,]
    r <- project_factor_features(v,if(side=="rows") Y else t(Y),side,"protein",ebnm_fn=p,maxiter=5)
    expect_identical(r$representation$activity,a)
    b <- t(qr.coef(qr(a),Y))
    native <- flashier::flash_init(if(side=="rows") Y else t(Y),var_type=if(side=="rows") 2 else 1)
    native <- flashier::flash_factors_init(native,if(side=="rows") list(a,b) else list(b,a),ebnm_fn=p)
    native <- flashier::flash_factors_fix(native,1:2,which_dim=if(side=="rows") "loadings" else "factors")
    native <- flashier::flash_backfit(native,maxiter=5,extrapolate=FALSE,verbose=0)
    expect_equal(unname(r$representation$effects),unname(if(side=="rows") native$F_pm else native$L_pm))
  }
})
test_that("transfer rejects integrity, units, rank and budget failures", {
  v <- fixture_view(); Y <- tcrossprod(v$activity,v$effects); priors <- list(ebnm::ebnm_normal,ebnm::ebnm_point_laplace)
  run <- function(q=Y,...) project_factor_activity(v,q,"rows","arcsinh_intensity",list(scope="identity"),ebnm_fn=priors,...)
  expect_error(run(Y[,1:2]),"Exact")
  expect_error(run(Y[,1,drop=FALSE],feature_match="intersection"),"Insufficient")
  r <- run(Y[,1:3],feature_match="intersection",maxiter=5)
  expect_identical(r$representation$effects,v$effects); expect_equal(r$analysis_metadata$coverage,0.75)
  expect_error(run(dense_budget_bytes=1),"allocations")
  expect_error(project_factor_activity(v,Y,"rows","logcounts",list(scope="identity"),ebnm_fn=priors),"incompatible")
  expect_error(project_factor_activity(v,Y,"rows","arcsinh_intensity",NULL,ebnm_fn=priors),"manifest")
  expect_error(project_factor_activity(v,Y,"rows","arcsinh_intensity",list(scope="identity")),"explicit ebnm_fn")
  expect_error(run(Matrix::Matrix(Y,sparse=TRUE),S=1),"Sparse")
  dup <- Y; rownames(dup)[1] <- rownames(dup)[2]; expect_error(run(dup),"unique")
  expect_error(project_factor_features(v,Y[1,,drop=FALSE],"rows","protein",ebnm_fn=priors),"Rank-deficient")
  bad <- Y; rownames(bad)[1] <- "absent"
  expect_error(project_factor_features(v,bad,"rows","protein",ebnm_fn=priors),"absent")
  mapped <- Y; colnames(mapped) <- paste0("q",1:4)
  m <- setNames(rownames(v$effects),colnames(mapped))
  expect_equal(run(mapped,gene_mapping=m,maxiter=5)$representation$effects,v$effects)
  m[2] <- m[1]; expect_error(run(mapped,gene_mapping=m),"unique")
})
test_that("inferred priors, subsets, nonzero offsets and S retain declared coordinates", {
  f <- fixture_fit(); v <- standardize_factors(f,scaling="raw")
  Y <- fixture_matrix()
  # fixture_fit deliberately leaves units unknown: declarations fail closed.
  expect_error(project_factor_activity(v,Y,"rows","unknown",list(scope="identity")),"non-unknown")
  v <- .make_representation(v$activity,v$effects,utils::modifyList(v$manifest,list(feature_scale="simulated")))
  inferred <- suppressWarnings(project_factor_activity(v,t(Y),"columns","simulated",list(scope="identity"),maxiter=3))
  expect_identical(inferred$representation$effects,v$effects)
  new <- suppressWarnings(project_factor_features(v,Y[12:1,1:3],"rows","protein",maxiter=3))
  expect_identical(new$representation$activity,v$activity[12:1,,drop=FALSE])
  offset <- setNames(rep(1,nrow(v$effects)),rownames(v$effects))
  shifted <- .make_representation(v$activity,v$effects,v$manifest,offset=offset)
  query <- sweep(Y,2,offset,"+")
  sparse <- Matrix::Matrix(query,sparse=TRUE)
  r <- suppressWarnings(project_factor_activity(shifted,sparse,"rows","simulated",list(scope="identity"),maxiter=3))
  expect_identical(r$representation$offset,offset)
  expect_identical(r$representation$effects,v$effects)
  expect_true(is.logical(r$native_iteration_diagnostics$iteration_limit_reached))
  expect_equal(r$native_iteration_diagnostics$iteration_limit,3)
  expect_error(project_factor_activity(shifted,sparse,"rows","simulated",list(scope="identity"),dense_budget_bytes=400),"allocations")
  selected <- suppressWarnings(project_factor_activity(v,Y,"rows","simulated",list(scope="identity"),factors=rev(colnames(v$effects)),maxiter=3))
  expect_equal(colnames(selected$representation$effects),colnames(v$effects))
  S <- matrix(1,nrow(Y),ncol(Y),dimnames=dimnames(Y))
  a <- suppressWarnings(project_factor_activity(v,Y[12:1,8:1],"rows","simulated",list(scope="identity"),S=S,maxiter=3))
  b <- suppressWarnings(project_factor_activity(v,Y[12:1,8:1],"rows","simulated",list(scope="identity"),S=S[12:1,8:1],maxiter=3))
  expect_equal(a$representation$activity,b$representation$activity)
  bad <- v$effects; bad[,2] <- 0
  zero <- .make_representation(v$activity,bad,v$manifest)
  expect_error(project_factor_activity(zero,Y,"rows","simulated",list(scope="identity")),"rank-deficient")
})

test_that("projected F2 uses new native coordinates, moments, PVE and draws", {
  for (orientation in c("rows", "columns")) {
    fit <- fixture_fit(orientation)
    for (scaling in c("raw", "ldf")) {
      view <- standardize_factors(fit, scaling = scaling, d_location = "activity")
      view <- .make_representation(view$activity, view$effects,
        utils::modifyList(view$manifest, list(feature_scale = "simulated")))
      Y <- fixture_matrix()
      query <- if (orientation == "rows") Y else t(Y)
      result <- suppressWarnings(project_factor_activity(view, query, orientation,
        "simulated", list(scope = "identity"), factors = "F2", maxiter = 200))
      projected <- result$representation
      expect_identical(projected$manifest$native_indices, 1L)
      expect_identical(projected$manifest$scaling, "raw")
      expect_equal(projected$manifest$activity_multiplier, 1)
      expect_equal(projected$manifest$effects_multiplier, 1)
      expect_equal(factor_pve(projected)$pve, unname(result$fit$pve), tolerance = 1e-8)
      for (side in c("activity", "effects")) {
        prefix <- if ((side == "activity") == (orientation == "rows")) "L" else "F"
        out <- factor_uncertainty(result$fit, side, representation = projected)
        expect_equal(out$posterior.sd, as.vector(result$fit[[paste0(prefix, "_psd")]]), tolerance = 1e-8)
      }
      intervals <- factor_intervals(result$fit, "activity", representation = projected,
        factors = "F2", ids = "s2", nsamp = 20, seed = 39)
      set.seed(39)
      draws <- result$fit$sampler(20)
      native_side <- if (orientation == "rows") 1L else 2L
      values <- vapply(draws, function(d) d[[native_side]][2, 1], numeric(1))
      expect_equal(intervals$mean, mean(values), tolerance = 1e-8)
      expect_equal(intervals$lower, unname(quantile(values, .025)), tolerance = 1e-8)
      restored <- standardize_factors(result$fit, scaling = "raw")
      expect_equal(restored$activity, projected$activity)
      expect_equal(restored$effects, projected$effects)
    }
  }
})

test_that("intersection moments match native IDs and absent features are unavailable", {
  view <- fixture_view()
  Y <- tcrossprod(view$activity, view$effects)[, c("g3", "g1", "g2")]
  priors <- list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)
  for (orientation in c("rows", "columns")) {
    result <- suppressWarnings(project_factor_activity(view,
      if (orientation == "rows") Y else t(Y), orientation, "arcsinh_intensity",
      list(scope = "identity"), ebnm_fn = if (orientation == "rows") priors else rev(priors),
      feature_match = "intersection", maxiter = 200))
    out <- factor_uncertainty(result$fit, "effects", representation = result$representation)
    expect_true(all(is.na(out$posterior.sd[out$feature == "g4"])))
    expect_true(all(out$status[out$feature == "g4"] == "unavailable"))
    expect_true(all(out$reason[out$feature == "g4"] == "not_in_native_fit"))
    native <- if (orientation == "rows") result$fit$F_psd else result$fit$L_psd
    measured <- out[out$feature != "g4", ]
    index <- cbind(match(measured$feature, rownames(native)), match(measured$factor, colnames(view$effects)))
    expect_equal(measured$posterior.sd, native[index], tolerance = 1e-8)
    expect_error(factor_uncertainty(result$fit, "effects", representation = result$representation,
      require_uncertainty = TRUE), "unavailable")
    expect_error(factor_intervals(result$fit, "effects", representation = result$representation,
      ids = "g4", nsamp = 20), "unavailable")
  }
})

test_that("projected native fits retain declared external offsets after serialization", {
  view <- fixture_view()
  offset <- setNames(c(1, 2, -1, .5), rownames(view$effects))
  shifted <- .make_representation(view$activity, view$effects, view$manifest, offset)
  query <- sweep(tcrossprod(view$activity, view$effects), 2, offset, "+")
  priors <- list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)
  result <- suppressWarnings(project_factor_activity(shifted, query, "rows", "arcsinh_intensity",
    list(scope = "identity"), ebnm_fn = priors, maxiter = 200))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(result$fit, path)
  restored <- standardize_factors(readRDS(path), scaling = "raw")
  expect_identical(restored$offset, offset)
  expect_equal(sweep(tcrossprod(restored$activity, restored$effects), 2, restored$offset, "+"),
    sweep(tcrossprod(result$representation$activity, result$representation$effects), 2, offset, "+"))
})
test_that("activity projection inherits reference var_type and checks matrix S", {
  X <- fixture_matrix()
  fit <- fit_ebmf(X, "rows", ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace),
    var_type = 1L, max_factors = 2, backfit = TRUE, nullcheck = FALSE, seed = 12, verbose = 0,
    feature_scale = "intensity")
  run <- function(...) suppressWarnings(project_factor_activity(fit, X, "rows", "intensity",
    list(scope = "identity"), maxiter = 20L, ...))
  inherited <- run()
  expect_identical(inherited$analysis_metadata$var_type, 1L)
  expect_equal(inherited$representation$activity, run(var_type = 1L)$representation$activity)
  view <- standardize_factors(fit, scaling = "raw")
  expect_identical(suppressWarnings(project_factor_activity(view, X, "rows", "intensity",
    list(scope = "identity"), maxiter = 5L))$analysis_metadata$var_type, 0L)
  S <- matrix(1, nrow(X), ncol(X), dimnames = dimnames(X)); S[1, 1] <- 0
  expect_error(run(S = S), "finite and positive")
})
test_that("fixed_reference projection reuses reference activity priors in both orientations", {
  X <- fixture_matrix()
  for (side in c("rows", "columns")) {
    priors <- list(ebnm::ebnm_point_normal, ebnm::ebnm_point_laplace)
    if (side == "columns") priors <- rev(priors)
    data <- if (side == "rows") X else t(X)
    fit <- fit_ebmf(data, side, ebnm_fn = priors, max_factors = 2, backfit = TRUE,
      nullcheck = FALSE, seed = 12, verbose = 0, feature_scale = "intensity")
    run <- function(pol, ...) suppressWarnings(project_factor_activity(fit, data, side, "intensity",
      list(scope = "identity"), maxiter = 50L, prior_policy = pol, ...))
    ghat <- if (side == "rows") "L_ghat" else "F_ghat"
    fixed <- run("fixed_reference")
    expect_equal(fixed$fit[[ghat]], fit[[ghat]])
    expect_false(isTRUE(all.equal(run("query")$fit[[ghat]], fit[[ghat]])))
    expect_identical(fixed$analysis_metadata$prior_policy, "fixed_reference")
    expect_match(fixed$analysis_metadata$prior_noise_policy, "fixed_to_reference")
    expect_false(any(fixed$representation$manifest$prior_identifiers == "custom_callback_required"))
    expect_identical(fixed$representation$manifest$factor_support$activity_support, rep("signed", 2))
  }
  expect_error(run("fixed_reference", ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_normal)), "ebnm family")
  expect_error(suppressWarnings(project_factor_activity(standardize_factors(fit), data, side, "intensity",
    list(scope = "identity"), prior_policy = "fixed_reference")), "native reference fit")
})
test_that("fixed priors with sample-wise noise make activity independent of the query set", {
  set.seed(1)
  X <- tcrossprod(matrix(rnorm(80), 40, 2), matrix(rnorm(40), 20, 2)) + matrix(rnorm(800, sd = .3), 40, 20)
  dimnames(X) <- list(paste0("s", 1:40), paste0("g", 1:20))
  fit <- fit_ebmf(X, "rows", ebnm_fn = list(ebnm::ebnm_point_normal, ebnm::ebnm_point_laplace),
    max_factors = 2, backfit = TRUE, nullcheck = FALSE, seed = 2, verbose = 0, feature_scale = "intensity")
  s1 <- function(rows) suppressWarnings(project_factor_activity(fit, X[rows, ], "rows", "intensity",
    list(scope = "identity"), maxiter = 2000L, var_type = 1L, prior_policy = "fixed_reference"))$representation$activity["s1", ]
  expect_equal(s1(1:20), s1(c(1, 21:40)), tolerance = 1e-5)
})
test_that("fixed_reference maps intercept fits and factor subsets to native reference priors", {
  X <- sweep(fixture_matrix(), 2, seq_len(8), "+")
  for (side in c("rows", "columns")) {
    data <- if (side == "rows") X else t(X)
    priors <- list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace)
    if (side == "columns") priors <- rev(priors)
    fit <- fit_ebmf(data, side, ebnm_fn = priors, intercept = "feature", max_factors = 2,
      backfit = TRUE, nullcheck = FALSE, verbose = 0, seed = 12, feature_scale = "intensity")
    meta <- attr(fit, "flashbridge_metadata")
    expect_identical(meta$native_indices[1], 2L)
    ghat <- if (side == "rows") "L_ghat" else "F_ghat"
    run <- function(...) suppressWarnings(project_factor_activity(fit, data, side, "intensity",
      list(scope = "identity"), maxiter = 50L, prior_policy = "fixed_reference", ...))
    all <- run()
    expect_equal(all$fit[[ghat]], fit[[ghat]][meta$native_indices])
    last <- tail(meta$factor_ids, 1)
    one <- run(factors = last)
    expect_equal(one$fit[[ghat]], fit[[ghat]][meta$native_indices[match(last, meta$factor_ids)]])
    expect_false(isTRUE(all.equal(one$fit[[ghat]][[1]], fit[[ghat]][[1]])))
    expect_equal(one$representation$manifest$source_coordinates$native_indices,
      meta$native_indices[match(last, meta$factor_ids)])
  }
})
