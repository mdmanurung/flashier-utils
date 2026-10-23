# FU-2A (continuous membership burden, opt-in) and FU-2C (declared-baseline
# group-mean centering + global cross-factor feature selection) gate tests.
# Oracles are independent hand/table calculations on the pinned fixtures.

# ---- FU-2A -----------------------------------------------------------------

test_that("membership burden matches hand calculation on the exact fixture", {
  v <- fixture_view()  # A: F1=(1,2,4), F2=(2,0,-1); m_F1=4, m_F2=2
  d <- factor_diagnostics(v, membership_burden = TRUE)
  expect_identical(colnames(d$samples), c("sample", "programs", "membership_burden"))
  # b = pmax(A,0) %*% 1/m: s1 = 1/4 + 2/2, s2 = 2/4 + 0, s3 = 4/4 + 0 (clip -1)
  expect_equal(d$samples$membership_burden, c(1.25, 0.5, 1.0), tolerance = 1e-12)
  md <- attr(d, "analysis_metadata")$membership_burden
  expect_equal(unname(md$maxima), c(4, 2))
  expect_identical(md$zero_max_factors, character(0))
  expect_identical(md$reference_population, c("s1", "s2", "s3"))
})

test_that("membership burden is off by default and input is unchanged", {
  v <- fixture_view()
  before <- v
  d <- factor_diagnostics(v)
  expect_identical(colnames(d$samples), c("sample", "programs"))
  expect_identical(v, before)
})

test_that("zero-max factors contribute 0 and are recorded", {
  v <- fixture_view_zero_activity()  # F2 activity all zero
  d <- factor_diagnostics(v, membership_burden = TRUE)
  md <- attr(d, "analysis_metadata")$membership_burden
  expect_identical(md$zero_max_factors, "F2")
  # only F1 contributes: b = pmax(A_F1,0)/max(A_F1)
  A <- factor_activity(v)
  expect_equal(d$samples$membership_burden, unname(pmax(A[, "F1"], 0) / max(A[, "F1"])), tolerance = 1e-12)
})

test_that("extreme sample rescales burdens via m_k and maxima are recorded", {
  v <- fixture_view()
  A <- factor_activity(v)
  A2 <- rbind(A, s4 = c(100, 0))  # extreme F1 sample
  v2 <- .make_representation(A2, factor_effects(v), list(fit_id = "perturbed", schema_version = "1.0.0",
        sample_side = "rows", native_indices = 1:2, scaling = "raw", orientation = "as_fit",
        constraints = c(activity = "signed", effects = "signed"), feature_scale = "arcsinh_intensity"))
  d <- factor_diagnostics(v2, membership_burden = TRUE)
  md <- attr(d, "analysis_metadata")$membership_burden
  expect_equal(unname(md$maxima), c(100, 2))
  # s1 burden shrinks from 1.25 to 1/100 + 2/2 = 1.01
  expect_equal(d$samples$membership_burden[1], 1.01, tolerance = 1e-12)
  expect_equal(d$samples$membership_burden[4], 1.0, tolerance = 1e-12)  # s4: 100/100 + 0
})

test_that("burden respects an explicit factor subset", {
  v <- fixture_view()
  d <- factor_diagnostics(v, factors = "F2", membership_burden = TRUE)
  expect_equal(d$samples$membership_burden, c(1, 0, 0), tolerance = 1e-12)  # pmax(F2,0)/2
  expect_identical(attr(d, "analysis_metadata")$membership_burden$factors, "F2")
})

# ---- FU-2C centering ---------------------------------------------------------

test_that("center_factor_activity matches an independent table calculation", {
  v <- fixture_view()
  md <- data.frame(sample = c("s1", "s2", "s3"), grp = c("a", "a", "b"))
  A <- factor_activity(v)
  # independent oracle: group means from the raw matrix
  gm <- tapply(seq_len(nrow(A)), md$grp, function(i) colMeans(A[i, , drop = FALSE]))
  for (baseline in c("equal_group", "observation_weighted")) {
    tab <- center_factor_activity(v, md, "grp", baseline = baseline)
    expect_true(all(c("baseline", "centered") %in% colnames(tab)))
    for (f in c("F1", "F2")) {
      sub <- tab[tab$factor == f, ]
      base_expected <- if (baseline == "equal_group") {
        mean(c(gm$a[[f]], gm$b[[f]]))
      } else {
        weighted.mean(c(gm$a[[f]], gm$b[[f]]), c(2, 1))
      }
      expect_equal(unique(sub$baseline), base_expected, tolerance = 1e-12)
      expect_equal(sub$centered[sub$grp == "a"], gm$a[[f]] - base_expected, tolerance = 1e-12)
      expect_equal(sub$centered[sub$grp == "b"], gm$b[[f]] - base_expected, tolerance = 1e-12)
      # raw means and counts preserved
      expect_equal(sub$mean[sub$grp == "a"], gm$a[[f]], tolerance = 1e-12)
      expect_equal(sub$n[sub$grp == "a"], 2L)
    }
    expect_identical(attr(tab, "analysis_metadata")$centering$baseline, baseline)
    expect_false(attr(tab, "analysis_metadata")$centering$inferential)
  }
})

test_that("equal_group and observation_weighted baselines differ under imbalance", {
  v <- fixture_view()
  md <- data.frame(sample = c("s1", "s2", "s3"), grp = c("a", "a", "b"))
  eq <- center_factor_activity(v, md, "grp", baseline = "equal_group")
  ow <- center_factor_activity(v, md, "grp", baseline = "observation_weighted")
  expect_false(identical(eq$baseline[eq$factor == "F1"], ow$baseline[ow$factor == "F1"]))
  # centering never changes within-factor ordering
  o_eq <- order(eq$centered[eq$factor == "F1"])
  o_raw <- order(eq$mean[eq$factor == "F1"])
  expect_identical(o_eq, o_raw)
})

# ---- FU-2C global feature selection -------------------------------------------

test_that("factor_features_global matches independent max-abs calculation", {
  v <- fixture_view()  # B: g1=(2,3), g2=(-1,-2), g3=(0,4), g4=(1,0)
  out <- factor_features_global(v, c("F1", "F2"), n = 3)
  B <- factor_effects(v)
  score <- apply(abs(B), 1, max)
  expected_order <- rownames(B)[order(-score, rownames(B))][1:3]
  expect_identical(out$feature, expected_order)  # g3 (4), g1 (3), g2 (2)
  expect_identical(out$factor, c("F2", "F2", "F2"))
  expect_equal(out$estimate, c(4, 3, -2), tolerance = 1e-12)  # signed at argmax
  expect_identical(attr(out, "analysis_metadata")$display_tie_rule, "score_then_feature_ID")
})

test_that("global selection respects factor subset and direction", {
  v <- fixture_view()
  out <- factor_features_global(v, "F1", n = 4)
  B <- factor_effects(v)
  elig <- abs(B[, "F1"]) != 0
  expected <- rownames(B)[elig][order(-abs(B[elig, "F1"]), rownames(B)[elig])]
  expect_identical(out$feature, expected)  # zero-effect g3 excluded
  expect_true(all(out$factor == "F1"))
  neg <- factor_features_global(v, c("F1", "F2"), n = 2, direction = "negative")
  expect_true(all(neg$estimate < 0))
  expect_identical(neg$feature[1], "g2")  # most negative entry overall (-2 at F2... g2 F1=-1,F2=-2)
})

test_that("global selection has stable ties and plot mode composes", {
  A <- matrix(c(1, 2, 4, 2, 0, -1), 3, dimnames = list(c("s1", "s2", "s3"), c("F1", "F2")))
  B <- matrix(c(2, 0, -2, 1, 0, 0, 0, 0.5), 4, dimnames = list(paste0("g", 1:4), c("F1", "F2")))
  v <- .make_representation(A, B, list(fit_id = "tie_fixture", schema_version = "1.0.0",
        sample_side = "rows", native_indices = 1:2, scaling = "raw", orientation = "as_fit",
        constraints = c(activity = "signed", effects = "signed"), feature_scale = "arcsinh_intensity"))
  out <- factor_features_global(v, c("F1", "F2"), n = 2)
  # g1 and g3 tie at max|.| = 2 (g1: F1=2; g3: F1=-2); feature-ID tiebreak -> g1, g3
  expect_identical(out$feature, c("g1", "g3"))
  p <- plot_factor_features(v, c("F1", "F2"), select = "global", n = 2)
  expect_s3_class(p, "ggplot")
  expect_identical(attr(p, "analysis_metadata")$selected_features, c("g1", "g3"))
})
