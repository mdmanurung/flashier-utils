test_that("descriptive summaries match IDs, units, weights and first-minus-second", {
  view <- fixture_view()
  md <- data.frame(sample=c("s3","s1","s2"),group=c("B","A","A"),donor=c("d2","d1","d1"),age=c(5,1,2))
  out <- summarize_factor_activity(view,md,"group",unit_col="donor")
  expect_equal(out$mean[out$group=="A"],c(1.5,1))
  expect_equal(out$n_units[out$group=="A"],c(1L,1L))
  weighted <- summarize_factor_activity(view,md,"group",weights=c(s1=1,s2=3,s3=1))
  expect_equal(weighted$mean[weighted$group=="A"],c(1.75,0.5))
  expect_true(all(is.na(weighted$q25)))
  contrast <- compare_factor_groups(view,md,"group",c("A","B"),unit_col="donor")
  expect_equal(contrast$factor_contrast$estimate,c(-2.5,2))
  expect_equal(contrast$feature_contrast$estimate,c(1, -1.5, 8, -2.5))
  expect_false("p.value" %in% names(contrast$factor_contrast))
  reverse <- compare_factor_groups(view,md,"group",c("B","A"))
  expect_equal(reverse$feature_contrast$estimate,-contrast$feature_contrast$estimate)
  rank <- rank_factors_by_metadata(view,md,"age")
  native <- lm(view$activity[,1] ~ md$age[match(rownames(view$activity),md$sample)])
  expect_equal(rank$value[1],summary(native)$r.squared)
})

test_that("metadata tuples remain distinct when display separators collide", {
  view <- fixture_view()
  metadata <- data.frame(sample = c("s1", "s2", "s3"),
    first = c("x.y", "x", "x.y"), second = c("z", "y.z", "z"))
  out <- summarize_factor_activity(view, metadata, c("first", "second"))
  expect_equal(nrow(out), 4)
  first <- out[out$factor == "F1", ]
  expect_equal(first$n[first$first == "x"], 1L)
  expect_equal(first$mean[first$first == "x"], 2)
  expect_equal(first$n[first$first == "x.y"], 2L)
  expect_equal(first$mean[first$first == "x.y"], 2.5)
})

## GB-005 / M6: activity prevalence summaries (strict >, raw counts, weighted proportion).
test_that("activity prevalence thresholds add raw and weighted active summaries", {
  view <- fixture_view()
  metadata <- data.frame(sample = c("s1", "s2", "s3"), group = c("g1", "g1", "g2"))
  out <- summarize_factor_activity(view, metadata, "group", activity_threshold = 1.5,
    activity_direction = "positive")
  f1 <- out[out$factor == "F1", ]
  f2 <- out[out$factor == "F2", ]
  expect_equal(f1$active_n[f1$group == "g1"], 1L)
  expect_equal(f1$active_proportion[f1$group == "g1"], 0.5)
  expect_equal(f1$active_n[f1$group == "g2"], 1L)
  expect_equal(f2$active_n[f2$group == "g1"], 1L)
  expect_equal(f2$active_n[f2$group == "g2"], 0L)
  edge <- summarize_factor_activity(view, metadata, "group", activity_threshold = 2)
  expect_equal(edge$active_n[edge$factor == "F1" & edge$group == "g1"], 0L)
  neg <- summarize_factor_activity(view, metadata, "group", activity_threshold = 0.5,
    activity_direction = "negative")
  expect_equal(neg$active_n[neg$factor == "F2" & neg$group == "g2"], 1L)
  expect_equal(neg$active_n[neg$factor == "F1" & neg$group == "g2"], 0L)
  absolute <- summarize_factor_activity(view, metadata, "group", activity_threshold = 1.5,
    activity_direction = "absolute")
  expect_equal(absolute$active_n[absolute$factor == "F2" & absolute$group == "g2"], 0L)
  expect_equal(absolute$active_n[absolute$factor == "F1" & absolute$group == "g1"], 1L)
  weights <- c(s1 = 1, s2 = 3, s3 = 2)
  weighted <- summarize_factor_activity(view, metadata, "group", weights = weights,
    activity_threshold = 1.5)
  pick <- weighted$factor == "F1" & weighted$group == "g1"
  expect_equal(weighted$active_n[pick], 1L)
  expect_equal(weighted$active_proportion[pick], 0.5)
  expect_equal(weighted$active_weighted_proportion[pick], 3 / 4)
  plain <- summarize_factor_activity(view, metadata, "group")
  expect_false(any(grepl("^active", names(plain))))
  expect_error(summarize_factor_activity(view, metadata, "group", activity_threshold = -1), "nonnegative")
  expect_error(summarize_factor_activity(view, metadata, "group", activity_threshold = c(1, 2)), "single nonnegative")
  expect_error(summarize_factor_activity(view, metadata, "group", activity_threshold = "1"), "single nonnegative")
  expect_error(summarize_factor_activity(view, metadata, "group", activity_threshold = 1, activity_direction = "up"), "positive")
})
test_that("cohort-wide metadata with extra samples is dropped, not refused", {
  view <- fixture_view()
  md <- data.frame(sample=c("s3","s1","s2"),group=c("B","A","A"))
  wide <- rbind(md, data.frame(sample="s99", group="A"))
  expect_equal(summarize_factor_activity(view,wide,"group")$mean, summarize_factor_activity(view,md,"group")$mean)
  expect_equal(compare_factor_groups(view,wide,"group",c("A","B"))$factor_contrast$estimate,
    compare_factor_groups(view,md,"group",c("A","B"))$factor_contrast$estimate)
})
