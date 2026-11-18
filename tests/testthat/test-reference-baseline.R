# flashier-starcat M0 (R001/R002): legacy alignment and projection outputs must
# not change while the reference layer is built. Snapshot saved at b03ba4c by
# validation/starcat-20261118/benchmark/m0-baseline-generate.R.
test_that("legacy replicability and projection outputs match the b03ba4c snapshot", {
  skip_if_not_installed("psych")
  skip_if_not_installed("clue")
  expected <- readRDS(test_path("fixtures", "m0-legacy-baseline.rds"))
  observed <- legacy_baseline_outputs()
  expect_identical(
    observed$replicability$matches[c("reference_factor", "target_factor", "status")],
    expected$replicability$matches[c("reference_factor", "target_factor", "status")]
  )
  expect_equal(observed$replicability, expected$replicability, tolerance = 1e-6)
  expect_equal(observed$activity, expected$activity, tolerance = 1e-6)
  expect_equal(observed$fixed_reference, expected$fixed_reference, tolerance = 1e-6)
  expect_equal(observed$features, expected$features, tolerance = 1e-6)
})

test_that("planted reference fixtures are deterministic and declare their support", {
  set.seed(1)
  a <- planted_programs(noise = 0.1)
  set.seed(1)
  b <- planted_programs(noise = 0.1)
  expect_identical(a, b)
  expect_named(a$views, c("D1", "D2", "D3"))
  expect_identical(unname(a$views$D1$manifest$constraints), c("signed", "signed"))
  nonneg <- planted_programs(support = c(activity = "nonnegative", effects = "nonnegative"), noise = 0.2)
  expect_true(all(vapply(nonneg$views, function(v) all(v$effects >= 0), logical(1))))
  semi <- planted_programs(support = c(activity = "nonnegative", effects = "signed"))
  expect_identical(unname(semi$views$D2$manifest$constraints), c("nonnegative", "signed"))
})
