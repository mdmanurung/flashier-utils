test_that("basis identifies numerical coordinates and survives serialization", {
  view <- fixture_view()
  path <- tempfile(fileext=".rds"); on.exit(unlink(path)); saveRDS(view,path)
  expect_identical(readRDS(path),view)
  expect_identical(.make_representation(view$activity,view$effects,view$manifest)$manifest$basis_id,view$manifest$basis_id)
  changed <- view; changed$activity[1,1] <- 3
  expect_error(.validate_view(changed),"modified")
  expect_false(identical(.make_representation(changed$activity,changed$effects,changed$manifest)$manifest$basis_id,view$manifest$basis_id))
})

test_that("native factor maps cannot recycle, duplicate or use invalid indices", {
  view <- fixture_view()
  for (indices in list(1L, c(1L, 1L), c(0L, 2L), c(1.5, 2), c(NA, 2))) {
    manifest <- utils::modifyList(view$manifest, list(native_indices = indices))
    expect_error(.make_representation(view$activity, view$effects, manifest), "native_indices")
  }
  for (key in c("activity_multiplier", "effects_multiplier")) {
    manifest <- view$manifest
    manifest[[key]] <- 1
    expect_error(.make_representation(view$activity, view$effects, manifest), "multiplier")
  }
})
