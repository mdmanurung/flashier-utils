test_that("basis identifies numerical coordinates and survives serialization", {
  view <- fixture_view()
  path <- tempfile(fileext=".rds"); on.exit(unlink(path)); saveRDS(view,path)
  expect_identical(readRDS(path),view)
  expect_identical(.make_representation(view$activity,view$effects,view$manifest)$manifest$basis_id,view$manifest$basis_id)
  changed <- view; changed$activity[1,1] <- 3
  expect_error(.validate_view(changed),"modified")
  expect_false(identical(.make_representation(changed$activity,changed$effects,changed$manifest)$manifest$basis_id,view$manifest$basis_id))
})
