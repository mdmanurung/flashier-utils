test_that("every actual callable export has a source record", {
  rows <- provenance(include_dependencies = TRUE)
  expect_setequal(rows$`function`, getNamespaceExports("flashier.utils"))
  expect_true(all(nzchar(rows$implementation_path)))
  deps <- attr(rows, "dependencies")
  expect_equal(deps$version[deps$package == "flashier"], as.character(packageVersion("flashier")))
  expect_equal(provenance("provenance")$classification, "NEW_QOL")
  expect_error(provenance("missing"), "Unknown")
})

test_that("runtime accessor preserves result origin separately from its own call", {
  result<-factor_activity(fixture_view())
  out<-provenance(result=result)
  expect_equal(attr(out,"analysis_metadata"),attr(result,"analysis_metadata"))
  expect_equal(attr(out,"accessor_metadata")[["function"]],"provenance")
})
