test_that("every actual callable export has a source record", {
  rows <- provenance(include_dependencies = TRUE)
  expect_true(all(getNamespaceExports("flashier.utils") %in% rows$`function`))
  expect_true(all(c("fit_multistart_ebmf", "factor_specificity") %in% rows$`function`))
  expect_true(all(nzchar(rows$implementation_path)))
  deps <- attr(rows, "dependencies")
  expect_equal(deps$version[deps$package == "flashier"], as.character(packageVersion("flashier")))
  expect_equal(provenance("provenance")$classification, "NEW_QOL")
  expect_error(provenance("missing"), "Unknown")
})

test_that("runtime accessor preserves result origin separately from its own call", {
  result <- factor_activity(fixture_view())
  out <- provenance(result = result)
  expect_equal(attr(out, "analysis_metadata"), attr(result, "analysis_metadata"))
  expect_equal(attr(out, "accessor_metadata")[["function"]], "provenance")
})

test_that("package code calls no private flashier function", {
  # Tests may inspect flashier internals; R/ must use exported functions only.
  sources <- list.files(test_path("..", "..", "R"), pattern = "[.]R$", full.names = TRUE)
  skip_if(!length(sources), "package sources are not available in this check")
  private <- grep(":::", unlist(lapply(sources, readLines)), fixed = TRUE, value = TRUE)
  expect_identical(private, character())
})
