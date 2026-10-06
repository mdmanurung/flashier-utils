test_that("core namespace is available without optional engines", {
  expect_true(requireNamespace("flashier.utils", quietly = TRUE))
  expect_false(any(c("brms", "Seurat", "variancePartition") %in% names(getNamespaceImports("flashier.utils"))))
})
