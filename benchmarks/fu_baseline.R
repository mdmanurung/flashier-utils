pkgload::load_all("/workspace/flashier-utils", quiet = TRUE)
files <- list.files("/workspace/flashier-utils/tests/testthat", pattern = "^test-.*[.]R$", full.names = TRUE)
out <- data.frame(file = character(), passed = integer(), failed = integer(), error = integer(), skipped = integer())
for (f in files) {
  r <- tryCatch(as.data.frame(testthat::test_file(f, package = "flashier.utils", reporter = "silent")),
                error = function(e) NULL)
  if (is.null(r)) {
    out <- rbind(out, data.frame(file = basename(f), passed = NA, failed = NA, error = 1L, skipped = NA))
    next
  }
  out <- rbind(out, data.frame(file = basename(f), passed = sum(r$passed), failed = sum(r$failed),
                               error = sum(r$error), skipped = sum(r$skip)))
}
write.csv(out, "/tmp/fu_tests.csv", row.names = FALSE)
cat("TOTAL passed:", sum(out$passed, na.rm = TRUE), " failed:", sum(out$failed, na.rm = TRUE),
    " error-files:", sum(out$error, na.rm = TRUE), " skipped:", sum(out$skipped, na.rm = TRUE), "\n")
cat("--- files with failures/errors ---\n")
print(out[!is.na(out$failed) & out$failed > 0 | !is.na(out$error) & out$error > 0, ])
