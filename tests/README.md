# Tests

Run everything from the repository root:

```sh
Rscript -e 'devtools::load_all(); testthat::test_dir("tests/testthat")'
```

A full run takes about 5 minutes. Most tests use the small fixtures in
`testthat/helper-fixtures.R` (`fixture_matrix()`, `fixture_fit()`): 12 samples by
8 features with two planted factors.

Tests of optional engines (harmony, limma, variancePartition, mashr, CorShrink,
fastTopics, gbcd and others) call `skip_if_not_installed()` and are skipped when
the package is missing. Under `R CMD check` the test that scans `R/` for `:::`
skips, because `R/` is not present in the check directory.

Tests set seeds explicitly and write only to `tempdir()`.
