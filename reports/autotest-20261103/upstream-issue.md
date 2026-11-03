# Draft issue for ropensci-review-tools/autotest (not posted)

**Title:** `autotest_package()` aborts on functions that return ggplot2 4.x (S7) objects

autotest 0.2.0.3 stops with an internal error when a package function with an
integer argument returns a ggplot built by ggplot2 >= 4.0, whose objects are S7
(`typeof(p) == "object"`).

```
Error in as.vector(x, mode) :
  cannot coerce type 'object' to vector of type 'any'
Calls: <Anonymous> ... test_int_as_dbl.autotest_obj -> setdiff -> as.vector
```

Reproduction: `repro-s7-abort.R` builds a one-function package (`bar_plot(n = 3L)`
returns `ggplot2::ggplot(...) + geom_col()`), installs it and runs
`autotest::autotest_package(<path>, test = TRUE)`.

Environment: R 4.5.1, autotest 0.2.0.3, ggplot2 4.0.3, Linux.

`setdiff()` on the two return values calls `as.vector()`, which has no method
for S7 objects. Comparing with `identical()`, or skipping the comparison when
`typeof()` is `"object"`, would avoid the abort.
