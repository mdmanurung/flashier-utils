# autotest pass, 2026-11-03

autotest 0.2.0.3 cannot complete on this package as shipped, for reasons inside
autotest. With the seven `plot_*` exports excluded it completes. Every error it
then reports traces to autotest's own call reconstruction, not to package code,
and it found one real input-validation defect, now fixed.

Run: R 4.5.1, autotest 0.2.0.3, ggplot2 4.0.3, flashier 1.0.60. Source copied
from the working tree into a private library first on the path (the stale-library
artefact of the 2026-10-27 review is gone).

## Aborts inside autotest

1. **ggplot2 4.x return values.** Any function with an integer argument that
   returns a ggplot aborts autotest in `test_int_as_dbl` with
   `cannot coerce type 'object' to vector of type 'any'`: ggplot2 4.x objects are
   S7, `typeof(p) == "object"`. `repro-s7-abort.R` reproduces it with a
   one-function package; `upstream-issue.md` is a draft report (not posted).
   Workaround: `exclude` the seven `plot_*` exports.
2. **`get_Rd_param` on character tests.** After the 2026-11-03 example additions,
   autotest aborts in `test_single_char_as_random -> char_param_is_arbitrary ->
   get_Rd_param` with `argument is of length zero`. Bisecting by function points at
   `standardize_factors`, but the abort also reproduces on the source from before the
   example additions and the air reformat when that function is tested alone, while a
   full run on the same source completed. The trigger is state inside autotest, not
   a package change. `get_Rd_param` fails when the function under test has no Rd
   alias match, which `standardize_factors` does have.

## Errors reported with plot exports excluded

Run on dfd2190 (before the example additions): 42 diagnostics, 69 warnings,
8 errors.

| Function | Errors | Cause | Evidence |
|---|---|---|---|
| `fit_ebmf` | 5, "Invalid var.type." | autotest state from another example | `autotest_package(functions = "fit_ebmf")` alone reports 0 errors; a trace on flashier's validator saw only `0L` and `2L`. The errors return only when `validate_backprojection` is in the same run. |
| `validate_backprojection` | 1, "top_n must contain positive integers" | autotest passes a quoted call as an argument | the error content shows `quote(structure(list(...`; the example runs cleanly in R CMD check. |
| `rank_features` | 2, "could not find function" | internal page whose example is wrapped in `if (FALSE)` | autotest runs the guarded example without the package namespace. |

The same session found a real defect while checking the `fit_ebmf` errors:
`var_type` accepted `c(0, 1)` and `c(2, 1)`, which fail inside flashier, and
`NULL` without `S`, which fails late in ebnm. Fixed in 14239c2 for `fit_ebmf()`,
`fit_initialized_ebmf()` and `project_factor_activity()`, with tests.

## Documentation diagnostics

- `par_is_demonstrated`: 50 parameters no example used, concentrated in
  `standardize_factors`, `backproject_contrast`, `project_factor_activity`,
  `fit_ebmf`, `validate_backprojection`, `factor_activity` and `factor_effects`.
  Examples for all seven now set those options (dfd2190); each page's examples run
  without errors or warnings. Whether autotest now counts them is unconfirmed
  while its run aborts (item 2 above). Run alone, the new examples add autotest
  errors for `backproject_contrast` (2) and `fit_ebmf` (5) whose messages show the
  same quoted-object reconstruction; replaying each mutated call directly
  succeeds.
- `return_desc_includes_class`: `factor_activity`, `factor_effects` and
  `factor_features` now name `matrix` or `data.frame` in `@return`. The
  `fit_ebmf` and `rank_features` rows of this test are the artefacts above.
- `fn_without_example`: 18 internal (`@keywords internal`) functions plus the two
  deprecated aliases. Internal pages are not part of the exported API; listed for
  the srr report.

## Not done

- A clean `autotest_package(test = TRUE)` with no `error` rows: blocked by the two
  aborts inside autotest.
- The upstream issue is drafted, not filed.
