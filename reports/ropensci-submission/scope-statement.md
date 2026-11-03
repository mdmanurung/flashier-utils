# Scope statement for the rOpenSci software-review issue

Draft, 2026-11-03. Paste into the "Scope" and "Overlap" fields of the submission
template.

## Scope

flashier.utils is statistical software in the **Unsupervised Learning** category
(standards G and UL). It does not fit a new model. It takes an empirical Bayes
matrix factorization fitted by flashier's public functions and returns a named
*view*: a samples × factors activity matrix and a features × factors effects
matrix that share factor IDs and reconstruct the fit. Results record the basis
they were computed on (scaling, sign, factor order, feature units), and steps
that combine results, such as back-projection and shrinkage, refuse input from
a different basis. On that view the package
ranks defining features, tests factors against sample metadata
(`test_factor_association()`), back-projects a group contrast into feature
units, and checks posterior intervals, resampling stability and held-out
prediction. `factor_provenance()` reports which upstream code and package
versions produced each result.

## Overlap with other software

- **flashier** fits the factorization. flashier.utils calls only flashier's
  exported functions (a development test checks `R/` for `:::`) and adds named,
  orientation-aware views, basis checks and provenance on top of the native
  object. It does not reimplement or change the fit.
- **MOFA2** (Bioconductor) fits a different multi-view factor model. The display
  functions follow MOFA2's downstream plot types and are independently written
  on ggplot2; no MOFA2 code is reused.
- **cNMF** (Python) and its R wrappers produce consensus NMF programs. The
  package compares saved fits, views or declared program matrices, including
  cNMF programs, through `compare_prior_fits()`, but runs no cNMF itself.
- **gbcd** and **fastTopics** are called through optional adapters
  (`fit_gbcd_programs()`, `test_gom_de()`) that check inputs and record
  provenance; their methods are unchanged.

## Names checked against CRAN

`pkgcheck::fn_names_on_cran()` on 2026-11-03 found two clashes among the 30
exports: `provenance` (AgriFusionR, dataset, metasurvey, provenance) and
`test_factors` (plspm). Both were renamed; the old names warn for one release.
