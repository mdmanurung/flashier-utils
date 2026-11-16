#' srr_stats
#'
#' Evidence for the rOpenSci statistical software standards (general and
#' unsupervised learning) lives in `@srrstats` tags at the end of the R file
#' it concerns. This block holds the standards that are not met, do not apply,
#' or are delegated to flashier or ebnm. Each reason starts "Not met:",
#' "Not applicable:" or "Delegated:".
#'
#' @srrstatsVerbose TRUE
#'
#' @noRd
NULL

#' NA_standards
#'
#' @srrstatsNA {G1.4a} Not met: 82 internal dot-functions (e.g. R/ids.R .validate_ids(), .validate_matrix(), .match_ids()) mostly carry plain # comments or none. R/ has only 39 @noRd or @keywords internal tags, so at least 43 have no roxygen block.
#' @srrstatsNA {G1.5} Not applicable: there is no associated publication with performance claims. Benchmark scripts and results live in benchmarks/ and benchmarks/REPORT.md, outside R/.
#' @srrstatsNA {G1.6} Not applicable: no performance claim against other R implementations. The introduction of benchmarks/REPORT.md says performance belongs to the native engine. benchmarks/ecosystem.R compares with mashr only for numerical agreement.
#' @srrstatsNA {G2.7} Not met: X accepts only a base numeric matrix or double sparse Matrix (R/ids.R .validate_matrix()). No data.frame, SummarizedExperiment or Seurat input for X. Seurat objects are accepted only as a container for adding a reduction (R/seurat.R add_factors_to_seurat()).
#' @srrstatsNA {G2.9} Not met: id_policy = 'generate' in R/ids.R .validate_matrix() inserts row/col names with no message or warning. R/fit.R fit_ebmf() only records it in fit metadata. Warnings exist only in R/batch-correction.R correct_feature_batch() and R/membership.R. Factor-to-character drops are also silent.
#' @srrstatsNA {G2.11} Not met: no handling or test for metadata columns with non-standard classes (units etc.). R/ids.R .match_ids() passes them straight to complete.cases, which may error without a clear message.
#' @srrstatsNA {G2.12} Not met: list columns in metadata are neither removed, converted nor rejected with a clear error. They reach complete.cases in R/ids.R .match_ids(). No test covers this.
#' @srrstatsNA {G2.14} Not met: NA options in R/ids.R .match_ids() are only fail or omit. R/correlation.R factor_cor() also offers use=. No imputation option, and no warning on omit, so the minimum set is incomplete. R/seurat.R factor_activity_percentile() allows only fail.
#' @srrstatsNA {G2.14b} Not met: na_action = 'omit' in R/ids.R .match_ids() drops rows and lists them in an exclusions table, but issues no warning or message. R/correlation.R factor_cor() use= is the only other ignore route.
#' @srrstatsNA {G2.14c} Not met: no imputation option exists in the package. NA in X is passed to flashier as observed missingness (R/fit.R fit_ebmf()), which is model-based handling, not an imputation choice.
#' @srrstatsNA {G2.15} Not met: mean()/sd() without na.rm in R/metadata.R summarize_factor_activity() and rank_factors_by_metadata() rely on upstream guarantees (finite activity from R/representation.R .make_representation(), complete cases from R/ids.R .match_ids()), not local guards. R/correlation.R factor_cor() does use na.rm = TRUE.
#' @srrstatsNA {G3.0} Not met: exact float tests remain, e.g. the sd(...) == 0 check in R/correlation.R factor_cor() and scale == 0 in R/comparison.R .program_unit(). Some paths use tolerances (R/backprojection.R backproject_contrast(), R/holdout.R .native_projection()), but this is not applied consistently.
#' @srrstatsNA {G3.1a} Not met: the only factor_cor() example uses "pearson" and is wrapped in @examplesIf FALSE, so it never runs (R/correlation.R factor_cor()). The engine choice is documented as an argument. No vignette or README shows the alternative engines.
#' @srrstatsNA {G5.0} Not met: tests use simulated matrices (fixture_matrix() in tests/testthat/helper-fixtures.R) and no standard reference datasets. Pinned upstream helper functions are reused for comparison (tests/testthat/helper-upstream.R, used in tests/testthat/test-features.R), but they are not datasets.
#' @srrstatsNA {G5.1} Not met: there is no data/ directory. Test datasets are built in unexported helpers (fixture_matrix() in tests/testthat/helper-fixtures.R), so users cannot rerun them.
#' @srrstatsNA {G5.2} Not met: 335 expect_error/expect_warning calls cover only part of 481 stop() calls. Errors are untyped bare stop() with only 8 call.=FALSE (e.g. R/ids.R .validate_matrix() style checks and R/fit.R). Error behaviour is not tested for every function.
#' @srrstatsNA {G5.2a} Not met: 23 stop() messages appear more than once across R/, for example "No factors selected" (twice in R/operations.R) and "level must be between zero and one" (R/backprojection.R backproject_contrast() and R/intervals.R). Messages are not unique.
#' @srrstatsNA {G5.2b} Not met: no test triggers every stop(), warning() or message(). Coverage of messages is partial.
#' @srrstatsNA {G5.3} Not met: no systematic check that returns hold no NA, NaN or Inf. Only isolated checks exist, such as the factor_uncertainty() posterior.sd check in tests/testthat/test-core-review.R and the CorShrink eigenvalue check in tests/testthat/test-correlation.R. No NA/Inf contract is documented.
#' @srrstatsNA {G5.4c} Not met: no stored values from published paper outputs appear in tests. Correctness comparisons use live oracles and pinned code, not stored published results.
#' @srrstatsNA {G5.6b} Not met: the recovery test in tests/testthat/test-fu3a-rna-adt.R uses one simulation seed (seed = 7). No recovery test repeats across simulation seeds. test-multistart.R varies init_seeds without testing recovery. test-srr-edges.R fits seeds 12 and 13 on one matrix and compares factors, not recovery of a truth.
#' @srrstatsNA {G5.7} Not met: no test shows behaviour changing with data size or convergence threshold. The tests in tests/testthat are correctness, edge and regression tests. No scaling or performance tests were found.
#' @srrstatsNA {G5.8} Not met: R/ids.R .validate_matrix() guards inputs, but tests cover only some edge cases. tests/testthat/test-ids.R covers NaN and duplicate IDs. test-srr-edges.R covers zero-extent, non-numeric and wide matrices. All-NA and constant columns in fit_ebmf() are untested, so there is no complete edge-condition suite.
#' @srrstatsNA {G5.8c} Not met: no test passes an all-NA or all-identical column to fit_ebmf() in R/fit.R. Constant-column tests exist only for downstream helpers, in tests/testthat/test-programs.R (estimate_marker_thresholds) and test-correlation.R.
#' @srrstatsNA {G5.10} Not met: the package has no extended-test flag or environment variable. Grepping tests/ and .github/ for EXTENDED_TESTS found nothing. The env block in .github/workflows/R-CMD-check.yaml sets only _R_CHECK_FORCE_SUGGESTS_.
#' @srrstatsNA {G5.11} Not applicable: no extended tests exist, and no test fetches large data. Tests use in-package fixtures from tests/testthat/helper-fixtures.R. The pancreas download lives in the pkgdown script .github/scripts/prepare-pancreas.R.
#' @srrstatsNA {G5.11a} Not applicable: no test downloads data, so no skip-on-download-failure logic is needed. The only download is download.file() in .github/scripts/prepare-pancreas.R, for the website build.
#' @srrstatsNA {UL1.4} Not met: R/fit.R fit_ebmf() documents "No preprocessing" and says feature_scale is "Recorded, never applied". It does not state assumptions about centering or scale, and no test checks what happens when they are violated.
#' @srrstatsNA {UL2.0} Not met: the only input diagnostic is rejection of Inf/NaN in .validate_matrix() (R/ids.R). No pre-processing step diagnoses scale or centering violations, and no parameter offers a transformation.
#' @srrstatsNA {UL2.2} Not met: fit_ebmf() in R/fit.R treats NA as observed missingness and rejects Inf/NaN through .validate_matrix() in R/ids.R. allow_missing exists only on that internal function. No public parameter controls NA, NaN or Inf handling.
#' @srrstatsNA {UL2.3} Not met: no collinearity check on input data X. Rank checks cover only regression designs in R/association.R and protected covariates in R/batch-correction.R (the per-cluster qr() rank check).
#' @srrstatsNA {UL3.0} Not applicable: fit_ebmf() in R/fit.R is a matrix factorization, not a clustering or partitioning method. It assigns no sequential group labels to input data. Factor IDs are F1..FK.
#' @srrstatsNA {UL3.1} Not met: fit_ebmf() in R/fit.R labels factors F1..FK in flashier greedy order, and standardize_factors() in R/representation.R has no ordering step. factor_pve() in R/extract.R reports variance explained, but ordering by it is not enforced.
#' @srrstatsNA {UL3.4} Not applicable: results of fit_ebmf() in R/fit.R are continuous factor activities and effects, not discrete groups. No clustering or partition output exists to carry intra- or inter-group variances.
#' @srrstatsNA {UL4.3} Not met: the package defines print only for prior_comparison (print.prior_comparison() in R/comparison.R; NAMESPACE has no other print method). factor_representation, gbcd_program_result and fit_ebmf() metadata have no package print method. Flash objects use flashier's own print, not verified here.
#' @srrstatsNA {UL4.3a} Not met: no package print method for factor_representation or gbcd_program_result, so default list printing can dump full matrices. The one print method, print.prior_comparison() in R/comparison.R, prints only counts and evidence state.
#' @srrstatsNA {UL4.4} Not met: no S3 summary method exists. NAMESPACE registers only print and plot for prior_comparison. R/program-diagnostics.R summarize_factors() is a plain function, not a summary() method, and is marked @keywords internal.
#' @srrstatsNA {UL6.1} Not met: no plot dispatch for the main fit or view classes points users to plot_factor_* functions (R/plots.R plot_factor_activity(), R/plot-suite.R plot_factor_weights()). Only plot.prior_comparison() in R/comparison.R exists.
#' @srrstatsNA {UL7.1} Not met: no test shows that violated input assumptions (scale, centering, NA/Inf) give unreliable output or how that shows in returned objects. Some inputs such as S fail with an error (R/projection.R .projection_S()), but no NA/Inf contract is documented.
#' @srrstatsNA {UL7.2} Not applicable: the package does not return sequential cluster labels. Factors keep flashier's order and outputs are not group partitions, so UL3.0 does not apply.
#' @srrstatsNA {UL7.4} Not met: projection of new data exists (R/projection.R project_factor_activity()), but no test or benchmark compares its runtime against refitting. No timing code appears in tests/testthat. Tests only check values (tests/testthat/test-projection.R).
#' @srrstatsNA {UL7.5} Not applicable: the package has no batch-processing routine for fitting. Row-chunked projection is deliberately not offered (R/projection.R project_factor_features() documentation). R/batch-correction.R handles technical batch effects, not batch processing.
#' @srrstatsNA {UL7.5a} Not applicable: no batch-processing routine exists to compare with direct processing. R/projection.R project_factor_features() documentation refuses chunking as an approximation.
#' @noRd
NULL
