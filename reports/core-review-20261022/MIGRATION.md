# Migration to the 2.0 development API

This is a breaking development API reduction. The user approved both the export and argument tables on 2026-10-23 ("approve both"). Both are applied in 2.0.0.9000; external positional calls can break. Scientific defaults and numerical algorithms are unchanged.

## Retained core signatures

`proposed-arguments.csv` records pre-migration and approved order for every retained function. It does not add arguments or change defaults. Required scientific declarations precede optional controls; native extraction consistently uses side, sample_side, representation. Single-factor operations retain their required factor ID immediately after fit. Cross-model pairs retain reference/target.

Package examples use named scientific arguments, manuals are regenerated, the site promotes only the 21 public functions, and API contract tests cover the new argument bindings. Final installed-package/native qualification passes: 1,475 assertions, native matrix 514/514, exact pre/post coordinate smoke, and R CMD check Status OK. Named calls keep their meaning. Older positional callers should stay on the prior version until migrated; do not guess their intent from argument values.

## Functions becoming internal

| Function | User migration |
| --- | --- |
| `add_factors_to_seurat()` | Optional integration; use the upstream engine explicitly. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `center_factor_activity()` | Descriptive convenience; use grouped summary output and explicit arithmetic/model. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `correlate_factor_markers()` | Specialist diagnostics; keep implementation internal during migration. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `embed_factors()` | Optional integration; use the upstream engine explicitly. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `enrich_factors()` | Optional integration; use the upstream engine explicitly. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `estimate_marker_thresholds()` | Specialist diagnostics; keep implementation internal during migration. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `factor_activity_percentile()` | Optional integration; use the upstream engine explicitly. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `factor_cor()` | Optional integration; use the upstream engine explicitly. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `factor_diagnostics()` | Specialist diagnostics; keep implementation internal during migration. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `factor_distinctiveness()` | Feature internals/display; factor_features and matrix operations cover the core. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `factor_features_global()` | Feature internals/display; factor_features and matrix operations cover the core. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `factor_gate_alignment()` | Specialist diagnostics; keep implementation internal during migration. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `factor_network()` | Specialist diagnostics; keep implementation internal during migration. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `factor_specificity()` | Specialist diagnostics; keep implementation internal during migration. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `fit_initialized_ebmf()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `fit_multistart_ebmf()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `fit_nonnegative_ebmf()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `fit_seminonnegative_ebmf()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `fix_factors()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `plot_factor_contrasts()` | Specialist diagnostics; keep implementation internal during migration. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `plot_factor_diagnostics()` | Specialist diagnostics; keep implementation internal during migration. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `plot_factor_embedding()` | Optional integration; use the upstream engine explicitly. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `plot_factor_gates()` | Specialist diagnostics; keep implementation internal during migration. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `plot_factor_network()` | Specialist diagnostics; keep implementation internal during migration. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `plot_factor_signature()` | Feature internals/display; factor_features and matrix operations cover the core. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `plot_feature_rankings()` | Feature internals/display; factor_features and matrix operations cover the core. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `plot_marker_thresholds()` | Specialist diagnostics; keep implementation internal during migration. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `rank_factors_by_metadata()` | Descriptive convenience; use grouped summary output and explicit arithmetic/model. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `rank_features()` | Feature internals/display; factor_features and matrix operations cover the core. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `refit_factors()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `remove_factors()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `reorder_factors()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `shrink_factor_effects()` | Optional integration; use the upstream engine explicitly. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |
| `summarize_factors()` | Feature internals/display; factor_features and matrix operations cover the core. The old version remains available for specialist workflows; internal calls are not a supported replacement public API. |

Preserve all origin records and licenses for retained internal source. Mark its manual pages internal and remove them from promoted navigation, rather than claiming that their historical authorship disappeared. Remove no implementation needed by the core. The provenance test must accept exported functions as a subset of historical source records.

## Documentation migration

The former `credits` article is consolidated in CREDITS.md. Old task tutorials are mapped in vignette-dispositions.json. Four current articles are first-fit, advanced-use, compare-priors and pancreas-programs.
