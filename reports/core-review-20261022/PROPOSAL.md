# Core API migration — approved 2026-10-23

Baseline: `63f197843cd3115c072fbc8038bf4027f69a141f`; clean initial tree.

The user approved both changes ("approve both", 2026-10-23): 21 public functions, down from 55, and all nine proposed argument reorders. Keep other implementations internal so core dependencies and regression tests remain intact. Removing exports breaks external calls: provide a migration table and archived source snapshot. Both changes are applied in development version 2.0.0.9000; final local qualification is complete (1,475 source and installed assertions; R CMD check Status OK).

## Retain

`backproject_contrast()`, `compare_factor_groups()`, `factor_activity()`, `factor_effects()`, `factor_features()`, `factor_holdout()`, `factor_intervals()`, `factor_pve()`, `factor_replicability()`, `factor_stability()`, `factor_uncertainty()`, `fit_ebmf()`, `plot_factor_activity()`, `plot_factor_features()`, `project_factor_activity()`, `project_factor_features()`, `provenance()`, `standardize_factors()`, `summarize_factor_activity()`, `test_factors()`, `validate_backprojection()`.

## Remove from the public API

| Function | Replacement / rationale |
| --- | --- |
| `add_factors_to_seurat()` | Optional integration; use the upstream engine explicitly |
| `center_factor_activity()` | Descriptive convenience; use grouped summary output and explicit arithmetic/model |
| `correlate_factor_markers()` | Specialist diagnostics; keep implementation internal during migration |
| `embed_factors()` | Optional integration; use the upstream engine explicitly |
| `enrich_factors()` | Optional integration; use the upstream engine explicitly |
| `estimate_marker_thresholds()` | Specialist diagnostics; keep implementation internal during migration |
| `factor_activity_percentile()` | Optional integration; use the upstream engine explicitly |
| `factor_cor()` | Optional integration; use the upstream engine explicitly |
| `factor_diagnostics()` | Specialist diagnostics; keep implementation internal during migration |
| `factor_distinctiveness()` | Feature internals/display; factor_features and matrix operations cover the core |
| `factor_features_global()` | Feature internals/display; factor_features and matrix operations cover the core |
| `factor_gate_alignment()` | Specialist diagnostics; keep implementation internal during migration |
| `factor_network()` | Specialist diagnostics; keep implementation internal during migration |
| `factor_specificity()` | Specialist diagnostics; keep implementation internal during migration |
| `fit_initialized_ebmf()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors |
| `fit_multistart_ebmf()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors |
| `fit_nonnegative_ebmf()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors |
| `fit_seminonnegative_ebmf()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors |
| `fix_factors()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors |
| `plot_factor_contrasts()` | Specialist diagnostics; keep implementation internal during migration |
| `plot_factor_diagnostics()` | Specialist diagnostics; keep implementation internal during migration |
| `plot_factor_embedding()` | Optional integration; use the upstream engine explicitly |
| `plot_factor_gates()` | Specialist diagnostics; keep implementation internal during migration |
| `plot_factor_network()` | Specialist diagnostics; keep implementation internal during migration |
| `plot_factor_signature()` | Feature internals/display; factor_features and matrix operations cover the core |
| `plot_feature_rankings()` | Feature internals/display; factor_features and matrix operations cover the core |
| `plot_marker_thresholds()` | Specialist diagnostics; keep implementation internal during migration |
| `rank_factors_by_metadata()` | Descriptive convenience; use grouped summary output and explicit arithmetic/model |
| `rank_features()` | Feature internals/display; factor_features and matrix operations cover the core |
| `refit_factors()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors |
| `remove_factors()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors |
| `reorder_factors()` | Native fitting/editing; use public flashier calls or fit_ebmf with explicit priors |
| `shrink_factor_effects()` | Optional integration; use the upstream engine explicitly |
| `summarize_factors()` | Feature internals/display; factor_features and matrix operations cover the core |

## Argument conventions

Primary fit/data first; metadata/scientific declarations next; consistent factor/group/contrast/ID selections; computation/display controls last. Reference/target comparison functions keep their paired first arguments. The compatible audit retained existing positions; the 2.0 migration applies the approved table. A bare native fit enters through standardize_factors(sample_side=...) once; keep the native object for uncertainty. Bare-fit feature-unit/preprocessing declarations were repaired in the compatible audit.

## Four short vignettes

1. first-fit: direct native fit, standardization, extraction, one plot, features.
2. advanced-use: coordinate preservation, conditional uncertainty, one group contrast and backprojection.
3. compare-priors: same-data priors, held-out metrics, program alignment, signal/subspace and known-truth recovery.
4. pancreas-programs: one compact attributed real-data interpretation and transfer case study.

Merge essential scientific explanations; archive all original tutorials locally. Preserve CREDITS.md and installed notices, plus the existing credits.Rmd source as a non-vignette credits source if necessary. It must not become a fifth vignette.

## Progress

- Pre-migration source and installed suites: 1,459 assertions each; R CMD check Status OK.
- Complete export/caller/test/vignette inventory in exports.csv.
- Compatibility/fidelity checks, bounded prior comparison, repairs, four vignettes and final local qualification complete. The approved migration is applied; current-source qualification is complete: 1,475 assertions each in source/installed suites, native matrix 514/514, site 21 exports/four articles, R CMD check Status OK.

Detailed reviewable signature orders: proposed-arguments.csv. External caller migration: MIGRATION.md.
