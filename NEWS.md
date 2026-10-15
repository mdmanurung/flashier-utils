# flashier.utils 1.1.0.9000

- Add opt-in `standardize_factors(scaling = "rms")`: rescales each activity
  column to unit root-mean-square and compensates on the effects side,
  preserving the reconstruction exactly. This is a factor-coordinate change
  (units), not a change to the fitted model or its uncertainty; the per-factor
  scale is recorded in the basis manifest as `rms_scale`. Defaults
  (`raw`/`ldf`) are unchanged. Added to support downstream `cellprograms`
  canonicalization parity (SYNC-1A).

- Add native `fit_multistart_ebmf()` with distinct initializer seeds, explicit
  sign constraints, candidate failure/ELBO/digest records, deterministic ties,
  optional candidate retention and callback-required stability replay. Selection
  measures optimization sensitivity and does not establish biological certainty.
- Intercept, initialized and projected views now retain validated per-factor
  support and roles in their basis hashes. Ordinary legacy manifests remain
  unchanged; projection preserves fixed-side support and records estimated priors.
- Zero-factor selection, association, backprojection, network, display and
  stability consumers fail with a clear biological-factor error. Empty extraction
  and summaries remain available, as does default intercept-only backfitting.
- Automatic biological bootstrap replays recorded feature-intercept recipes in
  both orientations; incomplete intercept recipes require an explicit callback.
- Add `fit_initialized_ebmf()` for named sample/feature initialization and
  optional one-side anchoring. Initial IDs are retained, additional IDs avoid
  collisions, and automatic nullchecking only removes new adaptive factors.
  Input hashes are recorded; stability refits require an explicit callback.
- Add opt-in feature intercepts to `fit_ebmf()`. Their sample-side ones are
  protected, offsets remain separate from semantic factors, and adaptive
  selection/reordering uses explicit native maps. Ordinary fitting is unchanged.
- New-measurement projection estimates signed feature offsets for affine
  sources using fixed activities and an augmented ones column. Intercept-only
  sources are supported; dependent columns and excessive allocations are refused.
- Recognize audited generalized-binary and nonnegative-unimodal EBNM support;
  zero-factor native interval extraction returns a typed empty table.
- Add `fit_seminonnegative_ebmf()` with semantic sample/feature priors, native
  nonnegative activity and signed effects, both input orientations, and explicit
  support/refit handling for custom effect priors.
- Projected fits now use their own native factor indices, scaling, PVE,
  posterior SDs and sampler coordinates. Intersection projection reports
  unmeasured-feature uncertainty as unavailable, and returned native fits retain
  external offsets through standardization and serialization.
- Metadata tuples remain distinct when labels contain punctuation; the same
  identity rule protects donor bootstrap strata.
- Automatic stability fitting now requires an explicit callback for edited,
  projected or custom fits. Ordinary fits and pure reorders retain their recipes.
- Validate native factor maps and return typed empty summaries, rankings,
  specificity and marker correlations for valid zero-factor fits.
- `reorder_factors()` now returns a fit with a new identity, so a view made before
  the reorder is refused instead of silently pairing posterior SDs with the wrong
  factor. Provenance `engine` and `engine_version` are single strings (a partial
  `$` match had turned them into lists).
- Fix `fix_factors(ids = )`: partial IDs now fix exactly the named samples or
  features in both orientations (a single ID previously fixed the whole axis and
  several IDs errored). `fix_factors()` and `remove_factors()` now require
  `factors` instead of treating `NULL` as every factor.
- Tutorials state where their data and objects come from: a visible "Data and objects"
  section with the command to obtain the pancreas data, explicit reads of the cached
  files, a visible simulation in the protein article, and a README table of data
  sources. `restore-programs.R` now takes the cache as an argument.
- Every tutorial ends with a "Function credits" table built from `provenance()`,
  listing the functions it calls and their original authors.
- Internal simplification, no user-visible change: one engine check for optional
  packages, `.data` imported from ggplot2 instead of rlang, and the pancreas
  restore code shared across tutorials.
- Add tutorial coverage for every exported function: a new "Reorder, remove, fix
  and refit factors" article (`reorder_factors`, `remove_factors`, `fix_factors`,
  `refit_factors`, `factor_pve`) and examples for `fit_nonnegative_ebmf` (GTEx),
  `rank_features`, `factor_distinctiveness`, `rank_factors_by_metadata` and
  `factor_cor` (pancreas programs).
- `factor_gate_alignment()` gains `direction = c("high", "low")`; `"low"` compares
  a gate with the lowest-activity samples. The default is unchanged and the
  summary gains a `direction` column.
- Document that the `eligible_above = 0.5` default of `estimate_marker_thresholds()`
  is assay-specific (function manual and protein tutorial); no behaviour change.

# flashier.utils 1.1.0

- Add signed program networks, annotation comparisons, descriptive diagnostics,
  native ROC specificity, protein thresholds and equal-population gates.
- Add native fixed-activity and fixed-program EBMF projection with explicit
  units, preprocessing, coverage, rank and allocation checks.
- Add workflow tutorials and a published-data pancreas example adapted from
  Peter Carbonetto's single-cell Jamboree analysis. Tutorials now use real data
  (flashier's GTEx matrix; pancreas caches from `.github/scripts/prepare-pancreas.R`,
  built by the pkgdown workflow and gated off during package checks).

# flashier.utils 1.0.0

* Added signed decoupleR source/target/mor networks through pathways, with native
  ULM/MLM scoring, no default centering, explicit one-to-one feature mapping,
  network hashes and retained annotation losses. Duplicate edges and invalid
  weights are refused; weighted networks cannot use fgsea. Scores annotate
  estimated feature programs and condition on the EBMF fit.
* Implemented the 34 planned exports around public flashier and optional engines.
* Enforced named coordinates, biological units and fixed-basis backprojection.
* Added donor resampling/checkpoints, independent-fit matching and leakage-safe
  projection/masked scoring with retained failures and explicit denominators.
* Added Seurat/uwot and display-only ECDF; original inferential scores stay signed.
* Added source receipts, complete helper notices and actual runtime records.
* Explicitly narrowed unsupported NNLM/cycling, SuSiE, CorShrink dosym/bootstrap,
  frozen-noise projection and brms/uncertain-response modes. No fallback engines.
* Confirmed package identity flashier.utils and maintainer Mikhael Manurung.
  GitHub distribution with Linux release/devel, macOS and Windows checks.
* Full-engine CI verifies every suggested dependency loads; installs GSL for mashr.
