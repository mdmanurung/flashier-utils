# Implementation tracker

Legend: `[ ]` not started · `[-]` in progress · `[x]` complete · `[!]` blocked

## M0 — Baseline

- [ ] Record baseline commit `75a1709e...`.
- [ ] Record current package version.
- [ ] Save export list.
- [ ] Confirm CI green.
- [ ] Snapshot ordinary EBMF tests.
- [ ] Snapshot projection tests.
- [ ] Snapshot stability/holdout tests.
- [ ] Confirm rows/columns fixtures.

**Status:** [ ]

## M1 — Foundation invariants

### Prior support

- [ ] Remove stale `ebnm_point_gamma`.
- [ ] Remove stale `ebnm_exponential`.
- [ ] Add generalized binary.
- [ ] Add nonnegative unimodal.
- [ ] Add audited signed families.
- [ ] Keep unknown/custom handling conservative.

### Refit policy

- [ ] Add `fit_recipe`.
- [ ] Add `refit_requires_callback`.
- [ ] Mark custom greedy initializer callback-required.
- [ ] Add stability guard before prior/config reconstruction.
- [ ] Mark `refit_factors()` results callback-required.
- [ ] Mark `remove_factors()` results callback-required.
- [ ] Mark fix/unfix results callback-required.
- [ ] Preserve policy for pure reorder.
- [ ] Mark projection fits callback-required.

### Projection manifest repair

- [ ] Separate source native indices from projection native indices.
- [ ] Store source basis provenance.
- [ ] Reset projection scaling to raw.
- [ ] Reset projection multipliers to one.
- [ ] Set projection native PVE from returned native fit.
- [ ] Add `feature_offset` native metadata.
- [ ] Teach `standardize_factors()` external offset path.
- [ ] Add uncertainty regression test on projection fit.
- [ ] Add PVE regression test on projection fit.

**Status:** [ ]

## M2 — Semi-NMF

- [ ] Add `fit_seminonnegative_ebmf()`.
- [ ] Require native greedy v1.
- [ ] Validate activity prior support.
- [ ] Map priors for sample rows.
- [ ] Map priors for sample columns.
- [ ] Reject `ebnm_fn` collision.
- [ ] Update fit recipe.
- [ ] Add native-equivalence tests.
- [ ] Add provenance/docs.

**Status:** [ ]

## M3 — Intercept -> offset

### Fit path

- [ ] Add `intercept`.
- [ ] Preserve old no-intercept fast path.
- [ ] Use native `flash_add_intercept()`.
- [ ] Handle `max_factors=0` without greedy call.
- [ ] Protect intercept from nullcheck.
- [ ] Reject user nullcheck `kset` in intercept mode.
- [ ] Record `factor_roles`.
- [ ] Record `offset_factor_index`.
- [ ] Keep Intercept native index 1.

### Representation

- [ ] Extract offset from RAW native intercept.
- [ ] Validate fixed activity equals one.
- [ ] Remove intercept from semantic factors.
- [ ] Subset LDF to adaptive native indices.
- [ ] Keep raw offset unchanged under LDF scaling.
- [ ] Set adaptive `native_indices`.
- [ ] Set adaptive native PVE.
- [ ] Set offset definition.
- [ ] Test intercept-only representation.

### Native posterior mapping

- [ ] Update `factor_uncertainty()`.
- [ ] Update `factor_intervals()`.
- [ ] Update semantic `factor_pve()`.
- [ ] Update stability semantic rank.
- [ ] Update stability semantic PVE.

### Edit operations

- [ ] Reject intercept removal.
- [ ] Reject intercept fix/unfix.
- [ ] Require intercept first in reorder.
- [ ] Test refit preserves fixed all-ones side.

### Projection

- [ ] Activity projection standardized native fit retains external offset.
- [ ] Add nonzero-offset feature-projection branch.
- [ ] Estimate new-feature native intercept.
- [ ] Return new feature offset.
- [ ] Correct feature-projection native indices/PVE.

**Status:** [ ]

## M4 — Initialized EBMF

- [ ] Add public function.
- [ ] Require both initialization matrices.
- [ ] Validate sample IDs.
- [ ] Validate feature IDs.
- [ ] Validate factor IDs.
- [ ] Reorder initialization rows to X.
- [ ] Reject nonfinite values.
- [ ] Validate nonnegative-prior initial values.
- [ ] Map semantic matrices to native axes.
- [ ] Map priors to native axes.
- [ ] Implement fixed activity.
- [ ] Implement fixed effects.
- [ ] Mandatory initial backfit.
- [ ] Add greedy factors optionally.
- [ ] Generate collision-free added IDs.
- [ ] Nullcheck new factors descending.
- [ ] Preserve initialized names after removals.
- [ ] Store initialization digests.
- [ ] Mark callback-required.
- [ ] Add both-orientation tests.

**Status:** [ ]

## M5 — GBCD adapter

- [ ] Add GBCD Suggests entry.
- [ ] Add pinned GBCD Remote.
- [ ] Update full-engine CI.
- [ ] Add `fit_gbcd_programs()`.
- [ ] Validate no missing values.
- [ ] Require explicit preprocessing.
- [ ] Require non-unknown feature scale.
- [ ] Use `kmax`, not `max_factors`.
- [ ] Validate `kmax >= 2`.
- [ ] Preserve sparse transpose.
- [ ] Implement memory preflight.
- [ ] Delegate to `gbcd::fit_gbcd()`.
- [ ] Validate L/F dimensions and names.
- [ ] Convert log2 effects with `log(2)`.
- [ ] Keep Baseline as factor.
- [ ] Set global activity support unknown.
- [ ] Add per-factor support/role metadata.
- [ ] Build signature-statistics table.
- [ ] Record conditional LM+ash uncertainty definition.
- [ ] Set native PVE unavailable.
- [ ] Update provenance dependency reporting.
- [ ] Test package without GBCD installed.
- [ ] Add docs/credits.

**Status:** [ ]

## M6 — Group prevalence

- [ ] Add activity threshold.
- [ ] Add direction.
- [ ] Use strict `>`.
- [ ] Add raw active count.
- [ ] Add raw active proportion.
- [ ] Add weighted active proportion when weights exist.
- [ ] Record basis/scaling.
- [ ] Add manual-calculation tests.

**Status:** [ ]

## M7 — Multistart

- [ ] Add public function.
- [ ] Require >=2 starts.
- [ ] Generate/validate unique init seeds.
- [ ] Normalize native prior list.
- [ ] Derive sign constraints.
- [ ] Build seeded default initializer preserving constraints.
- [ ] Restrict softImpute to signed-support priors.
- [ ] Restrict IRLBA to signed-support priors.
- [ ] Enforce upstream missing/dense restrictions.
- [ ] Reject caller greedy init callback.
- [ ] Preserve remaining controls.
- [ ] Capture candidate errors.
- [ ] Record candidate status/reason.
- [ ] Extract ELBO publicly.
- [ ] Select maximum finite ELBO.
- [ ] Mark selected fit callback-required.
- [ ] Test reproducibility.
- [ ] Test nonnegative initializer signs.

**Status:** [ ]

## M8 — STRUCTURE ordering

- [ ] Add `order_within_group`.
- [ ] Keep input default.
- [ ] Require group for MDS ordering.
- [ ] Order after current subsampling.
- [ ] Preserve first-appearance group order.
- [ ] Handle singleton groups.
- [ ] Handle zero-distance groups.
- [ ] Tie-break by original displayed order.
- [ ] Store displayed order metadata.
- [ ] Add deterministic plot tests.

**Status:** [ ]

## Release

- [ ] All original tests pass.
- [ ] All new tests pass.
- [ ] Linux release/devel pass.
- [ ] macOS pass.
- [ ] Windows pass.
- [ ] Full optional-engine job passes.
- [ ] pkgdown builds.
- [ ] provenance.csv complete.
- [ ] CREDITS/NOTICE reviewed.
- [ ] NEWS records projection metadata fix.
- [ ] roxygen regenerates NAMESPACE/Rd.
- [ ] No legacy/private FLASH workflow API copied.
