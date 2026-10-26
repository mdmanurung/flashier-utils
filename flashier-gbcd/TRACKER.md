# Implementation tracker

> **Reconciled 2026-10-26 against HEAD bc920ec.** Work shipped under the audited plan
> (docs/plans/2026-10-07-gbcd-starcat-audited-implementation-plan.md, GB-001..006),
> not by ticking this tracker. Each row was checked against R/, tests, NEWS, CI and
> local receipts: 142 done, 8 superseded by a differently named design, 13 partial
> (marked `[-]` with the gap), 1 open (custom greedy initializer). M0 evidence lives
> in local, untracked validation/ receipts. Scientific calibration is not claimed.

Legend: `[ ]` not started · `[-]` in progress · `[x]` complete · `[!]` blocked

## M0 — Baseline

- [x] Record baseline commit `75a1709e...`.
- [x] Record current package version.
- [x] Save export list.
- [x] Confirm CI green.
- [x] Snapshot ordinary EBMF tests.
- [x] Snapshot projection tests.
- [x] Snapshot stability/holdout tests.
- [x] Confirm rows/columns fixtures.

**Status:** [x]

## M1 — Foundation invariants

### Prior support

- [x] Remove stale `ebnm_point_gamma`.
- [x] Remove stale `ebnm_exponential`.
- [x] Add generalized binary.
- [x] Add nonnegative unimodal.
- [x] Add audited signed families.
- [x] Keep unknown/custom handling conservative.

### Refit policy

- [x] Add `fit_recipe`. (superseded: recipe_version/refit_policy)
- [x] Add `refit_requires_callback`. (superseded: refit_policy = "callback_required")
- [ ] Mark custom greedy initializer callback-required. (open 2026-10-26: control$greedy$init_fn is not marked callback-required)
- [x] Add stability guard before prior/config reconstruction.
- [-] Mark `refit_factors()` results callback-required. (partial: no stability-refusal test for refit results)
- [x] Mark `remove_factors()` results callback-required.
- [-] Mark fix/unfix results callback-required. (partial: unfix not tested)
- [x] Preserve policy for pure reorder.
- [x] Mark projection fits callback-required.

### Projection manifest repair

- [x] Separate source native indices from projection native indices.
- [x] Store source basis provenance.
- [x] Reset projection scaling to raw.
- [x] Reset projection multipliers to one.
- [x] Set projection native PVE from returned native fit.
- [x] Add `feature_offset` native metadata. (superseded: offset/offset_definition/offset_model)
- [x] Teach `standardize_factors()` external offset path.
- [x] Add uncertainty regression test on projection fit.
- [x] Add PVE regression test on projection fit.

**Status:** [-]

## M2 — Semi-NMF

- [x] Add `fit_seminonnegative_ebmf()`. (superseded: shipped; internal since dcb5c47)
- [x] Require native greedy v1.
- [x] Validate activity prior support.
- [x] Map priors for sample rows.
- [x] Map priors for sample columns.
- [x] Reject `ebnm_fn` collision.
- [x] Update fit recipe.
- [x] Add native-equivalence tests.
- [x] Add provenance/docs.

**Status:** [x]

## M3 — Intercept -> offset

### Fit path

- [x] Add `intercept`.
- [x] Preserve old no-intercept fast path.
- [x] Use native `flash_add_intercept()`.
- [x] Handle `max_factors=0` without greedy call.
- [x] Protect intercept from nullcheck.
- [x] Reject user nullcheck `kset` in intercept mode.
- [x] Record `factor_roles`.
- [x] Record `offset_factor_index`. (superseded: intercept_index)
- [x] Keep Intercept native index 1.

### Representation

- [x] Extract offset from RAW native intercept.
- [-] Validate fixed activity equals one. (partial: no test perturbs the ones)
- [x] Remove intercept from semantic factors.
- [x] Subset LDF to adaptive native indices.
- [x] Keep raw offset unchanged under LDF scaling.
- [x] Set adaptive `native_indices`.
- [x] Set adaptive native PVE.
- [-] Set offset definition. (partial: no test asserts offset_definition)
- [x] Test intercept-only representation.

### Native posterior mapping

- [x] Update `factor_uncertainty()`.
- [x] Update `factor_intervals()`.
- [x] Update semantic `factor_pve()`.
- [x] Update stability semantic rank.
- [-] Update stability semantic PVE. (partial: no stability-level PVE assertion)

### Edit operations

- [x] Reject intercept removal.
- [x] Reject intercept fix/unfix.
- [x] Require intercept first in reorder. (superseded: reorder pins the intercept first automatically)
- [x] Test refit preserves fixed all-ones side.

### Projection

- [x] Activity projection standardized native fit retains external offset.
- [x] Add nonzero-offset feature-projection branch.
- [x] Estimate new-feature native intercept.
- [x] Return new feature offset.
- [x] Correct feature-projection native indices/PVE.

**Status:** [-]

## M4 — Initialized EBMF

- [x] Add public function. (superseded: fit_initialized_ebmf shipped; internal since dcb5c47)
- [x] Require both initialization matrices.
- [x] Validate sample IDs.
- [x] Validate feature IDs.
- [x] Validate factor IDs.
- [x] Reorder initialization rows to X.
- [x] Reject nonfinite values.
- [x] Validate nonnegative-prior initial values.
- [x] Map semantic matrices to native axes.
- [x] Map priors to native axes.
- [x] Implement fixed activity.
- [x] Implement fixed effects.
- [x] Mandatory initial backfit.
- [x] Add greedy factors optionally.
- [x] Generate collision-free added IDs.
- [x] Nullcheck new factors descending.
- [x] Preserve initialized names after removals.
- [-] Store initialization digests. (partial: digest values not asserted)
- [x] Mark callback-required.
- [x] Add both-orientation tests.

**Status:** [-]

## M5 — GBCD adapter

- [x] Add GBCD Suggests entry.
- [x] Add pinned GBCD Remote.
- [x] Update full-engine CI.
- [x] Add `fit_gbcd_programs()`.
- [x] Validate no missing values.
- [x] Require explicit preprocessing.
- [x] Require non-unknown feature scale.
- [x] Use `kmax`, not `max_factors`.
- [x] Validate `kmax >= 2`.
- [x] Preserve sparse transpose.
- [x] Implement memory preflight.
- [x] Delegate to `gbcd::fit_gbcd()`.
- [x] Validate L/F dimensions and names.
- [x] Convert log2 effects with `log(2)`.
- [x] Keep Baseline as factor.
- [x] Set global activity support unknown.
- [x] Add per-factor support/role metadata.
- [x] Build signature-statistics table.
- [-] Record conditional LM+ash uncertainty definition. (partial: not asserted in a test)
- [x] Set native PVE unavailable.
- [x] Update provenance dependency reporting.
- [x] Test package without GBCD installed.
- [x] Add docs/credits.

**Status:** [-]

## M6 — Group prevalence

- [x] Add activity threshold.
- [x] Add direction.
- [x] Use strict `>`.
- [x] Add raw active count.
- [x] Add raw active proportion.
- [x] Add weighted active proportion when weights exist.
- [-] Record basis/scaling. (partial: no scaling field recorded)
- [x] Add manual-calculation tests.

**Status:** [-]

## M7 — Multistart

- [x] Add public function. (superseded: fit_multistart_ebmf shipped; internal since dcb5c47)
- [x] Require >=2 starts.
- [x] Generate/validate unique init seeds.
- [x] Normalize native prior list.
- [x] Derive sign constraints.
- [x] Build seeded default initializer preserving constraints.
- [x] Restrict softImpute to signed-support priors.
- [x] Restrict IRLBA to signed-support priors.
- [x] Enforce upstream missing/dense restrictions.
- [x] Reject caller greedy init callback.
- [-] Preserve remaining controls. (partial: no test passes other controls)
- [x] Capture candidate errors.
- [x] Record candidate status/reason.
- [x] Extract ELBO publicly.
- [x] Select maximum finite ELBO.
- [x] Mark selected fit callback-required.
- [x] Test reproducibility.
- [x] Test nonnegative initializer signs.

**Status:** [-]

## M8 — STRUCTURE ordering

- [x] Add `order_within_group`.
- [x] Keep input default.
- [x] Require group for MDS ordering.
- [x] Order after current subsampling.
- [x] Preserve first-appearance group order.
- [x] Handle singleton groups.
- [x] Handle zero-distance groups.
- [-] Tie-break by original displayed order. (partial: fixture has no MDS ties)
- [x] Store displayed order metadata.
- [x] Add deterministic plot tests.

**Status:** [-]

## Release

- [x] All original tests pass.
- [x] All new tests pass.
- [x] Linux release/devel pass.
- [x] macOS pass.
- [x] Windows pass.
- [x] Full optional-engine job passes.
- [x] pkgdown builds.
- [x] provenance.csv complete.
- [-] CREDITS/NOTICE reviewed. (partial: no review receipt)
- [x] NEWS records projection metadata fix.
- [-] roxygen regenerates NAMESPACE/Rd. (partial: no CI drift check)
- [-] No legacy/private FLASH workflow API copied. (partial: no test or receipt)

**Status:** [-]
