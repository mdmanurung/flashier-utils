# Iterative review notes — v3

## Pass 1 — Review of v2 bundle itself

### Issue 1: intercept-only fit was invalid

v2 pseudocode always called `flash_greedy(Kmax=max_factors)` after adding the
intercept.

Upstream `flash_greedy()` requires `Kmax >= 1`.

Therefore:

```text
intercept="feature", max_factors=0
```

would fail.

**v3:** explicit intercept-only branch; skip greedy.

### Issue 2: LDF offset extraction was underspecified

An LDF decomposition rescales both sides of every native factor.

The intercept activity column therefore need not remain one after LDF
normalization.

**v3:** always derive the semantic offset from raw native posterior means before
LDF transformation.

---

## Pass 2 — Native-index audit

### Issue 3: uncertainty functions currently assume native K equals semantic K

`factor_uncertainty()` requires native matrix dimensions to equal the semantic
view dimensions.

`factor_intervals()` similarly assumes the sampler draw factor count equals the
view factor count.

Hiding a native intercept would break both.

**v3:** `manifest$native_indices` is the mapping authority.

### Issue 4: stability rank would count the intercept

Current stability stores:

```r
rank = fitted$n_factors
```

which would be one too large for an intercept fit.

**v3:** use semantic `ncol(view$activity)`.

---

## Pass 3 — Projection audit

This was the largest new finding.

### Issue 5: projection result inherits source native indices

A projection creates a new native FLASH fit, but `_projection_output()` copies
the source manifest.

Today this is often invisible because source and projection factors have the
same count/order.

With an intercept source, source native indices may be:

```text
2,3,...
```

while the new projection fit native indices are:

```text
1,2,...
```

This would break native uncertainty mapping.

**v3:** source indices are provenance only; projection native indices are
rebuilt from the projection fit.

### Issue 6: projection result inherits source PVE

Source PVE is not PVE of the returned projection fit.

**v3:** projection manifest stores PVE from the returned native fit.

### Issue 7: projection result inherits source LDF multipliers

The projection-native fit is initialized/fixed directly in the source numerical
basis. Its returned A/B matrices are therefore raw coordinates of the new
native fit.

Copying source LDF multipliers makes posterior uncertainty transformations
incorrect.

**v3:** projection representation uses raw projection-native scaling and unit
multipliers; source coordinate provenance is stored separately.

### Issue 8: projection native fit loses external offset

`project_factor_activity()` subtracts source offset outside FLASH.

Its returned native fit alone therefore does not remember the full model.

**v3:** attach the fixed feature offset to native projection metadata and teach
`standardize_factors()` to restore it.

### Issue 9: new-feature projection ignored source offsets

`project_factor_features()` currently returns zero offset.

Once intercept support exists, projecting new measurements/features needs to
estimate their feature-specific baseline.

**v3:** add a native intercept in the nonzero-offset feature-projection branch.

---

## Pass 4 — Refit-safety audit

### Issue 10: the problem predates the new fitting functions

Post-hoc native edits preserve the original fitting configuration.

For example, after:

```text
remove_factors
fix_factors
refit_factors
```

automatic `factor_stability()` can reconstruct the old ordinary fit rather than
the edited procedure.

Projection fits also do not carry an ordinary reproducible configuration.

**v3:** `refit_requires_callback` is assigned systematically to existing edited
and projected fits, not only to new initialized/multistart fits.

Pure factor reordering remains safe because alignment is permutation-invariant
and the model itself is unchanged.

---

## Pass 5 — GBCD audit

### Issue 11: v2 used misleading `max_factors`

Upstream `fit_gbcd()` states final rank may reach:

```text
2*Kmax - 1
```

**v3:** wrapper parameter is `kmax`.

### Issue 12: v2 globally labeled activity nonnegative

The GBCD source treats the Baseline component with a signed point-Laplace
covariance factor while non-baseline GEP memberships use the nonnegative prior.

**v3:** global activity support is `"unknown"` plus per-factor support metadata.

### Issue 13: GBCD z/lfsr do not describe the exact returned EB-SNMF point estimate

Upstream:

```text
F$lfc -> EB-SNMF posterior means
F$z_score and F$lfsr -> separate gene-wise LM + ashr procedure conditional on L
```

**v3:** preserve both but explicitly record their different statistical origins.

### Issue 14: memory risk should be surfaced

The pinned GBCD implementation can materialize a dense cell × cell covariance
matrix depending on its sparsity heuristic.

**v3:** add a preflight tied to the pinned implementation and explicit memory
budget.

---

## Pass 6 — Initialized-fit audit

### Issue 15: `backfit=FALSE` could return an unoptimized initialization

`flash_factors_init()` is initialization, not a completed optimized fit.

**v3:** initialized wrapper always runs an initial backfit.

### Issue 16: nonnegative priors need compatible initialized values

v2 did not require this.

**v3:** negative initialization under known nonnegative support errors.

### Issue 17: nullchecking several newly added factors can lose generated-ID mapping

Native removal changes subsequent indices.

**v3:** nullcheck newly added factors one at a time from highest native index
downward.

---

## Pass 7 — Multistart audit

### Issue 18: changing initializer seed via a custom init function can lose sign constraints

FLASH normally infers sign constraints and passes them to its default
initializer only when `init_fn=NULL`.

A custom seeded wrapper bypasses that behavior.

**v3:** derive and pass sign constraints explicitly.

### Issue 19: softImpute and IRLBA initializers do not accept sign constraints

Using them for NMF/semi-NMF starts could violate the intended constrained
initialization.

**v3:** signed-prior fits only for these initializers in v1.

### Issue 20: failures need to be first-class

A single bad start should not discard all candidates.

**v3:** candidate status/reason table and all-fail behavior are specified.

---

# Items intentionally removed or still rejected

- arbitrary EBNM closure serialization;
- copied GBCD covariance code;
- copied DRIFT source;
- automatic softImpute preprocessing;
- GBCD-as-FLASH masquerading;
- GBCD Baseline-as-offset;
- generic external-factor S3 framework;
- one-sided automatic initialized-factor least squares;
- nullchecking user anchors;
- public split-signed initialization without benchmarks;
- summed factor PVE as a multistart "best model" criterion.

---

# Residual risks to verify during implementation

1. LDF normalization with a native intercept for every supported LDF type.
2. Native sampler dimensions after fixed/intercept factor operations.
3. GBCD memory heuristic if upstream pin is changed.
4. Projection result backward compatibility: corrected manifests intentionally
   change projection basis IDs.
5. Sign-support classification should remain conservative; do not broaden it
   merely to make multistart accept more priors.
