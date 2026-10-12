# Acceptance tests

These tests are release gates.

## A — Existing behavior

### A001 Ordinary `fit_ebmf()` unchanged

With `intercept="none"` or omitted, fixed seed/input must reproduce the current
baseline/native posterior means.

### A002 Full NMF unchanged

`fit_nonnegative_ebmf()` must reproduce its current native
point-exponential recipe.

### A003 Ordinary basis identity unchanged

For a pre-existing no-offset fit, `standardize_factors()` manifest and
`basis_id` remain unchanged.

---

# B — Foundation invariants

## B001 Prior support: generalized binary

```r
expect_identical(
  .prior_support(ebnm::ebnm_generalized_binary),
  "nonnegative"
)
```

## B002 Prior support: nonnegative unimodal

Expected `"nonnegative"`.

## B003 Stale names removed

No implementation path should claim support from nonexistent audited exports
such as `ebnm_point_gamma`.

## B004 Unknown callback remains unknown

Anonymous/parameterized callbacks are not falsely identified.

## B005 Ordinary automatic stability refit works

Reference created by standard `fit_ebmf()` with reconstructable priors may use
the default stability callback.

## B006 Edited fit refuses unsafe automatic refit

Test at least:

```text
remove_factors result
fix_factors result
refit_factors result
```

with no explicit `fit_fn`.

Expected explicit callback error.

## B007 Pure reorder remains eligible

A pure permutation does not need a different fitting recipe.

## B008 Projection fit refuses unsafe automatic refit

Projection-native fit has:

```text
refit_requires_callback = TRUE
```

and `factor_stability()` without `fit_fn` refuses.

## B009 Projection manifest native indices

For K projected factors without native intercept:

```r
expect_equal(
  result$representation$manifest$native_indices,
  seq_len(K)
)
```

regardless of source native indices.

## B010 Projection PVE belongs to projection fit

```r
expect_equal(
  result$representation$manifest$native_pve,
  result$fit$pve[result$representation$manifest$native_indices]
)
```

within tolerance.

## B011 Projection multipliers

Projection representation uses multipliers consistent with its returned native
fit. For current fixed-side projections these are ones.

## B012 Projection uncertainty

`factor_uncertainty(result$fit, representation=result$representation, ...)`
agrees with direct native PSD values.

---

# C — Semi-NMF

## C001 Rows orientation

Wrapper equals direct FLASH call with:

```text
row prior = point exponential
column prior = point Laplace
```

## C002 Columns orientation

Native prior order reverses; semantic output agrees after transposition.

## C003 Activity nonnegative

Default activity posterior means >= `-tol`.

## C004 Constraints manifest

```r
c(activity="nonnegative", effects="signed")
```

## C005 Unknown activity prior refused

A custom prior whose nonnegative support cannot be verified is not certified by
the convenience wrapper.

---

# D — Intercept/offset

## D001 Native public-pipeline equivalence

Compare to direct:

```text
flash_init
flash_add_intercept
flash_greedy when K > 0
optional backfit
adaptive-only nullcheck
```

## D002 Intercept-only model

```r
fit_ebmf(..., max_factors=0, intercept="feature")
```

must succeed and have one native factor but zero semantic factors.

## D003 Native role metadata

Exactly one native factor has role `"offset"` and is index 1.

## D004 Semantic intercept removal

`"Intercept"` is absent from semantic factor IDs.

## D005 Raw offset exactness

The semantic offset equals the raw native feature-side posterior mean of factor
1.

## D006 Fixed-one validation

Native semantic activity for offset factor is all ones.

## D007 Raw reconstruction

```r
tcrossprod(A, B) + offset
```

equals native fitted posterior mean.

## D008 LDF reconstruction, D on effects

Passes with raw offset unchanged.

## D009 LDF reconstruction, D on activity

Passes with raw offset unchanged.

## D010 Native indices

Adaptive semantic factors map to native indices 2..K_native.

## D011 PVE consistency

```r
factor_pve(fit)
factor_pve(view)
```

contain the same adaptive factor IDs/PVE and exclude Intercept.

## D012 Native posterior SD mapping

Posterior SD uses mapped native indices.

## D013 Native sampler mapping

Sampler intervals use mapped native indices and tolerate extra native offset
columns.

## D014 Stability rank

Replicate rank equals semantic factor count, not native count including offset.

## D015 Stability PVE

Replicate PVE excludes offset.

## D016 Intercept removal blocked

Generic remove operation errors.

## D017 Intercept fix/unfix blocked

Generic fix API errors.

## D018 Reorder protection

Adaptive factors may move, but Intercept must remain first.

## D019 Activity projection offset

Projection subtracts the reference offset for inference and adds it back exactly
once for reconstruction.

## D020 Projection native-fit standardization

`standardize_factors(activity_projection$fit, scaling="raw")` recovers the
external fixed offset declared in native metadata.

## D021 New-feature projection with source offset

When source offset is nonzero:

- returned new-feature offset is nonzero/estimated;
- biological activities remain fixed;
- representation reconstructs the returned native fitted mean;
- biological native indices start after the native intercept.

## D022 Zero-offset feature projection regression

The existing zero-offset branch remains numerically unchanged.

---

# E — Initialized EBMF

## E001 Exact axis validation

Reject missing/extra/duplicated sample, feature, or factor IDs.

## E002 Reordering by names

Initialization rows supplied in different order are aligned exactly to X.

## E003 Nonfinite values refused

NA/NaN/Inf initialization entries error.

## E004 Prior support consistency

Negative initialized activity errors under known nonnegative activity prior.

## E005 Direct public-pipeline equivalence

Compare to direct current FLASH public functions.

## E006 Fixed activity invariant

Fixed semantic activity remains unchanged.

## E007 Fixed effects invariant

Fixed semantic effects remain unchanged.

## E008 Both sample orientations

All fixed-side mappings pass in rows/columns mode.

## E009 Initial factors protected

Initialized factors are never automatically nullchecked.

## E010 New factor removal preserves IDs

Construct a fixture where one new factor is nullchecked away; initialized names
and surviving generated names remain correctly mapped.

## E011 Callback policy

Automatic stability refit refuses without explicit callback.

---

# F — GBCD

## F001 Optional dependency

Package loads without GBCD.

Wrapper gives clear unavailable-engine error if invoked without it.

## F002 Missing data refused

NA input errors before entering GBCD.

## F003 `kmax` validation

Values <2 error.

## F004 Orientation

Rows/columns inputs produce the same semantic cell × feature Y.

## F005 Memory preflight

A mocked/predictable dense branch exceeding budget errors before allocating the
covariance matrix.

## F006 Output dimension validation

Malformed L/F axes error.

## F007 Reconstruction scale

```r
view$effects == log(2) * res$F$lfc
```

## F008 Baseline remains factor

`"Baseline"` is retained and offset stays zero.

## F009 Global support conservative

```r
view$manifest$constraints["activity"] == "unknown"
```

## F010 Component support

Baseline marked signed; GEP components marked nonnegative.

## F011 Feature statistics preserved

`log2_fc`, `z_score`, `lfsr` exactly preserve upstream output.

## F012 Statistical provenance

Metadata explicitly distinguishes:

```text
EB-SNMF effect estimate
conditional LM + ashr z/lfsr
```

## F013 Basis linkage

Every feature-statistics row carries the representation `basis_id`.

## F014 PVE unavailable

No fake native PVE is generated.

## F015 Standard interpretation works

At least:

```text
factor_activity
factor_effects
factor_features
summarize_factor_activity
plot_factor_activity
plot_factor_features
```

work on the representation.

---

# G — Group prevalence

## G001 Positive threshold

Manual raw counts/proportions agree.

## G002 Negative threshold

Manual raw counts/proportions agree.

## G003 Absolute threshold

Manual raw counts/proportions agree.

## G004 Weighted prevalence

`active_weighted_proportion` equals weighted mean of binary active indicator.

## G005 No display scaling

Threshold uses resolved representation coordinates exactly.

---

# H — Multistart

## H001 Unique initializer seeds

Generated/provided initializer seeds are unique and recorded.

## H002 Actual initializer receives seed

Instrument initializer and verify each declared seed reaches the upstream helper.

## H003 Sign constraints preserved

For a nonnegative native dimension, all default-initializer values respect the
required sign constraint.

## H004 Unknown support refused

Native-random multistart does not silently treat unknown-support priors as
signed.

## H005 softImpute constraint rule

Nonnegative-prior multistart with softImpute initializer is refused in v1.

## H006 IRLBA missing-data rule

Missing input is refused for IRLBA initializer.

## H007 Init callback conflict

Caller-supplied greedy init function errors rather than being overwritten.

## H008 Candidate failure capture

One failing start does not discard successful candidates.

## H009 All-fail behavior

All candidate failures produce a clear aggregate error/result.

## H010 ELBO accessor

Recorded ELBO matches `flash_fit_get_elbo(flash_fit(candidate))`.

## H011 ELBO selection

Maximum finite successful ELBO selected.

## H012 Reproducibility

Same input/master seed/config produces same candidate summary.

## H013 Stability policy

Selected multistart fit requires an explicit stability callback.

---

# I — STRUCTURE ordering

## I001 Default input order

Current default output order unchanged.

## I002 MDS requires group

Ungrouped MDS request errors.

## I003 Group order

Groups appear in first-display-order sequence.

## I004 Within-group only

Samples never cross group blocks.

## I005 Singleton fallback

One-sample group remains unchanged.

## I006 Zero-distance fallback

Identical within-group factor profiles do not make plotting fail.

## I007 Determinism

Repeated call returns same displayed order.

## I008 Display-only

Underlying basis, sample values and selected sample set are unchanged.

---

# Global test invariant

No new test should depend on legacy old-API calls from DRIFT or GBCD workflow
scripts when an audited current public upstream function exists.
