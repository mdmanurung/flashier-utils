# Patch map

This file translates the specification into concrete repository edits.

## `R/fit.R`

### Modify `.prior_support()`

- remove stale nonexistent exports;
- add current audited signed/nonnegative families;
- remain conservative for unknown callbacks/nonpositive families.

### Modify `fit_ebmf()`

Add:

```r
intercept = c("none", "feature")
```

Required changes:

- old direct `flash()` branch remains for `intercept == "none"`;
- intercept branch uses `flash_init()` + `flash_add_intercept()`;
- skip greedy when `max_factors == 0L`;
- exclude native factor 1 from nullcheck;
- protect `control$nullcheck$kset` in intercept mode;
- populate `factor_roles`, `offset_factor_index`, `fit_recipe`,
  `refit_requires_callback`;
- include `intercept` in stored `configuration`;
- mark callback-required if custom `control$greedy$init_fn` is used.

### Add helper(s)

Recommended internal helpers:

```r
.native_prior_list()
.fit_recipe_default()
.set_refit_policy()
```

Keep them small.

---

## `R/operations.R`

### Add `fit_seminonnegative_ebmf()`

Thin wrapper around `fit_ebmf()`.

### Modify `.native_changed()`

Preserve all offset/recipe metadata.

Set callback requirement by operation:

```text
reorder -> preserve existing policy
backfit -> TRUE
remove  -> TRUE
fix     -> TRUE
unfix   -> TRUE
```

### Modify `remove_factors()`

If intercept fit:

```text
reject "Intercept"
```

### Modify `fix_factors()`

If intercept fit:

```text
reject "Intercept"
```

### Modify `reorder_factors()`

If intercept fit:

```text
require first requested ID == "Intercept"
```

Adaptive factors may still be permuted.

### Add `fit_initialized_ebmf()`

Could live here or in a new `R/initialized.R`.

Required helpers:

```r
.validate_initialized_basis()
.semantic_to_native_pair()
.semantic_to_native_priors()
.next_factor_ids()
.nullcheck_new_factors_descending()
```

---

## `R/representation.R`

### Modify `standardize_factors()`

Support two offset sources:

1. native intercept:
   ```r
   metadata$offset_factor_index
   ```
2. fixed external offset:
   ```r
   metadata$feature_offset
   ```

Precedence:

```text
native intercept > external feature offset > zero
```

Reject simultaneous conflicting definitions.

For native intercept:

- derive offset from raw semantic posterior means;
- validate fixed activity equals one;
- remove the native intercept from factor axes;
- retain actual adaptive native indices.

For external offset:

- validate finite named vector exactly matches semantic features;
- retain all native factors.

Only add offset-specific manifest keys when such an offset exists.

### Review `.make_representation()`

No schema rewrite is needed.

It already hashes the offset into `basis_id`.

Ensure:

```r
length(manifest$native_indices) == ncol(activity)
```

is explicitly validated if not already.

---

## `R/extract.R`

### Modify `factor_uncertainty()`

Replace native/semantic dimension identity assumption.

Algorithm:

```r
native_indices <- view$manifest$native_indices
native_mat <- native_mat[, native_indices, drop = FALSE]
```

then apply semantic ID selection and multipliers.

Validate native row axis independently.

### Modify `factor_pve()`

For wrapped native fits with `offset_factor_index`:

- exclude offset factor;
- return adaptive semantic IDs/PVE only.

For a representation:

- preserve existing behavior;
- error clearly when `native_pve` unavailable.

---

## `R/intervals.R`

### Modify `factor_intervals()`

Sampler draws may contain more native factors than the semantic view.

Use:

```r
native_indices <- view$manifest$native_indices
```

for sampler column lookup.

Replace exact full-dimension equality check with:

```text
correct native row count
sufficient native factor count
```

then select native columns.

---

## `R/stability.R`

### Automatic refit guard

Before prior/config reconstruction:

```r
if (isTRUE(native_metadata$refit_requires_callback))
    require fit_fn
```

### Semantic replicate rank/PVE

After refit:

```r
view <- standardize_factors(fitted, ...)
rank <- ncol(view$activity)
pve <- factor_pve(view)
```

### Projection/edited fits

Tests must prove automatic refit refuses when those fits are used as references.

---

## `R/projection.R`

This file needs a targeted correctness refactor before intercept support.

### Modify `_projection_output()`

Do not inherit projection-native metadata blindly from the source view.

Store source fields separately:

```text
source_basis_id
source_factor_ids
source_native_indices
source_scaling
source_orientation
```

Build projection-fit native fields from the returned native fit:

```text
native_indices
native_pve
scaling = raw
multipliers = 1
orientation = as_fit
```

Attach:

```r
refit_requires_callback = TRUE
fit_recipe$method = "projection"
```

### `project_factor_activity()`

Keep current mathematical offset subtraction.

On returned native fit attach:

```r
feature_offset = source$offset
offset_definition = "fixed_external_projection_offset"
```

so `standardize_factors(result$fit)` can reconstruct the full semantic fit.

### `project_factor_features()`

Keep current zero-offset branch unchanged.

Add nonzero-offset branch:

```text
flash_init
-> flash_add_intercept
-> residual extraction
-> initialize biological effects
-> fix biological activity
-> joint backfit
-> offset = native intercept effect
-> biological native_indices = 2:(K+1)
```

Returned representation and native fit must agree through the new offset-aware
`standardize_factors()` logic.

---

## New `R/gbcd.R`

Add:

```r
fit_gbcd_programs()
```

Responsibilities:

- optional dependency check;
- validate no missing data;
- validate `kmax >= 2`;
- orient sample × feature;
- memory preflight;
- call `gbcd::fit_gbcd`;
- validate output axes;
- convert log2 effects to natural-log reconstruction scale;
- construct conservative mixed-constraint representation;
- build signature-statistics table;
- attach provenance.

Do not copy upstream implementation functions.

---

## New `R/multistart.R`

Add:

```r
fit_multistart_ebmf()
```

Recommended helpers:

```r
.multistart_seeds()
.native_sign_constraints()
.multistart_init_fn()
.candidate_elbo()
```

Important:

- custom initializer must explicitly preserve sign constraints;
- candidate failures are recorded, not fatal unless all fail.

---

## `R/metadata.R`

Extend:

```r
summarize_factor_activity()
```

with:

```text
activity_threshold
activity_direction
active_n
active_proportion
active_weighted_proportion
```

No p-values.

---

## `R/plots.R`

Extend:

```r
plot_factor_activity()
```

with:

```r
order_within_group = c("input", "mds1")
```

Only apply MDS ordering for grouped structure plots.

---

## `R/runtime.R`

### Provenance engine metadata

Ensure newly exported wrappers identify their numerical engine correctly.

### Dependency reporting

When provenance engine/result is GBCD, include:

```text
gbcd
ashr
```

where relevant, in addition to the current core dependencies.

Do not globally require optional packages.

---

## `DESCRIPTION`

Add:

```text
gbcd
```

to Suggests.

Add pinned remote:

```text
stephenslab/gbcd@670e3afb50236303ddf99db33a3394abd4a711bd
```

Retain the existing flashier pin.

---

## `.github/workflows/R-CMD-check.yaml`

Full-engine job:

- install GBCD optional dependency;
- verify namespace loads;
- keep default check job functional with `_R_CHECK_FORCE_SUGGESTS_=false`.

Add tests that do not require GBCD to the ordinary job.

GBCD integration tests may use skips unless full-engine job is active.

---

## `inst/provenance.csv`

Add/update rows for:

```text
fit_seminonnegative_ebmf
fit_initialized_ebmf
fit_gbcd_programs
fit_multistart_ebmf
fit_ebmf intercept support
projection manifest/offset adaptation
```

Use exact source pins from `UPSTREAM_AUDIT.md`.

---

## Credits / notices

Review/update:

```text
CREDITS.md
vignettes/credits.Rmd
inst/NOTICE
```

GBCD package code is not copied; this is call-only integration.

DRIFT remains conceptual inspiration only.

---

## Tests

Expected new/modified files:

```text
tests/testthat/test-fit.R
tests/testthat/test-operations.R
tests/testthat/test-representation.R
tests/testthat/test-intervals.R
tests/testthat/test-projection.R
tests/testthat/test-resampling.R
tests/testthat/test-metadata.R
tests/testthat/test-plots.R
tests/testthat/test-gbcd.R
tests/testthat/test-multistart.R
```

---

## Documentation

Add one advanced fitting vignette/article covering:

```text
ordinary EBMF
full NMF
semi-NMF
feature intercept/offset
initialized factors
multistart optimization sensitivity
GBCD interoperability
```

Explicitly separate:

```text
model constraint
coordinate representation
initialization
optimization
projection
validation
external engine
```
