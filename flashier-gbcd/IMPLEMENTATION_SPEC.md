# Implementation specification

## 1. Non-negotiable design constraints

All implementation work must preserve these package invariants.

### 1.1 Numerical ownership

`flashier-utils` delegates numerical algorithms to their maintained packages.

Do not reimplement:

- FLASH greedy fitting;
- FLASH backfitting;
- EBNM solvers;
- softImpute;
- IRLBA;
- GBCD covariance decomposition.

### 1.2 Semantic orientation

Public package semantics are always:

```text
samples × factors = activity
features × factors = effects
```

`sample_side` maps those semantic axes to native FLASH loadings/factors.

### 1.3 Representation ownership

`factor_representation` is the interpretation/transfer object:

```r
list(
  activity,
  effects,
  offset,
  factors,
  manifest
)
```

A native FLASH intercept or external fixed baseline must be represented through
`offset`, not as an ordinary biological factor.

### 1.4 Provenance and resampling

No automatic stability refit may silently execute a fitting procedure different
from the one that created the reference result.

If the recipe is not reproducible through the package's standard automatic
callback:

```r
refit_requires_callback = TRUE
```

### 1.5 No hidden preprocessing

All normalization, transformation, scaling and GBCD input preparation remain
caller-controlled and recorded.

---

# 2. M1 — Foundation invariants

M1 must be completed before adding new fitting wrappers.

## 2.1 Correct `.prior_support()`

Replace the current stale/nonexistent family assumptions with audited current
`ebnm` exports.

At minimum:

```r
signed <- c(
  "ebnm_normal",
  "ebnm_point_normal",
  "ebnm_point_laplace",
  "ebnm_unimodal",
  "ebnm_unimodal_symmetric",
  "ebnm_normal_scale_mixture",
  "ebnm_horseshoe",
  "ebnm_npmle",
  "ebnm_flat",
  "ebnm_ash"
)

nonnegative <- c(
  "ebnm_point_exponential",
  "ebnm_unimodal_nonnegative",
  "ebnm_generalized_binary"
)
```

Remain conservative for:

```text
ebnm_unimodal_nonpositive
ebnm_group
ebnm_point_mass
ebnm_deconvolver
arbitrary closures
```

unless their orientation/sign behavior is explicitly audited.

Remove stale names not exported by the audited `ebnm` version, including:

```text
ebnm_point_gamma
ebnm_exponential
```

Do not invent a fourth support state in this milestone.

## 2.2 Minimal fitting-recipe metadata

Add to wrapped native fits:

```r
fit_recipe = list(
  method = "greedy",
  intercept = "none",
  initialization = "native_greedy",
  selection = NULL
)

refit_requires_callback = FALSE
```

This is a safety declaration, not a general workflow language.

## 2.3 Systematic automatic-refit policy

Audit every native-fit-producing path.

### Automatically reproducible

```text
fit_ebmf() with recognized priors and standard native initialization
fit_nonnegative_ebmf()
fit_seminonnegative_ebmf()
fit_ebmf(intercept = "feature") using the standard recipe
reorder_factors() when it is only a permutation
```

### Explicit callback required

Set:

```r
refit_requires_callback = TRUE
```

for at least:

```text
fit_initialized_ebmf()
selected fit from fit_multistart_ebmf()
project_factor_activity() native fit
project_factor_features() native fit
refit_factors()
remove_factors()
fix_factors()
unfix through fix_factors()
```

Reason: their current numerical/model state is not reproduced by the original
`fit_ebmf()` configuration alone.

Also require an explicit callback when:

```r
control$greedy$init_fn
```

is caller-supplied.

Custom/parameterized EBNM priors already trigger the existing
`custom_callback_required` behavior and should retain that rule.

## 2.4 Stability guard

Before automatic prior/configuration reconstruction in `factor_stability()`:

```r
if (isTRUE(native_metadata$refit_requires_callback) && is.null(fit_fn)) {
  stop("This fitting recipe requires an explicit fit_fn for stability refitting")
}
```

Run this guard before accessing `prior_identifiers` or `configuration`.

## 2.5 Projection manifest repair

This is a prerequisite for intercept support.

A projected native fit is a **new native fit**. Its semantic representation may
use values inherited from a source basis, but its native factor indices,
multipliers and PVE belong to the projection fit.

Current source metadata must not be reused blindly.

### Preserve source provenance separately

Before rebuilding projection metadata, retain:

```r
source_basis_id
source_factor_ids
source_native_indices
source_scaling
source_orientation
```

### Projection result representation

For the returned native projection fit:

```r
manifest$native_indices
```

must index the returned native fit, not the source fit.

Unless the projection native fit itself has an intercept:

```r
native_indices = seq_len(K)
```

Projection activity/effect matrices are native posterior coordinates of that
projection fit because the fixed side is initialized at exactly those supplied
values.

Therefore reset:

```r
scaling = "raw"
ldf_type = NULL
d_location = NULL
activity_multiplier = rep(1, K)
effects_multiplier = rep(1, K)
orientation = "as_fit"
orientation_sign = rep(1, K)
orientation_anchor = rep(NA_character_, K)
```

and set:

```r
native_pve = projection_native_fit$pve[native_indices]
```

Do not copy source PVE into a projection result.

## 2.6 Fixed external feature offsets on native projection fits

`project_factor_activity()` estimates activity on:

```text
query - source offset
```

so its returned native fit models the residual matrix, not the complete
measurement matrix.

To allow `standardize_factors(result$fit)` to recover the complete semantic
representation, attach to its `flashbridge_metadata`:

```r
feature_offset = source$offset
offset_definition = "fixed_external_projection_offset"
```

and set:

```r
refit_requires_callback = TRUE
```

`standardize_factors()` must use `feature_offset` when there is no native
offset-factor declaration.

Do not add these metadata fields to ordinary zero-offset fits; preserve existing
basis hashes for unchanged old workflows.

---

# 3. M2 — Semi-nonnegative EBMF

## 3.1 Public API

```r
fit_seminonnegative_ebmf(
    X,
    sample_side,
    max_factors = 50L,
    activity_prior = ebnm::ebnm_point_exponential,
    effect_prior = ebnm::ebnm_point_laplace,
    initialization = "native_greedy",
    seed = NULL,
    control = list(),
    ...
)
```

Only `"native_greedy"` is supported initially.

## 3.2 Semantic model

```text
activity >= 0
effects signed
```

### Samples in rows

```r
ebnm_fn <- list(activity_prior, effect_prior)
```

### Samples in columns

```r
ebnm_fn <- list(effect_prior, activity_prior)
```

Delegate to `fit_ebmf()`.

## 3.3 Validation

- reject `ebnm_fn` supplied through `...`;
- require `.prior_support(activity_prior) == "nonnegative"`;
- effect prior may be `"signed"` or `"unknown"`;
- do not copy GBCD split-signed initialization;
- do not transform X.

A parameterized activity-prior closure whose support cannot be verified should
be fitted through explicit lower-level `fit_ebmf()` rather than falsely
certified by this convenience wrapper.

## 3.4 Metadata

Update the returned native metadata:

```r
fit_recipe$method <- "semi_nonnegative"
fit_recipe$initialization <- "native_greedy"
```

For default priors:

```r
prior_support == c(
  activity = "nonnegative",
  effects = "signed"
)
```

in both orientations.

---

# 4. M3 — Native intercept → semantic offset

## 4.1 Public fit API

Extend:

```r
fit_ebmf(
    ...,
    intercept = c("none", "feature")
)
```

Do not expose custom intercept priors in v1.

`"feature"` means:

```text
feature-specific baseline shared across all samples
```

and uses public:

```r
flashier::flash_add_intercept()
```

## 4.2 Preserve the old default path

When:

```r
intercept = "none"
```

leave the current direct `flashier::flash()` fast path untouched wherever it is
currently used.

This is required for numerical regression compatibility.

## 4.3 Intercept-enabled pipeline

```r
object <- flashier::flash_init(...)

object <- flashier::flash_add_intercept(
  object,
  rowwise = sample_side == "rows"
)

if (max_factors > 0L) {
  object <- flashier::flash_greedy(
    object,
    Kmax = max_factors,
    ebnm_fn = ebnm_fn,
    ...
  )
}

if (backfit && object$n_factors > 0L) {
  object <- flashier::flash_backfit(object, ...)
}

if (nullcheck && object$n_factors > 1L) {
  adaptive <- setdiff(seq_len(object$n_factors), 1L)
  object <- flashier::flash_nullcheck(
    object,
    kset = adaptive,
    ...
  )
}
```

### Critical edge case

`max_factors = 0L` is valid.

It produces an intercept-only model and must **not** call:

```r
flash_greedy(Kmax = 0)
```

because upstream requires `Kmax >= 1`.

### Control protection

For intercept-enabled fitting, caller-supplied:

```r
control$nullcheck$kset
```

must be rejected so the fixed intercept cannot accidentally enter the
nullcheck set.

Do not break existing no-intercept behavior solely to enforce this rule.

## 4.4 Native metadata

After fitting and adaptive nullchecking:

```r
factor_ids <- c(
  "Intercept",
  paste0("F", seq_len(object$n_factors - 1L))
)

factor_roles <- c(
  "offset",
  rep("factor", object$n_factors - 1L)
)

offset_factor_index <- 1L
```

Store:

```r
configuration$intercept <- "feature"

fit_recipe$intercept <- "feature"
```

The intercept is always kept as native factor 1 in v1.

## 4.5 Standardization algorithm

### Step 1 — obtain raw semantic native means

Use `.raw_pair()` and map native dimensions according to `sample_side`.

### Step 2 — derive the offset from the RAW intercept

Do this **before LDF normalization**.

Validate:

```r
semantic_activity[, offset_factor_index]
```

is all ones within strict numerical tolerance.

Set:

```r
offset <- semantic_effects[, offset_factor_index]
```

The LDF-normalized intercept column must never define the offset because LDF
rescales both sides.

### Step 3 — build the adaptive factor set

```r
adaptive_native_indices <- setdiff(
  seq_len(native_K),
  offset_factor_index
)
```

Remove the native intercept from both semantic factor matrices.

### Step 4 — LDF if requested

When `scaling = "ldf"`:

- compute native LDF for the whole fit;
- subset only `adaptive_native_indices`;
- place D on the requested adaptive side;
- compute multipliers only for adaptive factors.

The raw offset is unchanged by LDF scaling.

### Step 5 — manifest

Set:

```r
manifest$native_indices <- adaptive_native_indices
manifest$native_pve <- fit$pve[adaptive_native_indices]
manifest$offset_definition <- "native_flash_intercept"
```

Add offset-specific manifest fields only for offset-bearing fits so old
zero-offset basis hashes remain unchanged.

## 4.6 Semantic PVE

For an intercept-bearing wrapped fit:

```r
factor_pve(fit)
factor_pve(standardize_factors(fit))
```

should both report **adaptive semantic factors only**.

Do not expose `"Intercept"` as an ordinary factor through `factor_pve()`.

If native offset PVE is useful for diagnostics, record it separately in
metadata rather than mixing it into the factor table.

## 4.7 Native-index consumers

All native posterior consumers must use:

```r
view$manifest$native_indices
```

instead of assuming semantic column k equals native column k.

### `factor_uncertainty()`

- verify native row axis matches semantic IDs;
- subset native posterior matrices by `native_indices`;
- then apply representation multipliers.

### `factor_intervals()`

Sampler draw matrices retain all native factors.

Do not require:

```r
dim(draw[[native_side]]) == dim(view[[side]])
```

when native offset factors exist.

Instead validate:

```text
same observation/feature row count
native factor count covers max(native_indices)
```

then index sampler columns through `native_indices`.

## 4.8 Stability semantics

For each replicate:

```r
view <- standardize_factors(fitted, ...)
```

Use:

```r
rank = ncol(view$activity)
pve = factor_pve(view)
```

not native `fitted$n_factors` or `factor_pve(fitted)`.

The offset is not a biological program and must not inflate semantic rank.

## 4.9 Editing invariants

For fits with `offset_factor_index = 1L`:

### `remove_factors()`

Reject removal of `"Intercept"`.

### `fix_factors()`

Reject generic fixing/unfixing of `"Intercept"`.

### `reorder_factors()`

Allow adaptive factors to be reordered, but require:

```r
order[1] == "Intercept"
```

so the native offset remains factor 1 in v1.

### `refit_factors()`

May backfit the intercept and adaptive factors because the all-ones side remains
natively fixed.

The result is nevertheless marked:

```r
refit_requires_callback = TRUE
```

because its post-hoc refit recipe is not represented by the original automatic
fit configuration.

## 4.10 Offset-aware projection of new features

`project_factor_activity()` is already offset-aware for fixed source effects.

`project_factor_features()` is not.

When the source representation has a nonzero offset, new features require a new
feature-specific baseline.

Do not silently return zero offsets.

### Required algorithm

For a nonzero-offset source:

1. create the new-feature native FLASH object;
2. add a native feature intercept using `flash_add_intercept()` with semantic
   orientation;
3. obtain residuals after this baseline initialization;
4. initialize biological effects from fixed source activities against those
   residuals;
5. initialize biological factors after the native intercept;
6. fix only the biological activity side;
7. jointly backfit the intercept effect and biological effects;
8. extract native factor 1 as the new feature offset;
9. return factors 2..K+1 as biological factors.

The result manifest uses:

```r
native_indices = 2:(K + 1L)
offset_definition = "native_flash_intercept"
```

For zero-offset sources, keep the current branch to preserve existing numerical
behavior.

---

# 5. M4 — Initialized / anchored EBMF

## 5.1 Public API

```r
fit_initialized_ebmf(
    X,
    sample_side,
    activity_init,
    effect_init,
    activity_prior = ebnm::ebnm_point_normal,
    effect_prior = ebnm::ebnm_point_laplace,
    fixed_side = c("none", "activity", "effects"),
    add_factors = 0L,
    nullcheck_new = TRUE,
    S = NULL,
    var_type = 0L,
    id_policy = "require",
    feature_scale = "unknown",
    independent_unit = NULL,
    seed = NULL,
    verbose = 1L,
    control = list()
)
```

Intercept + initialized fitting is deferred.

The helper always performs an initial backfit; raw unoptimized initialization is
not exposed as a completed fit.

## 5.2 Required initialization matrices

```text
activity_init: samples  × K_init
effect_init:   features × K_init
```

Require:

- exact sample IDs and order-independent matching;
- exact feature IDs;
- equal factor counts;
- identical factor names and order;
- nonempty unique factor names;
- finite values;
- `K_init >= 1`.

Reorder rows to the validated X axes before native fitting.

## 5.3 Prior-support consistency

If an initialized side uses a recognized nonnegative prior, reject negative
initial values on that side.

Do not claim support for a custom prior whose support is unknown.

## 5.4 Native ordering

```r
native_init <- if (sample_side == "rows") {
  list(activity_init, effect_init)
} else {
  list(effect_init, activity_init)
}

native_priors <- if (sample_side == "rows") {
  list(activity_prior, effect_prior)
} else {
  list(effect_prior, activity_prior)
}
```

## 5.5 Fixed-side mapping

```text
semantic activity:
  rows    -> loadings
  columns -> factors

semantic effects:
  rows    -> factors
  columns -> loadings
```

Only one semantic side may be fixed.

## 5.6 Pipeline

```text
flash_init
  -> flash_factors_init
  -> optional flash_factors_fix on initialized factors
  -> mandatory initial flash_backfit
  -> optional flash_greedy(add_factors)
  -> final flash_backfit when new factors were added
  -> optional nullcheck of NEW factors only
```

Initialized factors are never automatically nullchecked.

## 5.7 Stable IDs

Preserve all initialized factor IDs.

Generate added factor IDs deterministically by scanning:

```text
F1, F2, ...
```

and skipping collisions with user names.

## 5.8 Nullcheck while preserving initialized IDs

If new factors are nullchecked with `remove = TRUE`, do not nullcheck all new
indices in one opaque call and then guess which IDs survived.

Process new factors from highest native index downward.

For each new factor:

1. record native factor count;
2. nullcheck that single native index;
3. if factor count drops, remove exactly that generated ID;
4. continue downward.

This preserves initialized indices and user names.

## 5.9 Refit policy

Store initialization matrix digests and factor IDs, but not full matrices.

Set:

```r
fit_recipe$method <- "initialized"
fit_recipe$initialization <- "user_matrices"
refit_requires_callback <- TRUE
```

Digests establish identity; they do not recreate the fitting inputs.

---

# 6. M5 — GBCD interoperability

## 6.1 Dependency pin

Add optional dependency:

```text
Suggests: gbcd
```

and GitHub remote:

```text
stephenslab/gbcd@670e3afb50236303ddf99db33a3394abd4a711bd
```

Do not import its symbols into NAMESPACE.

Update full-engine CI to install/load it.

## 6.2 Public API

Use GBCD's own terminology:

```r
fit_gbcd_programs(
    X,
    sample_side,
    kmax,
    feature_scale,
    preprocessing,
    prior = ebnm::ebnm_generalized_binary,
    dense_budget_bytes = 2 * 1024^3,
    keep_native = FALSE,
    ...
)
```

Do **not** call this parameter `max_factors`.

Upstream documents:

```text
Kmax >= 2
final rank often near Kmax
final rank can be as high as 2*Kmax - 1
```

Validate:

```r
kmax >= 2L
```

## 6.3 Input requirements

Use:

```r
X <- .validate_matrix(X, allow_missing = FALSE)
```

GBCD input is explicitly caller-prepared normalized/log-transformed
single-cell expression.

Require:

- non-unknown `feature_scale`;
- explicit `preprocessing` list with nonempty scope/description;
- no hidden normalization;
- no missing values.

Orient without unnecessary densification:

```r
Y <- if (sample_side == "rows") X else Matrix::t(X)
```

## 6.4 Memory preflight

The audited GBCD implementation may either use a low-rank-plus-sparse
representation or materialize a cell × cell covariance matrix.

Before fitting:

1. estimate the same sparsity branch used by the pinned GBCD implementation;
2. record the predicted strategy;
3. if the dense covariance branch is predicted, estimate at least:

```r
8 * nrow(Y)^2
```

bytes plus a safety factor;
4. refuse execution if this exceeds `dense_budget_bytes`.

The preflight is tied to the pinned GBCD source and must be updated if that pin
changes.

## 6.5 Delegate fitting

Call only:

```r
gbcd::fit_gbcd(
  Y = Y,
  Kmax = kmax,
  prior = prior,
  ...
)
```

Do not copy:

```text
fit_ebmf_to_YY
fit_ebmf_to_Y
init_cov_ebnmf
```

from the source package/workflow.

## 6.6 Output validation

Require:

```text
res$L                  sample × component
res$F$lfc              feature × component
res$F$z_score          same dimensions as lfc
res$F$lfsr             same dimensions as lfc
identical factor axes
```

Set/verify row names against the input axes.

Reject malformed upstream output rather than guessing.

## 6.7 Reconstruction scale

GBCD returns feature effects as log2 fold changes.

Use:

```r
A <- as.matrix(res$L)
B_log2 <- as.matrix(res$F$lfc)
B <- log(2) * B_log2
```

so:

```r
tcrossprod(A, B)
```

matches the GBCD reconstruction scale.

## 6.8 Mixed activity support

Do not set global activity support to `"nonnegative"`.

The audited implementation treats:

```text
Baseline -> signed point-Laplace covariance component
GEPs     -> nonnegative membership prior
```

Therefore set conservative global representation support:

```r
constraints = c(
  activity = "unknown",
  effects = "signed"
)
```

and add component-level manifest metadata:

```r
factor_constraints <- data.frame(
  factor = factor_ids,
  activity_support = c("signed", rep("nonnegative", K - 1L)),
  effects_support = rep("signed", K)
)

factor_roles <- c(
  "baseline",
  rep("program", K - 1L)
)
```

The returned `"Baseline"` remains a factor. It is not an all-ones offset.

## 6.9 Representation

Build:

```r
factor_representation(
  activity = A,
  effects = B,
  offset = 0
)
```

with:

```text
engine = gbcd
native_indices = seq_len(K)
native_pve = NULL
feature_scale = caller declaration
preprocessing = caller declaration
```

## 6.10 Signature statistics

Return a separate long table:

```text
feature
factor
effect
log2_fc
z_score
lfsr
basis_id
```

Do not imply that `z_score` and `lfsr` are posterior uncertainty summaries for
the returned EB-SNMF `effect` estimate.

The upstream code uses two different stages:

```text
effect / log2_fc:
  EB-SNMF posterior mean

z_score / lfsr:
  gene-wise linear regression conditional on estimated L
  followed by ashr shrinkage
```

Record in `analysis_metadata`:

```text
effect_source = gbcd_eb_snmf
z_score_definition = ashr_posterior_mean / posterior_sd
lfsr_source = ashr
uncertainty_conditioning = conditional_on_estimated_memberships
```

## 6.11 Result object

Return a package-owned result:

```r
structure(
  list(
    representation = view,
    feature_statistics = table,
    native = if (keep_native) res else NULL,
    analysis_metadata = metadata
  ),
  class = "gbcd_program_result"
)
```

No generic external-factor S3 framework is introduced yet.

## 6.12 Native-only operations

A GBCD representation does not support native FLASH posterior operations:

```text
factor_uncertainty
factor_intervals
native factor PVE
```

`factor_pve(view)` must clearly report native PVE unavailable.

Extraction, ranking, metadata summaries, plotting and enrichment may operate on
the representation normally.

## 6.13 Provenance dependency reporting

`provenance(..., include_dependencies = TRUE)` must include `gbcd` when the
result engine is GBCD, along with its runtime version/RemoteSha when available.

---

# 7. M6 — Group prevalence summaries

Extend:

```r
summarize_factor_activity()
```

with:

```r
activity_threshold = NULL
activity_direction = c("positive", "negative", "absolute")
```

A single nonnegative scalar threshold is enough for v1.

The threshold is applied to the resolved representation coordinates using a
strict:

```text
>
```

operator.

No display scaling is performed.

When a threshold is present, add:

```text
active_n
active_proportion
```

per factor × group.

If observation weights are supplied, additionally return:

```text
active_weighted_proportion
```

defined as the weighted mean of the binary active indicator.

Retain:

```text
active_n
active_proportion
```

as raw observation summaries so weighted and unweighted quantities are never
confused.

Record the representation basis/scaling in analysis metadata.

---

# 8. M7 — Multistart optimization sensitivity

## 8.1 Public API

The first implementation covers ordinary non-intercept, non-initialized EBMF.

```r
fit_multistart_ebmf(
    X,
    sample_side,
    ebnm_fn = ebnm::ebnm_point_normal,
    S = NULL,
    var_type = 0L,
    max_factors = 50L,
    backfit = FALSE,
    nullcheck = TRUE,
    id_policy = "require",
    feature_scale = "unknown",
    independent_unit = NULL,
    n_starts = 10L,
    init_method = c("native_random", "softImpute", "irlba"),
    init_seeds = NULL,
    seed,
    select = c("elbo", "none"),
    keep_fits = TRUE,
    verbose = 1L,
    control = list()
)
```

Require either:

- explicit unique `init_seeds`, or
- a master `seed` from which unique initializer seeds are deterministically
  generated.

Require `n_starts >= 2`.

## 8.2 Native prior ordering

Normalize:

```r
native_priors <- if (is.list(ebnm_fn)) ebnm_fn else rep(list(ebnm_fn), 2L)
```

and require exactly two functions.

## 8.3 Constraint-preserving native random starts

Overriding `flash_greedy(init_fn=...)` disables the native initializer wrapper
that normally injects sign constraints.

Therefore derive native sign constraints explicitly:

```text
support == signed       -> 0
support == nonnegative  -> +1
support == unknown      -> unsupported for native_random multistart
```

Then each start uses:

```r
init_fn <- local({
  this_seed <- init_seed
  this_constraints <- sign_constraints
  function(f) {
    flashier::flash_greedy_init_default(
      f,
      sign_constraints = this_constraints,
      seed = this_seed
    )
  }
})
```

This is essential for NMF/semi-NMF-compatible priors.

## 8.4 softImpute / IRLBA restrictions

The upstream softImpute and IRLBA initializer helpers do not accept sign
constraints.

For v1:

```text
softImpute
irlba
```

are allowed only when both native prior dimensions are audited `"signed"`.

Also enforce upstream restrictions:

```text
softImpute:
  not lowrank/lrps input

irlba:
  no missing values
  intended for dense input
```

## 8.5 Control conflict

Reject caller-supplied:

```r
control$greedy$init_fn
```

because multistart owns initialization.

Preserve compatible remaining greedy/backfit/nullcheck controls.

## 8.6 Failure-aware execution

Each candidate records:

```text
start
init_seed
status
reason
rank
elbo
fit_id
elapsed_seconds
selected
```

A failed start does not abort all other starts.

If every start fails, error with the candidate failure summary attached or
returned in a structured failure object.

Do not make summed factor PVE a required selection statistic; it is not a
unique variance decomposition.

## 8.7 ELBO

Extract only through public:

```r
flashier::flash_fit_get_elbo(
  flashier::flash_fit(fit)
)
```

`select = "elbo"` chooses the maximal finite ELBO among successful starts.

Documentation must state:

```text
ELBO ranks optimization solutions under the fitted model.
It does not define biological interpretability.
```

## 8.8 Refit policy

The selected fit represents the **multistart selection procedure**, not one
ordinary greedy run.

Set:

```r
fit_recipe$method <- "multistart"
fit_recipe$selection <- "max_elbo"
refit_requires_callback <- TRUE
```

until `factor_stability()` can explicitly rerun the full multistart procedure.

---

# 9. M8 — STRUCTURE-style within-group ordering

Extend:

```r
plot_factor_activity()
```

with:

```r
order_within_group = c("input", "mds1")
```

Default `"input"` preserves current behavior.

`"mds1"` is valid only when:

```text
style == "structure"
group is supplied
```

Algorithm:

1. perform existing nonnegative check;
2. remove zero-sum rows;
3. row-sum normalize;
4. perform existing bounded sample subsampling;
5. determine group order by first appearance among displayed metadata rows;
6. for each group:
   - if fewer than 2 displayed samples, keep input order;
   - compute Euclidean distance on displayed factor proportions;
   - if all distances are zero, keep input order;
   - otherwise run `stats::cmdscale(..., k = 1)`;
   - order by MDS coordinate with original input position as deterministic tie
     breaker;
7. concatenate groups in the preserved group order;
8. use that vector as the sample factor levels in the plot.

Record:

```text
order_within_group
displayed_sample_order
group_order
ordering_distance = euclidean
ordering_embedding = cmdscale_k1
```

This remains display-only.

---

# 10. Deferred experimental work

## 10.1 Split-signed semi-NMF initialization

Do not add to the public API during M1–M8.

Benchmark only after native semi-NMF and multistart are implemented.

Compare:

```text
ELBO
held-out reconstruction where scientifically valid
factor congruence
biological-unit stability
runtime
program interpretability summaries
```

Promote only on reproducible benefit.

## 10.2 General external-factor adapter

Do not introduce a generic S3 conversion framework solely for GBCD.

Revisit after at least one second external engine needs the same abstraction.

## 10.3 Arbitrary prior serialization

Do not reverse-engineer arbitrary R closures.

Continue requiring explicit callbacks for unreconstructable priors.

---

# 11. Backward compatibility

The following default calls must remain numerically/semantically unchanged:

```text
fit_ebmf(..., intercept = "none")
fit_nonnegative_ebmf(...)
standardize_factors() on ordinary no-offset fits
factor_uncertainty() on ordinary fits
factor_intervals() on ordinary fits
plot_factor_activity(..., order_within_group = "input")
```

Existing ordinary-fit basis IDs should remain unchanged.

Projection result manifests/basis IDs **may intentionally change** because M1
repairs incorrect native-index/PVE/multiplier provenance in those results.
This must be called out in NEWS as a correctness fix.

---

# 12. Release rule

A feature is complete only when:

1. direct public-upstream equivalence tests pass;
2. both sample orientations pass where applicable;
3. semantic/native factor-index mapping is correct;
4. offset reconstruction identity passes;
5. resampling either reproduces the fitting procedure or refuses automatically;
6. provenance identifies the original numerical engine;
7. existing default workflows remain green.
