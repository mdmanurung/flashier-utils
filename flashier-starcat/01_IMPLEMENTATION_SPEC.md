# flashier-utils Consensus Reference and Projection Layer
## Revised implementation specification

**Target repository:** `mdmanurung/flashier-utils`  
**Primary conceptual reference:** Kotliar et al., *Nature Methods* (2025), DOI `10.1038/s41592-025-02793-1`  
**Inference engine:** `flashier`  
**EBNM engine:** `ebnm`

---

# 1. Product goal

Extend `flashier-utils` with a reusable cross-dataset factor-reference workflow that can:

1. accept independently fitted `flashier` factors from multiple datasets;
2. compare feature programs without changing the source inferential bases;
3. identify reproducible cross-dataset program groups;
4. construct an auditable consensus `factor_reference`;
5. project the fixed reference programs into a novel dataset using native `flashier` EBMF machinery;
6. quantify feature coverage, reconstruction and reference support;
7. optionally discover **query-residual programs** after accounting for the fixed reference, if a dedicated native-feasibility probe succeeds.

The design borrows the **independent discovery -> cross-dataset consensus -> fixed-reference projection** pattern from STARCAT while retaining `flashier`'s signed/semi-nonnegative factor models and posterior/shrinkage machinery.

---

# 2. Scope boundaries

## Core release

Core release includes:

- cross-dataset comparison profiles;
- all-pairs factor similarity;
- deterministic reference graph and clustering;
- consensus program estimation;
- `factor_reference` object and accessors;
- native EBMF fixed-reference projection;
- support, coverage and reconstruction diagnostics;
- provenance, tests and benchmarks.

## Conditional extension

Residual-program discovery is implemented only after the native feasibility milestone proves that reference feature programs remain fixed through greedy factor addition and backfitting.

## Experimental / deferred

The following are explicitly not required for the core release:

- hierarchical or tensor EBMF;
- joint multi-dataset factorization;
- batch correction or dataset integration;
- automatic biological annotation/classification;
- learned transportable consensus activity priors;
- frozen residual-noise parameters;
- a second production NNLS implementation;
- automatic singleton rescue based on biological interpretation;
- consensus offsets/intercepts.

---

# 3. Existing package contracts that must be preserved

The implementation must reuse and preserve current behavior of:

- `standardize_factors()`;
- `factor_replicability()` and `.align_programs()`;
- `factor_stability()`;
- `project_factor_activity()`;
- `project_factor_features()`;
- `.native_projection()`;
- representation `basis_id` integrity;
- explicit `feature_scale` declarations;
- preprocessing manifests;
- one-to-one feature mapping rules;
- dense-allocation guards;
- native warning retention;
- provenance registry and credit requirements.

The current source `factor_representation` remains the authoritative inferential representation of each independently fitted model.

---

# 4. Architecture

## Layer A — source factor representations

For dataset `d`:

\[
X_d \approx A_d B_d^T
\]

where `A_d` is activity and `B_d` is feature effect.

Source inputs may be native `flash` fits or `factor_representation` objects. Native fits are resolved with:

```r
standardize_factors(fit, scaling = "raw")
```

unless the caller explicitly provides a representation.

Source objects are never mutated.

## Layer B — comparison profiles

Comparison profiles are normalized feature vectors used only for cross-dataset matching. They are **not** returned as an inferential `factor_representation` and do not replace source activity/effect coordinates.

Internal conceptual object:

```r
list(
  profiles = <feature x member matrix or pairwise accessor>,
  members = <member metadata>,
  feature_sets = <per-dataset and common sets>,
  transform = <normalization/standardization manifest>,
  exclusions = <factor/feature exclusion records>
)
```

## Layer C — factor reference

A `factor_reference` contains final consensus programs on an explicitly defined reference feature universe plus the graph, memberships, diagnostics and projection contract.

## Layer D — reference projection

`project_reference()` fixes the consensus feature programs and estimates activity in a new query using the same native fixed-factor pattern already used by `project_factor_activity()` and `.native_projection()`.

---

# 5. Public API

## Core new exports

```r
build_factor_reference()
project_reference()
reference_programs()
reference_members()
reference_diagnostics()
print.factor_reference()
summary.factor_reference()
```

## Conditional export after feasibility milestone

```r
discover_residual_programs()
```

## Optional plotting exports

Add only if the corresponding diagnostic cannot be expressed cleanly through existing plotting utilities:

```r
plot_reference_similarity()
plot_reference_members()
plot_reference_program()
```

## Internal, not initially exported

```r
.reference_collect_sources()
.reference_make_profiles()
.reference_pair_features()
.reference_similarity_table()
.reference_build_graph()
.reference_cluster_graph()
.reference_estimate_consensus()
.validate_factor_reference()
.reference_projection_view()
```

This keeps the public surface substantially smaller than the first bundle.

---

# 6. Source eligibility and comparability

Before computing any similarity, build a member eligibility table.

Required columns:

```text
member_id
dataset_id
factor_id
basis_id
feature_scale
activity_support
effect_support
n_features
zero_activity
zero_effect
eligible
reason
```

## Schema-v1 hard requirements

A source factor is ineligible when:

- its feature effect vector is identically zero;
- its feature IDs are not unique;
- its feature scale is `unknown`;
- source offsets are nonzero;
- its declared constraints are internally inconsistent;
- it lacks enough comparable features for the requested pairwise operation.

Zero or otherwise unusable factors should be **recorded and excluded**, not silently discarded.

Global incompatibilities should fail the build, including incompatible feature semantics or impossible feature mappings.

---

# 7. Comparison contract

`build_factor_reference()` must require an explicit comparison contract whenever more than one dataset is supplied.

Recommended argument:

```r
comparison_manifest = list(
  mode = "generic",
  feature_scale = "...",
  preprocessing_compatibility = "user_declared",
  profile_normalization = "auto",
  feature_standardization = "none"
)
```

Allowed `mode` values:

```text
generic
starcat2025
```

The manifest is stored verbatim after validation and included in `reference_id` hashing.

`feature_scale` must agree with every eligible source representation.

The package does not infer that two differently processed matrices are biologically or numerically comparable simply because both are called `logcounts`, `TP10K`, etc.

---

# 8. Generic comparison profiles

For every member factor, create a profile on the requested feature set.

## Non-negative effects

If feature-effect support is declared non-negative:

\[
p_k = B_k / \sum_j B_{jk}
\]

after subsetting to the relevant feature universe.

Requirements:

- all retained effects must be non-negative within numerical tolerance;
- the retained sum must be positive;
- negative values outside tolerance fail rather than being clipped.

## Signed effects with signed activity

If both activity and effect supports are signed:

\[
p_k = B_k / ||B_k||_2
\]

The profile may be sign-flipped for comparison because flipping both source sides preserves the model constraints. Orientation is deterministic using:

1. explicit supplied anchor if valid;
2. largest absolute-effect feature;
3. lexical tie break.

## Semi-nonnegative factors

If activity is non-negative and feature effects are signed, sign is **not arbitrary**. Do not sign-flip the feature vector. Normalize with L2 norm while retaining the fitted sign.

## Unknown support

Unknown support may be compared without sign flipping. Operations that require deciding sign equivalence must fail closed.

---

# 9. STARCAT-2025 comparison mode

`mode = "starcat2025"` is a specialized compatibility mode for source programs that satisfy STARCAT-like assumptions.

It requires the caller to provide:

```r
source_feature_sets   # named list of discovery/overdispersed features per dataset
source_feature_sd     # named list of per-feature SD vectors in the declared linear expression scale
```

and to declare a non-negative linear program scale suitable for the requested transformation.

For a pair of source programs from datasets `a` and `b`:

1. restrict both spectra to features measured in both datasets;
2. renormalize each restricted non-negative spectrum to sum to `1e6`;
3. divide each spectrum feature-wise by the corresponding source-dataset SD;
4. restrict the standardized vectors to `union(source_feature_sets[[a]], source_feature_sets[[b]])`;
5. compute Pearson correlation.

Thus the final correlation feature set is:

```text
features measured in both datasets
INTERSECT
union(source_feature_sets[[a]], source_feature_sets[[b]])
```

All required source SD values must be finite and strictly positive on retained features.

The core generic EBMF workflow must **not** automatically apply this transformation to log-expression or signed factors.

---

# 10. Pairwise similarity

The reference graph needs all-pairs similarities, not one-to-one assignment.

Internal interface:

```r
.reference_similarity_table(
  sources,
  metric = c("congruence", "pearson", "cosine"),
  min_shared_features = 2L,
  allow_sign_flip = NULL,
  ...
)
```

Default generic metric: `congruence`, matching existing package semantics.

`starcat2025` mode forces Pearson correlation.

Output columns:

```text
member_i
dataset_i
factor_i
member_j
dataset_j
factor_j
raw_similarity
similarity
orientation
shared_features
coverage_i
coverage_j
eligible
reason
```

Rules:

- only one unordered row per pair;
- same-dataset similarities may be computed for diagnostics but are never eligible graph edges;
- `similarity = abs(raw_similarity)` only when sign equivalence is allowed by both source constraints;
- no LSAP assignment is applied here;
- feature coverage is retained for every comparison.

---

# 11. Reference graph

Internal graph construction receives the all-pairs table.

Conceptual interface:

```r
.reference_build_graph(
  similarity,
  min_similarity,
  rank_max,
  reciprocal = TRUE,
  rank_scope = "all_other",
  cross_dataset_only = TRUE
)
```

Return both:

```text
directed_candidates
accepted_edges
```

This separation is required for auditability.

## Generic defaults

Recommended generic defaults:

```text
min_similarity: caller-selected; no claim of universal biological threshold
rank_max: NULL = no rank filtering
reciprocal: FALSE unless requested
```

## STARCAT-2025 defaults

For `mode="starcat2025"`:

```text
metric = Pearson
min_similarity = 0.5 with strict > comparison
rank_max = number of source datasets when NULL
reciprocal = TRUE
rank_scope = all factors from other datasets
```

Directed-candidate construction follows the audited STARCAT implementation contract:

1. for each source member, rank all members from other datasets by similarity;
2. retain the first `rank_max` members overall;
3. discard candidates not strictly above `min_similarity`;
4. within the retained list, keep only the highest-ranked candidate from each target dataset;
5. create an accepted undirected edge only when the directed candidate relationship is reciprocal.

For exact TCAT reproduction, `rank_max=7`.

Tie order must be deterministic by similarity then member ID.

---

# 12. Deterministic consensus clustering

Internal interface:

```r
.reference_cluster_graph(
  graph,
  min_cluster_density = 2/3,
  max_members_per_dataset = 1L
)
```

Algorithm:

1. initialize every eligible member as its own set;
2. order accepted undirected edges by decreasing similarity;
3. consider merging the sets containing each edge endpoint;
4. reject a merge if the resulting set contains more than one member from any dataset;
5. otherwise compute density on the **accepted undirected graph**:

\[
density(C) = |E_C| / { |C| \choose 2 }
\]

6. merge only if the resulting density satisfies the declared threshold;
7. retain all unresolved members as singleton clusters.

For `starcat2025`, use `density >= 2/3`, matching the paper's "at least two-thirds" rule. This also matches the audited implementation's practical `> 0.666` criterion at exact two-thirds density. Freeze this choice in the algorithm version and parity tests.

The algorithm must be invariant to input row/factor/dataset ordering.

---

# 13. Reference feature universe

Pairwise comparison and final consensus estimation use different feature policies.

## Generic default

```text
reference_features = intersection of measured source effect features
```

The caller may supply an explicit subset.

## STARCAT-2025 policy

When `source_feature_sets` are supplied:

```text
union of source discovery/overdispersed feature sets
INTERSECT
features measured in every source dataset
```

This follows the published catalog-construction logic.

## Schema-v1 restriction

Do not construct a reference by averaging over a union with implicit zeros or missing-value imputation.

A future schema may add presence-threshold references, but v1 uses an explicit common feature basis.

---

# 14. Consensus program estimation

Internal interface:

```r
.reference_estimate_consensus(
  sources,
  clusters,
  reference_features,
  estimator = "mean"
)
```

## Default estimator

For each cluster:

1. regenerate each member comparison profile on the **final reference feature universe**;
2. align sign only when sign equivalence is permitted;
3. compute the element-wise arithmetic mean;
4. renormalize the final consensus:
   - L1 for non-negative programs;
   - L2 for signed/semi-nonnegative programs.

Do not average raw arbitrarily scaled EBMF effects.

## Constraint compatibility

All members in one consensus must have compatible effect and activity support.

Examples:

- all nonnegative/nonnegative -> valid;
- all nonnegative activity + signed effect -> valid;
- all signed/signed -> valid with allowed sign alignment;
- mixed nonnegative and signed activity support -> reject cluster as constraint-incompatible.

Rejected clusters remain in diagnostics with a reason; they do not enter the projection-ready program matrix.

## Singletons

All singletons remain in `clusters` and diagnostics.

By default, only clusters with:

```r
n_datasets >= min_datasets
```

enter `reference$programs`.

Default `min_datasets = 2L`.

Biologically justified singleton inclusion is an explicit caller decision, never an automatic package rule.

---

# 15. `factor_reference` class

```r
class(x) <- c("factor_reference", "list")
```

Required fields:

```r
list(
  programs,              # feature x reference-program matrix
  program_constraints,   # activity/effect support per reference program
  programs_table,
  members_table,
  feature_table,
  clusters,
  similarity,
  graph,
  excluded_members,
  feature_sets,
  comparison_manifest,
  projection_manifest,
  projection_ready,
  source_table,
  build_settings,
  reference_id,
  schema_version = "1.0.0",
  provenance
)
```

## No consensus offsets in schema v1

Reference offset is defined as zero on all reference features.

Source inputs with nonzero offsets are refused by the builder.

## Compact source provenance

Do not store source `flash` fits or complete sample-level manifests by default.

Store a compact `source_table` containing at least:

```text
dataset_id
basis_id
fit_id
feature_scale
preprocessing_hash
activity_support
effect_support
n_factors
n_features
```

Optional `keep_source_objects=TRUE` may be provided for debugging but is not the default serialized reference format.

---

# 16. Deterministic `reference_id`

`reference_id` is a SHA-256 digest of stable scientific content only:

- final programs;
- program constraints;
- membership table;
- reference feature IDs;
- source basis IDs;
- comparison manifest;
- clustering/consensus settings;
- projection manifest when present;
- schema version.

Exclude:

- wall-clock timestamps;
- object memory addresses;
- non-deterministic ordering;
- transient warnings.

Runtime/package versions remain in provenance but need not alter the scientific content hash unless they change an explicitly versioned algorithm contract.

---

# 17. Top-level builder

```r
build_factor_reference(
  fits,
  dataset_id = names(fits),
  representations = NULL,
  comparison_manifest,
  source_feature_sets = NULL,
  source_feature_sd = NULL,
  reference_features = NULL,
  min_shared_features = 2L,
  min_similarity = NULL,
  rank_max = NULL,
  reciprocal = NULL,
  min_cluster_density = 2/3,
  min_datasets = 2L,
  projection_manifest = NULL,
  keep_source_objects = FALSE
)
```

## Defaults by mode

`comparison_manifest$mode="generic"`:

- metric `congruence`;
- `min_similarity` is required from the caller;
- no rank restriction unless requested.

The package must not silently substitute STARCAT's `0.5` threshold for generic signed/EBMF programs.

`comparison_manifest$mode="starcat2025"`:

- Pearson;
- `min_similarity=0.5` if omitted;
- reciprocal top-rank filter;
- `rank_max=length(unique(dataset_id))` if omitted;
- source feature sets and source feature SD are required.

## Output

Returns a validated `factor_reference`.

---

# 18. Projection readiness

A consensus catalog can be scientifically useful even when no universal query preprocessing contract has been declared.

Therefore:

```text
projection_manifest = NULL
```

is allowed during reference construction.

Such an object has:

```text
projection_ready = FALSE
```

and `project_reference()` refuses it.

A projection-ready reference requires:

1. a manifest describing the expected query measurement basis; and
2. a single common activity-support contract across all selected reference programs, because schema-v1 native projection initializes all fixed programs jointly under one activity prior family.

References containing heterogeneous activity support remain valid for discovery/inspection but are not projection-ready in schema v1.

Example projection manifest:

```r
projection_manifest = list(
  feature_scale = "shifted_log_reference_library_normalized_measurements",
  scope = "external_reference",
  transformation = "...",
  compatibility = "exact_declared_contract"
)
```

If source datasets do not share a single projection transform, the caller must define the external reference transform explicitly rather than letting the package guess.

---

# 19. Native EBMF reference projection

```r
project_reference(
  reference,
  query,
  sample_side,
  feature_scale,
  preprocessing,
  programs = NULL,
  ebnm_fn = NULL,
  feature_match = c("exact", "intersection"),
  gene_mapping = NULL,
  min_coverage = 0.5,
  maxiter = 200L,
  S = NULL,
  var_type = 0L,
  normalize_activity = c("none", "l1"),
  dense_budget_bytes = 256 * 1024^2
)
```

## Implementation

Reuse the current native projection flow:

1. validate and orient query;
2. match/reference-map features;
3. subtract zero reference offset (schema v1);
4. initialize query activities;
5. `flash_factors_init()` with reference effects;
6. `flash_factors_fix()` on the feature side;
7. `flash_backfit()`;
8. verify reference feature programs did not change;
9. report raw query activity and reconstruction diagnostics.

## Priors

A consensus reference does not automatically inherit a fitted activity prior.

Core behavior:

- require all selected reference programs to share one compatible activity-support class in schema v1;
- if `ebnm_fn` is supplied, use it after validating its support when recognizable;
- otherwise use a projection prior family explicitly recorded in `projection_manifest`, if present;
- otherwise fail with an actionable error.

Per-program activity prior families are deferred because the current shared projection path initializes the fixed factors jointly. A later extension may implement sequential factor initialization if it is shown to preserve the intended joint model.

Do **not** silently select one source factor's fitted prior or pool fitted priors.

Query-side hyperparameters may re-estimate under the supplied prior family; this is recorded in metadata.

## L1 normalized activity

`normalize_activity="l1"` is allowed only when projected activity support is non-negative.

Raw activity remains primary and is always retained.

For zero-sum rows, normalized activity remains all-zero and a `zero_sum_activity` flag is returned.

---

# 20. Projection return object

Return class:

```r
class(x) <- c("reference_projection", "list")
```

Required fields:

```r
list(
  activity,
  activity_normalized,
  zero_sum_activity,
  representation,
  reconstruction,
  coverage,
  missing_reference_features,
  excluded_query_features,
  warnings,
  native_iteration_diagnostics,
  reference_id,
  analysis_metadata,
  provenance
)
```

`representation$effects` must equal the selected fixed reference program matrix exactly within the same tolerance used by the current projection implementation.

Reconstruction metrics are explicitly labeled **in-sample conditional fit diagnostics**, not independent prediction.

---

# 21. Residual-program discovery: feasibility gate first

Do not implement a public residual-program function until the tests in `06_FLASHIER_FEASIBILITY_PROBES.md` pass.

If they pass, expose:

```r
discover_residual_programs(
  reference,
  query,
  sample_side,
  feature_scale,
  preprocessing,
  max_residual = 10L,
  reference_ebnm_fn = NULL,
  residual_ebnm_fn = NULL,
  backfit = TRUE,
  ...
)
```

Model:

\[
X = A_{ref} B_{ref}^{T} + A_{res} B_{res}^{T} + E
\]

Rules:

- reference feature programs stay fixed;
- reference activities may update on query observations;
- residual factors are added with native `flash_greedy()` machinery;
- backfitting must preserve fixed reference feature programs;
- no public claim of biological novelty is made.

Output wording uses **residual program**, not novel program.

---

# 22. Known-versus-residual contribution

Do not report additive reference/residual factor PVE.

Report:

```text
known_fit_r2
full_fit_r2
delta_fit_r2_residual
known_sse
full_sse
delta_sse_residual
```

where the known model contains only fixed reference programs and the full model additionally contains residual programs.

The incremental metric is conditional on this model order and is labeled accordingly.

---

# 23. Residual-program evidence table

For every residual factor return:

```text
factor
best_reference_program
similarity
shared_features
native_pve_or_factor_metric
stability_if_available
status
```

Recommended descriptive statuses:

```text
reference_like_residual
distinct_residual
unstable_residual
ambiguous
```

No status is a biological annotation.

Independent replication or resampling is required before describing a residual factor as a new biological program.

---

# 24. Stability integration

Existing `factor_stability()` remains the within-dataset resampling engine.

Reference construction may accept optional source stability summaries and add them to `members_table`.

Consensus diagnostics report, without initially weighting by them:

```text
member_stability_median
member_stability_min
member_recovery_rate_mean
```

Stability-weighted consensus is experimental and benchmark-only until demonstrated to improve held-out transfer.

---

# 25. Diagnostics and visualization

Minimum diagnostics required before release:

1. pairwise similarity heatmap;
2. accepted reference graph;
3. consensus-by-dataset member matrix;
4. consensus versus member feature profiles;
5. support versus within-cluster similarity;
6. reference projection coverage;
7. raw activity and optional L1 activity comparison;
8. known-only versus known+residual reconstruction, if residual discovery is enabled.

Plot functions must never mutate the reference or projection object.

---

# 26. Provenance and credits

Every new public function must receive a row in `inst/provenance.csv` and corresponding documentation in:

- `CREDITS.md`;
- `vignettes/credits.Rmd`;
- roxygen author/details blocks;
- function-credit tables where relevant.

STARCAT should be credited for conceptual inspiration behind:

- cross-dataset consensus program catalog construction;
- reciprocal similarity graph logic;
- fixed-reference program projection concept.

`flashier` and `ebnm` authors remain credited for native inference and prior machinery.

No STARCAT source code needs to be copied for the planned implementation.

---

# 27. Recommended source-file layout

Prefer a small number of coherent files:

```text
R/
  reference-build.R        # profiles, similarity, graph, clustering, consensus
  reference-class.R        # factor_reference validation/accessors/print/summary
  reference-projection.R   # project_reference and shared native projection glue
  reference-residual.R     # conditional residual-program extension
  reference-plots.R        # only if needed
```

Avoid creating one R file per tiny helper.

Tests:

```text
tests/testthat/
  test-reference-profiles.R
  test-reference-similarity.R
  test-reference-graph.R
  test-reference-consensus.R
  test-reference-class.R
  test-reference-projection.R
  test-reference-residual.R     # conditional
```

---

# 28. Backward compatibility

Existing exported functions remain supported.

The new implementation should reuse shared internals where safe, but must preserve current regression outputs for:

- `factor_replicability()`;
- `project_factor_activity()`;
- `project_factor_features()`;
- `factor_holdout()`.

Do not deprecate these functions in the first reference release.

---

# 29. Core definition of done

Core reference/projection release is complete when:

1. multiple independent fits can be converted into a deterministic `factor_reference`;
2. comparison normalization never mutates the source inferential basis;
3. all-pairs similarities, graph exclusions and merge decisions are auditable;
4. no consensus cluster contains multiple factors from the same dataset;
5. final consensus programs are on an explicit common feature universe;
6. reference programs serialize and validate after reload;
7. projection keeps reference feature programs fixed;
8. feature mapping/coverage and preprocessing contracts are explicit;
9. legacy projection and alignment tests still pass;
10. STARCAT graph logic is validated in an external parity benchmark;
11. simulation and leave-one-dataset-out validation are documented;
12. provenance/credits and package checks pass.

Residual-program discovery has its own conditional definition of done and does not block the core consensus-reference release.
